//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file blast.cpp
//  \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//         cylindrical, and spherical coordinates.
//

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#define DIRECT_INPUT 0
#define WOOD_MODEL 1
#define COHEN_MODEL 2
#define JOHNSTONE_MODEL 3
#define MESQUITA_MODEL 4

// MACHINE UNIT VALUES:
// FIXED VALUES
// distance: 1 AU
// velocity: 52.483 km/s
// mass density: 2.8E13 kg/AU^3 (or 5 m_p/cm^3)
// mass: 2.8E13 kg
// DERIVED VALUES
// time: 33.0 days
// pressure: 2.30E11 N/m^2
// energy density: 2.3E11 N/m^2
// magnetic field: 5.38E-9 T
// momentum: 262.4 (m_p cm^-3)(km s^-1)
// temperature: 111110 K


int input_mode;

Real threshold, rin, rout, rin_mag;

Real vx_ISM, vy_ISM, vz_ISM, d_ISM, p_ISM;
Real b0, angle;

Real d_wind, mom_rad_wind, e_wind, b_rad_wind;

Real x1_0, x2_0, x3_0;

bool split_monopole, plasma_sheet;

Real hps_height, hps_factor;

int RefinementCondition(MeshBlock *pmb);

void WindSource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void ISMBoundary_ix1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  EnrollUserExplicitSourceFunction(WindSource);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, ISMBoundary_ix1);
  // if (adaptive) {
  //   EnrollUserRefinementCondition(RefinementCondition);
  //   threshold = pin->GetReal("problem","thr");
  // }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  x1_0 = 0.0;
  x2_0 = 0.0;
  x3_0 = 0.0;
  vx_ISM = pin->GetOrAddReal("problem", "vx_ISM", 0.0);
  vy_ISM = pin->GetOrAddReal("problem", "vy_ISM", 0.0);
  vz_ISM = pin->GetOrAddReal("problem", "vz_ISM", 0.0);
  d_ISM = pin->GetOrAddReal("problem", "d_ISM", 0.0);
  p_ISM = pin->GetOrAddReal("problem", "p_ISM", 0.0);
  rout = pin->GetReal("problem", "radius");
  rin  = pin->GetOrAddReal("problem", "radius_inner", 1.0);
  rin_mag = pin->GetOrAddReal("problem", "radius_inner_mag", 1.0);
  b_rad_wind = pin->GetOrAddReal("problem", "b_rad_wind", 0.0);
  e_wind = pin->GetOrAddReal("problem", "e_wind", 0.0);

  input_mode = pin->GetOrAddInteger("problem", "input_mode", DIRECT_INPUT);

  if(input_mode == DIRECT_INPUT){
    mom_rad_wind = pin->GetOrAddReal("problem", "mom_rad_wind", 0.0);
    d_wind = pin->GetOrAddReal("problem", "d_wind", 0.0);
  }
  else{
    Real v_wind = pin->GetOrAddReal("problem","v_wind",7.62); //in machine units
    Real M_dot_solar = 129.0; //solar mass loss rate, equal to 2E-14 M_solar yr^-1 (from Johnstone et al. 2015)
    if(input_mode == WOOD_MODEL){
      Real R_star = pin->GetOrAddReal("problem","R_star_over_solar",1.0); //stellar radius in solar radii
      Real L_x_star = pin->GetOrAddReal("problem","L_x_star_over_solar",1.0); //stellar total X-ray luminosity, in solar units
      d_wind = M_dot_solar/(v_wind*4.0*PI*rin*rin) * std::pow(R_star,-0.52)*std::pow(L_x_star,1.26);
      mom_rad_wind = d_wind*v_wind;
      printf("wood density calculated\n");
      printf("%e",d_wind);
    }
    else if(input_mode == COHEN_MODEL){
      Real G_grav = 6.674E-11; //in SI units
      Real mu_0 = 1.257E-6; //in SI units
      Real B_0 = pin->GetOrAddReal("problem","B_0",0.0005); //weak dipolar field component, input in tesla
      Real d_Alfven = pin->GetOrAddReal("problem","d_Alfven",6.6);
      Real R_star = pin->GetOrAddReal("problem","R_star",1.0)*6.96E8; //stellar radius, input in solar radii -> meters
      Real M_star = pin->GetOrAddReal("problem","M_star",1.0)*1.99E30; //stellar mass, input in solar mass -> kilograms
      d_wind = B_0*B_0*std::pow(R_star,2.5) / ((v_wind*52.483*1000.0)*mu_0*std::pow(rin*1.496E11,2.0)*std::pow(d_Alfven,4.0)*std::sqrt(2.0*G_grav*M_star));
      d_wind *= 1.196E20; //convert back to machine units from SI
      printf("cohen density calculated\n");
      printf("%e\n",d_wind);
      mom_rad_wind = d_wind*v_wind;
    }
    else if(input_mode == JOHNSTONE_MODEL){
      Real G_grav = 6.674E-11; //in SI units
      Real avg_particle_mass = 0.6*1.673E-27; //kg
      Real gamma = 5.0/3.0;
      Real k_B = 1.381E-23;

      Real M_solar = 1.989E30 / 2.8E13; //machine units
      Real one_year = 365.0 / 33.0; //machine units

      Real c_s_over_v_e = pin->GetOrAddReal("problem","c_s_over_v_e",0.329);
      Real R_star = pin->GetOrAddReal("problem","R_star",1.0); //stellar radius, input in solar radii
      Real M_star = pin->GetOrAddReal("problem","M_star",1.0)*1.99E30; //stellar mass, input in solar mass -> kilograms
      Real rho_0 = pin->GetOrAddReal("problem","rho_0",1.0E7*avg_particle_mass/0.6); //g cm^{-3}

      Real T_0 = 2.0*G_grav*avg_particle_mass*c_s_over_v_e*c_s_over_v_e*M_star / (gamma*k_B*R_star*6.96E8); //kelvin
      T_0 /= 1E6; //megakelvin
      printf("johnstone base temp calculated\n");
      printf("%e MK\n",T_0);

      Real log_parameter = 3.42+10.20*std::log10(T_0)-1.94*std::pow(std::log10(T_0),2.0)+3.79*std::pow(std::log10(T_0),3.0);
      Real M_dot = rho_0*R_star*R_star*std::pow(10.0,log_parameter)*M_solar/one_year; //machine units
      d_wind = M_dot / (4.0*PI*rin*rin*v_wind);
      printf("johnstone density calculated\n");
      printf("%e\n",d_wind);
      mom_rad_wind = d_wind*v_wind;
    }
    else if(input_mode == MESQUITA_MODEL){
      Real L_x_star = pin->GetOrAddReal("problem","L_x_star",5.7E26); //erg s^{-1}
      Real R_star = pin->GetOrAddReal("problem","R_star",0.437); //in solar radii
      rin = 300.0*R_star*0.00465;
      Real lambda_N = std::pow(10.0,-24.5); //erg cm^3 s^{-1}
      Real proton_mass = 1.673E-27*1000.0; //g
      Real k_B = 1.381E-23; //SI

      L_x_star /= 10.0; //assume 10% of X-ray flux is from open field corona

      Real rho_0 = proton_mass * std::sqrt(L_x_star/(2.0*PI*std::pow(R_star*6.96E10,3.0)*lambda_N)); //g cm^{-3}
      if(rho_0 < 1.0E-14){ //low beta
        v_wind = 1.59E-5*std::pow(rho_0,-0.55);//km s^{-1}
        d_wind = 1.44E15*std::pow(rho_0,2.74);//gm cm^{-3}
      } else { //high beta
        v_wind = 1.64E7*std::pow(rho_0,0.31);//km s^{-1}
        d_wind = 8.35E-3*std::pow(rho_0,1.49);//gm cm^{-3}
      }
      Real T_pl = 7.8E14*std::pow(rho_0,0.61);//K
      Real n_SI = d_wind/proton_mass*1.0E6; //m^-3
      Real thermal_energy = 2.0*n_SI*k_B*T_pl/(2.0/3.0); //N m^-2
      thermal_energy /= 2.3E11; //machine units
      v_wind /= 52.483;//machine units
      d_wind /= 5.0*proton_mass; //machine units
      mom_rad_wind = d_wind*v_wind; //machine units
      e_wind = 0.5*d_wind*v_wind*v_wind + thermal_energy;//machine units
      printf("mesquita model\n");
      printf("rho_0: %e\n",rho_0);
      printf("d_wind: %e\n",d_wind);
      printf("v_wind: %e\n",v_wind);
      printf("mom_rad_wind: %e\n",mom_rad_wind);
      printf("T_pl: %e\n",T_pl);
      printf("e_wind: %e\n",e_wind);
    }
    // mom_rad_wind = 8.75 * (M_dot/9.13E-14);
    // b_rad_wind = pin->GetOrAddReal("problem", "b_rad_wind", 0.0);
    // d_wind = 1.1 * (M_dot/9.13E-14);
    // e_wind = pin->GetOrAddReal("problem", "e_wind", 0.0);
  }
  
  split_monopole = pin->GetOrAddBoolean("problem", "split_monopole", false);
  plasma_sheet = pin->GetOrAddBoolean("problem", "plasma_sheet", false);
  hps_height = pin->GetOrAddReal("problem", "hps_height", 0.01);
  hps_factor = pin->GetOrAddReal("problem", "hps_factor", 4.0);
  Real pa  = p_ISM;
  Real da  = d_ISM;
  if (MAGNETIC_FIELDS_ENABLED) {
    b0 = pin->GetReal("problem", "b0");
    angle = (PI/180.0)*pin->GetReal("problem", "angle");
  }
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  // setup uniform ambient medium
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        
        Real den = da;

        phydro->u(IDN,k,j,i) = d_ISM;

        // MOMENTUM-INITIALIZED ISM
        phydro->u(IM1,k,j,i) = d_ISM*vx_ISM;
        phydro->u(IM2,k,j,i) = d_ISM*vy_ISM;
        phydro->u(IM3,k,j,i) = d_ISM*vz_ISM;

        //BOUNDARY-ONLY ISM MOMENTUM
        // phydro->u(IM1,k,j,i) = 0.0;
        // phydro->u(IM2,k,j,i) = 0.0;
        // phydro->u(IM3,k,j,i) = 0.0;

        if (NON_BAROTROPIC_EOS) {
          Real pres = pa;
          phydro->u(IEN,k,j,i) = pres/gm1; //thermal energy
          phydro->u(IEN,k,j,i) += 0.5*d_ISM*(vx_ISM*vx_ISM+vy_ISM*vy_ISM+vz_ISM*vz_ISM); //bulk kinetic energy
          if (RELATIVISTIC_DYNAMICS)  // this should only ever be SR with this file
            phydro->u(IEN,k,j,i) += den;
        }
      }
    }
  }
              // Real falloff = rin*rin/(rad*rad);
              // Real curr_b_rad_wind = b_rad_wind * falloff;
              // if(y<0.0) curr_b_rad_wind *= -1.0; //enforce split monopole
              // cons(IB1,k,j,i) = curr_b_rad_wind * x / rad;
              // cons(IB2,k,j,i) = curr_b_rad_wind * y / rad;
              // cons(IB3,k,j,i) = curr_b_rad_wind * z / rad;

  // initialize solar wind magnetic field
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          Real rad;
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            Real x = pcoord->x1v(i);
            Real y = pcoord->x2v(j);
            Real z = pcoord->x3v(k);
            rad = std::sqrt(SQR(x - x1_0) + SQR(y - x2_0) + SQR(z - x3_0));
            Real falloff = rin*rin/(rad*rad);
            Real curr_b_rad_wind = b_rad_wind * falloff;
            if(split_monopole && y<0.0) curr_b_rad_wind *= -1.0; //enforce split monopole
            // pfield->b.x1f(k,j,i) = b0;
            // pfield->b.x1f(k,j,i) = b0 * std::cos(angle);
            pfield->b.x1f(k,j,i) = curr_b_rad_wind * (x / rad);
          }
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real rad;
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            Real x = pcoord->x1v(i);
            Real y = pcoord->x2v(j);
            Real z = pcoord->x3v(k);
            rad = std::sqrt(SQR(x - x1_0) + SQR(y - x2_0) + SQR(z - x3_0));
            Real falloff = rin*rin/(rad*rad);
            Real curr_b_rad_wind = b_rad_wind * falloff;
            if(split_monopole && y<0.0) curr_b_rad_wind *= -1.0; //enforce split monopole
            pfield->b.x2f(k,j,i) = curr_b_rad_wind * (y / rad);
            // pfield->b.x2f(k,j,i) = 0.0;
            // pfield->b.x2f(k,j,i) = b0 * std::sin(angle);
          }
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real rad;
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            Real x = pcoord->x1v(i);
            Real y = pcoord->x2v(j);
            Real z = pcoord->x3v(k);
            rad = std::sqrt(SQR(x - x1_0) + SQR(y - x2_0) + SQR(z - x3_0));
            Real falloff = rin*rin/(rad*rad);
            Real curr_b_rad_wind = b_rad_wind * falloff;
            if(split_monopole && y<0.0) curr_b_rad_wind *= -1.0; //enforce split monopole
            pfield->b.x3f(k,j,i) = curr_b_rad_wind * (z / rad);
            // pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real rad;
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            Real x = pcoord->x1v(i);
            Real y = pcoord->x2v(j);
            Real z = pcoord->x3v(k);
            rad = std::sqrt(SQR(x - x1_0) + SQR(y - x2_0) + SQR(z - x3_0));
            Real falloff = rin*rin/(rad*rad);
            Real curr_b_rad_wind = b_rad_wind * falloff;
            phydro->u(IEN,k,j,i) += 0.5*curr_b_rad_wind*curr_b_rad_wind;
          }
          // phydro->u(IEN,k,j,i) += 0.5*b0*b0;
        }
      }
    }
  }
}

// Boundary condition definition
void ISMBoundary_ix1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set primitive variables in inlet ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,il-i) = d_ISM;
          prim(IVX,k,j,il-i) = vx_ISM;
          prim(IVY,k,j,il-i) = vy_ISM;
          prim(IVZ,k,j,il-i) = vz_ISM;
          prim(IPR,k,j,il-i) = p_ISM;
      }
    }
  }

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
            // b.x1f(k,j,il-i) = b0;
            // b.x2f(k,j,il-i) = 0.0;
            // b.x3f(k,j,il-i) = 0.0;
            b.x1f(k,j,il-i) = b0*std::cos(angle);
            b.x2f(k,j,il-i) = b0*std::sin(angle);
            b.x3f(k,j,il-i) = 0.0;
        }
      }
    }
  }
}

// User-defined source term
void WindSource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons) {
  // get coordinate location of the center, convert to Cartesian
  Real x0, y0, z0;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ParameterInput" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real rad, x, y, z;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          x = pmb->pcoord->x1v(i);
          y = pmb->pcoord->x2v(j);
          z = pmb->pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          x = pmb->pcoord->x1v(i)*std::cos(pmb->pcoord->x2v(j));
          y = pmb->pcoord->x1v(i)*std::sin(pmb->pcoord->x2v(j));
          z = pmb->pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
          x = pmb->pcoord->x1v(i)*std::sin(pmb->pcoord->x2v(j))*std::cos(pmb->pcoord->x3v(k));
          y = pmb->pcoord->x1v(i)*std::sin(pmb->pcoord->x2v(j))*std::sin(pmb->pcoord->x3v(k));
          z = pmb->pcoord->x1v(i)*std::cos(pmb->pcoord->x2v(j));
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        }
        if (rad < rout) {
          Real hps_height = 10.0;
          Real hps_factor = 4.0; //multiplicative factor by which density and momentum density increase at HPS
          if (rad < rin) {
            Real curr_mom_rad_wind = mom_rad_wind;
            //augments density and momentum (proportionally) at equator to account for heliospheric plasma sheet
            if(plasma_sheet && std::abs(z)<=hps_height){
              Real t = (hps_height - std::abs(z))/hps_height; //varies from 0 at outer edge to 1 at equator
              Real factor = (3.0*t*t - 2.0*t*t*t) * (hps_factor-1.0) + 1.0; //varies factor from 1 at outer edge to 4 at equator (cubic Hermite spline)
              curr_mom_rad_wind *= factor;
              cons(IDN,k,j,i) = d_wind*factor;
            } else {
              cons(IDN,k,j,i) = d_wind;
            }
            cons(IM1,k,j,i) = curr_mom_rad_wind * x / rad;
            cons(IM2,k,j,i) = curr_mom_rad_wind * y / rad;
            cons(IM3,k,j,i) = curr_mom_rad_wind * z / rad;
            cons(IEN,k,j,i) = e_wind;
            if (MAGNETIC_FIELDS_ENABLED && rad < rin_mag) {
              Real curr_b_rad_wind = b_rad_wind;
              if(split_monopole && y<0.0) curr_b_rad_wind *= -1.0; //enforce split monopole
              cons(IB1,k,j,i) = curr_b_rad_wind * (x / rad);
              cons(IB2,k,j,i) = curr_b_rad_wind * (y / rad);
              cons(IB3,k,j,i) = curr_b_rad_wind * (z / rad);
            }
          } else {   // ramp all variables down as 1/r^2 (???) but continuous with inner value
            Real falloff = rin*rin/(rad*rad);
            Real curr_mom_rad_wind = mom_rad_wind * falloff;
            //augments density and momentum (proportionally) at equator to account for heliospheric plasma sheet
            if(plasma_sheet && std::abs(z)<=hps_height){
              Real t = (hps_height - std::abs(z))/hps_height; //varies from 0 at outer edge to 1 at equator
              Real factor = (3.0*t*t - 2.0*t*t*t) * (hps_factor-1.0) + 1.0; //varies factor from 1 at outer edge to 4 at equator (cubic Hermite spline)
              curr_mom_rad_wind *= factor;
              cons(IDN,k,j,i) = d_wind * falloff * factor;
            } else {
              cons(IDN,k,j,i) = d_wind * falloff;
            }
            cons(IM1,k,j,i) = curr_mom_rad_wind * x / rad;
            cons(IM2,k,j,i) = curr_mom_rad_wind * y / rad;
            cons(IM3,k,j,i) = curr_mom_rad_wind * z / rad;
            cons(IEN,k,j,i) = e_wind * falloff;
            if (MAGNETIC_FIELDS_ENABLED && rad < rin_mag) {
              Real curr_b_rad_wind = b_rad_wind * falloff;
              if(split_monopole && y<0.0) curr_b_rad_wind *= -1.0; //enforce split monopole
              cons(IB1,k,j,i) = curr_b_rad_wind * x / rad;
              cons(IB2,k,j,i) = curr_b_rad_wind * y / rad;
              cons(IB3,k,j,i) = curr_b_rad_wind * z / rad;
            }
          }
        }
      }
    }
  }
  return;    
}
