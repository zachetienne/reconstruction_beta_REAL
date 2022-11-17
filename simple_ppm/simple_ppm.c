#include "stdio.h"
#include "math.h"

#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )


#define MINUS2 0
#define MINUS1 1
#define PLUS_0  2
#define PLUS_1  3
#define PLUS_2  4


//Eq. 60 in JOURNAL OF COMPUTATIONAL PHYSICS 123, 1-14 (1996)
//   [note the factor of 2 missing in the |a_{j+1} - a_{j}| term].
//   Recall that dU = U_{i} - U_{i-1}.
static double slope_limit(const double dU, const double dUp1) {
  // Set SLOPE_LIMITER_COEFF = 2.0 for MC, 1 for minmod
#define SLOPE_LIMITER_COEFF 2.0

  if(dU*dUp1 > 0.0) {
    //delta_m_U=0.5 * [ (u_(i+1)-u_i) + (u_i-u_(i-1)) ] = (u_(i+1) - u_(i-1))/2  <-- first derivative, second-order; this should happen most of the time (smooth flows)
    const double delta_m_U = 0.5*(dU + dUp1);
    // EXPLANATION OF BELOW LINE OF CODE.
    // In short, sign_delta_a_j = sign(delta_m_U) = (0.0 < delta_m_U) - (delta_m_U < 0.0).
    //    If delta_m_U>0, then (0.0 < delta_m_U)==1, and (delta_m_U < 0.0)==0, so sign_delta_a_j=+1
    //    If delta_m_U<0, then (0.0 < delta_m_U)==0, and (delta_m_U < 0.0)==1, so sign_delta_a_j=-1
    //    If delta_m_U==0,then (0.0 < delta_m_U)==0, and (delta_m_U < 0.0)==0, so sign_delta_a_j=0
    const int sign_delta_m_U = (0.0 < delta_m_U) - (delta_m_U < 0.0);
    //Decide whether to use 2nd order derivative or first-order derivative, limiting slope.
    return sign_delta_m_U*MIN(fabs(delta_m_U),MIN(SLOPE_LIMITER_COEFF*fabs(dUp1),SLOPE_LIMITER_COEFF*fabs(dU)));
  }
  return 0.0;
}

static void compute_UrUl_onevar(const double U[7], double *Ur, double *Ul) {
  const double slope_limited_dU_m1 = slope_limit(U[MINUS1] - U[MINUS2], U[PLUS_0] - U[MINUS1]);
  const double slope_limited_dU_p0 = slope_limit(U[PLUS_0] - U[MINUS1], U[PLUS_1] - U[PLUS_0]);
  const double slope_limited_dU_p1 = slope_limit(U[PLUS_1] - U[PLUS_0], U[PLUS_2] - U[PLUS_1]);

  *Ur = 0.5*(U[PLUS_1] + U[PLUS_0]) + (1.0/6.0)*(slope_limited_dU_p0 - slope_limited_dU_p1);
  *Ul = 0.5*(U[PLUS_0] + U[MINUS1]) + (1.0/6.0)*(slope_limited_dU_m1 - slope_limited_dU_p0);
}


// Compute ftilde, which is used for flattening left and right face values
// DEPENDENCIES: P(MINUS2,MINUS1,PLUS_1,PLUS_2) and v^m(MINUS1,PLUS_1), where m=flux_dirn={1,2,3}={x,y,z}.
#define OMEGA1   0.75
#define OMEGA2  10.0
#define EPSILON2 0.33
static double shock_detection__ftilde(const double P[7], const double v_flux_dirn[7]) {
  double dP1 = P[PLUS_1] - P[MINUS1];
  double dP2 = P[PLUS_2] - P[MINUS2];

  // MODIFICATION TO STANDARD PPM:
  // Cure roundoff error issues when dP1==0 or dP2==0 to 15 or more significant digits.
  const double avg1=0.5*(P[PLUS_1] + P[MINUS1]);
  const double avg2=0.5*(P[PLUS_2] + P[MINUS2]);
  if(fabs(dP1)/avg1<1e-15) dP1=0.0; /* If this is triggered, there is NO shock */
  if(fabs(dP2)/avg2<1e-15) dP2=0.0; /* If this is triggered alone, there may be a shock. Otherwise if triggered with above, NO shock. */

  double dP1_over_dP2=1.0;
  if (dP2 != 0.0) dP1_over_dP2 = dP1/dP2;

  const double q1 = (dP1_over_dP2-OMEGA1)*OMEGA2;
  const double q2 = fabs(dP1)/MIN(P[PLUS_1], P[MINUS1]);

  // w==0 -> NOT inside a shock
  double w=0.0;

  // w==1 -> inside a shock
  if (q2 > EPSILON2 && q2*( (v_flux_dirn[MINUS1]) - (v_flux_dirn[PLUS_1]) ) > 0.0) w = 1.0;

  return MIN(1.0, w*MAX(0.0,q1));
}

// standard Colella-Woodward parameters:
//    K0 = 0.1d0, eta1 = 20.0, eta2 = 0.05, epsilon = 0.01d0
#define K0      0.1
#define ETA1   20.0
#define ETA2    0.05
#define EPSILON 0.01
static void steepen_rho(const double rho[7],const double P[7], const double Gamma_eff,
                        double *rho_br_ppm, double *rho_bl_ppm) {

  // Next compute centered differences d RHOB and d^2 RHOB
  const double d1rho_b     = 0.5*(rho[PLUS_1] - rho[MINUS1]);
  const double d2rho_b_m1  = rho[PLUS_0] - 2.0*rho[MINUS1] + rho[MINUS2];
  const double d2rho_b_p1  = rho[PLUS_2] - 2.0*rho[PLUS_1] + rho[PLUS_0];

  // Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
  const double contact_discontinuity_check = Gamma_eff*K0*fabs(rho[PLUS_1]-rho[MINUS1])*
    MIN(P[PLUS_1],P[MINUS1])
    -fabs(P[PLUS_1]-P[MINUS1])*MIN(rho[PLUS_1],rho[MINUS1]);
  const double second_deriv_check = -d2rho_b_p1*d2rho_b_m1;
  const double relative_change_check = fabs(2.0*d1rho_b) - EPSILON*MIN(rho[PLUS_1],rho[MINUS1]);

  if(contact_discontinuity_check >= 0.0 && second_deriv_check >= 0.0
     && relative_change_check >= 0.0) {

    const double slope_limited_drho_m1 = slope_limit(rho[MINUS1] - rho[MINUS2], rho[PLUS_0] - rho[MINUS1]);
    const double slope_limited_drho_p1 = slope_limit(rho[PLUS_1] - rho[PLUS_0], rho[PLUS_2] - rho[PLUS_1]);

    double eta_tilde=0.0;
    if (fabs(d1rho_b) > 0.0) {
      eta_tilde = -(1.0/6.0)*(d2rho_b_p1-d2rho_b_m1)/(2.0*d1rho_b);
    }
    const double eta = MAX(0.0,MIN(ETA1*(eta_tilde - ETA2),1.0));
    // Next compute Urp1 and Ul for RHOB, using the MC prescription:
    // Ur_p1 = U_p1   - 0.5*slope_lim_dU_p1
    const double rho_br_mc_p1 = rho[PLUS_1] - 0.5*slope_limited_drho_p1;
    // Ul = U_m1 + 0.5*slope_lim_dU_m1
    // Based on this line of code, Ur[index] = a_j - \delta_m a_j / 2. (cf. Eq. 65 in Marti & Muller's "PPM Method for 1D Relativistic Hydro." paper)
    //    So: Ur[indexp1] = a_{j+1} - \delta_m a_{j+1} / 2. This is why we have rho_br_mc[indexp1]
    const double rho_bl_mc    = rho[MINUS1] + 0.5*slope_limited_drho_m1;

    *rho_bl_ppm = (*rho_bl_ppm)*(1.0-eta) + rho_bl_mc*eta;
    *rho_br_ppm = (*rho_br_ppm)*(1.0-eta) + rho_br_mc_p1*eta;

  }
}


static void flatten_r_and_l(const double U, const double ftilde, double *Ur, double *Ul) {
  *Ur = U*ftilde + (*Ur)*(1.0-ftilde);
  *Ul = U*ftilde + (*Ul)*(1.0-ftilde);
}


static inline void monotonize_r_and_l(const double U, double *Ur, double *Ul) {
  const double dU = (*Ur) - (*Ul);
  const double mU = 0.5*((*Ur)+(*Ul));

  if ( ((*Ur)-U)*(U-(*Ul)) <= 0.0) {
    (*Ur) = U;
    (*Ul) = U;
    return;
  }
  if ( dU*(U-mU) > (1.0/6.0)*(dU*dU)) {
    (*Ul) = 3.0*U - 2.0*(*Ur);
    return;
  }
  if ( dU*(U-mU) < -(1.0/6.0)*(dU*dU)) {
    (*Ur) = 3.0*U - 2.0*(*Ul);
    return;
  }
}


// Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
void simple_ppm_1D(const double rho[7], const double P[7],
                   const double vx[7], const double vy[7], const double vz[7],
                   const double other_vars[8][7], const int num_other_vars,
                   const double v_flux_dirn[7], const double Gamma_eff) {

  // Interpolate primitives to faces with a slope limiter.
  double rhor, rhol;  compute_UrUl_onevar(rho, &rhor, &rhol);
  double Pr  , Pl;    compute_UrUl_onevar(P,   &Pr,   &Pl);
  double vxr , vxl;   compute_UrUl_onevar(vx,  &vxr,  &vxl);
  double vyr , vyl;   compute_UrUl_onevar(vy,  &vyr,  &vyl);
  double vzr , vzl;   compute_UrUl_onevar(vz,  &vzr,  &vzl);
  double other_varsr[8], other_varsl[8];
  for(int var=0;var<num_other_vars;var++) {
    compute_UrUl_onevar(other_vars[var], &other_varsr[var], &other_varsl[var]);
  }

  // Steepen rho
  steepen_rho(rho, P, Gamma_eff, &rhor, &rhol);

  // Flatten and monotonize all variables
  // First detect shocks / steep gradients:
  {
    const double ftilde = shock_detection__ftilde(P, v_flux_dirn);
    flatten_r_and_l(   rho[PLUS_0],ftilde, &rhor,&rhol);
    monotonize_r_and_l(rho[PLUS_0],        &rhor,&rhol);
    flatten_r_and_l(   P[PLUS_0],  ftilde, &Pr,&Pl);
    monotonize_r_and_l(P[PLUS_0],          &Pr,&Pl);
    flatten_r_and_l(   vx[PLUS_0], ftilde, &vxr,&vxl);
    monotonize_r_and_l(vx[PLUS_0],         &vxr,&vxl);
    flatten_r_and_l(   vy[PLUS_0], ftilde, &vyr,&vyl);
    monotonize_r_and_l(vy[PLUS_0],         &vyr,&vyl);
    flatten_r_and_l(   vz[PLUS_0], ftilde, &vzr,&vzl);
    monotonize_r_and_l(vz[PLUS_0],         &vzr,&vzl);
    for(int var=0;var<num_other_vars;var++) {
      flatten_r_and_l(   other_vars[var][PLUS_0], ftilde, &other_varsr[var],&other_varsl[var]);
      monotonize_r_and_l(other_vars[var][PLUS_0],         &other_varsr[var],&other_varsl[var]);
    }
  }
}
