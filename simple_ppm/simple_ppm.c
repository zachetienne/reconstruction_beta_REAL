#include "stdio.h"
#include "math.h"

#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )


//Eq. 60 in JOURNAL OF COMPUTATIONAL PHYSICS 123, 1-14 (1996)
//   [note the factor of 2 missing in the |a_{j+1} - a_{j}| term].
//   Recall that dU = U_{i} - U_{i-1}.
static double slope_limit(const double dU,const double dUp1) {
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

  *Ur = 0.5*(U[PLUS_1] + U[PLUS_0]) + (1.0/6.0)*(slope_lim_dU_p0 - slope_lim_dU_p1);
  *Ul = 0.5*(U[PLUS_0] + U[MINUS1]) + (1.0/6.0)*(slope_lim_dU_m1 - slope_lim_dU_p0);
}

//typedef struct __ppm_
#define MINUS2 0
#define MINUS1 1
#define PLUS_0  2
#define PLUS_1  3
#define PLUS_2  4
void simple_ppm_1D(const double rho[7], const double P[7],
                   const double vx[7], const double vy[7], const double vz[7],
                   const double other_vars[8][7], const int num_other_vars) {
  
  
}
