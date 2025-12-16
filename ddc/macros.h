#ifndef MACROS_H
#define MACROS_H
/**
 * @brief this file defines governing equations of double-component systems like double-diffusive convection.
 * 
 * A system of double components has form:
 * Momentum equation (variable U, P):           dt(U) + div(grad(U)) = - grad(P) + p1*lap(U) + p2(p3*T-p4*S)*ey
 * First component's equation (variable T):     dt(T) + div(grad(T)) =             p5*lap(T)
 * Second component's equation (variable S):    dt(S) + div(grad(S)) =             p6*lap(S) + p7*lap(T)
 * where U, = Ubase + u, is total velocity.
 * 
 */

// Example 1: Stationary-wall bounded double-diffusive convection
// velocity is normalized by \kappa_T/H
// #define P1 Pr 
// #define P2 Pr*Ra
// #define P3 1.0
// #define P4 Rrho
// #define P5 1.0
// #define P6 1.0/Le

// Example 2: Stationary/Moving-wall bounded double-diffusive convection
// velocity is normalized by free-fall velocity (= boundary velocity U_b for case of Ri=1)
// here: Rho ~ Salinity/Temperature
#define P1 sqrt(Pr/Ra) 
#define P2 1.0
#define P3 1.0
#define P4 Rrho
#define P5 1.0/sqrt(Pr*Ra)
#define P6 1.0/(Le*sqrt(Pr*Ra))

// Example 3: Moving-wall bounded double-diffusive convection
// velocity is normalized by boundary velocity U_b
// #define P1 1.0/Rey
// #define P2 Ri/(Rrho-1)
// #define P3 1.0
// #define P4 Rrho
// #define P5 1.0/(Rey*Pr)
// #define P6 1.0/(Le*Rey*Pr)

// Example 4: Binary fluid convection [Mercader2013JFM]
// #define P1 Pr 
// #define P2 Pr*Ra
// #define P3 (1.0+Rsep)
// #define P4 Rsep
// #define P5 1.0
// #define P6 1.0/Le
// #define P7 1.0

/* Example 5: Stratified plane Couette flow [Langham2019JFM] */
// #define P1 1.0/Rey
// #define P2 (-1.0)
// #define P3 Ri
// #define P5 (1.0/(Rey*Pr))

// Example 6: Inclined layer convection [Zheng JFM 2024] / RBC
// #define P1 sqrt(Pr/Ra) 
// #define P2 1.0
// #define P3 1.0
// #define P5 1.0/sqrt(Pr*Ra)

// Example 7: Plane Couette flow 
// #define P1 1.0/Rey

// Example 8: Rotating RBC
// #define P1 sqrt(Pr/Ra) 
// #define P2 1.0
// #define P3 1.0
// #define P5 1.0/sqrt(Pr*Ra)
// #define P8 (1.0/Ek*sqrt(Pr/Ra))  // Coriolis term coefficient

// Example 9: Rotating salt-finger convection 
// here: Rho ~ Temperature/Salinity
// #define P1 sqrt(Pr/Ra) 
// #define P2 1.0
// #define P3 1.0
// #define P4 1.0/Rrho
// #define P5 1.0/sqrt(Pr*Ra)
// #define P6 1.0/(Le*sqrt(Pr*Ra))
// #define P8 1.0/Ek*sqrt(Pr/Ra)  // Coriolis term coefficient

// #define FREESLIP
// #define SAVESTATS
// #define FREEZEvelocity
#endif