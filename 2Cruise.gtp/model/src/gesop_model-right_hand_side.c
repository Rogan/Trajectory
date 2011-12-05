/* -------------------------------------------------------------------------------
 * Gesop_Model.Right_Hand_Side
 * -------------------------------------------------------------------------------
 *  Copyright
 *                        University of Stuttgart
 *                   Institute of Space Systems (IRS)
 *                          Pfaffenwaldring 31
 *                       70569 Stuttgart, Germany
 *
 *                        FAX: (+49)-711-685-63596
 *				     E-Mail: roeser@irs.uni-stuttgart.de
 * ---------------------------------------------------------------------------
 * Copyright (c) 2010, University of Stuttgart, IRS
 * All rights reserved.
 * ---------------------------------------------------------------------------
 * THE UNIVERSITY OF STUTTGART, IRS AND ASTOS GmbH DISCLAIM ALL WARRANTIES 
 * WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF 
 * MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THEY BE LIABLE FOR ANY 
 * SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER 
 * RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF 
 * CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN 
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 * ---------------------------------------------------------------------------
 * PROJECT:   BW-1
 * AUTHORS:   D. Fischer, C. Moellman, R. Shimmin
 * ---------------------------------------------------------------------------
 */

#include <stdio.h>
#include "gesop_model.h"
#include "nrutil.h"
#include "matrices.h"
#include "vectors.h"
#include "kepler.h"
#include "problem.h"
#include "mee2cart.h"
#include "disturbances.h"
#include "constants.h"

void __cdecl Right_Hand_Side(
	const int				*phase,		// Current phase number
	const int				*dimx,		// Dimension of state vector X
	const int				*dimu,		// Dimension of control vector U
	const int				*dimip,		// Dimension of integer parameter vector
	const int				*dimrp,		// Dimension of real parameter vector
	const Phase_Info_Type	*fazinf,	// Record with additional phase info
	const double			*t,			// Time t of evaluation
	const double			x[],		// State vector X at time t
	const double			u[],		// Control vector U at time t
	const int				ipar[],		// Integer parameter vector of phase
	const double			rpar[],		// Real parameter vector of phase
	double          		*dx,		// State derivative vector at time t
	int             		*error		// Error flag
	)
{
	/* computes the state derivative of the dynamic system (ODE)
	 * at the specified evaluation time of the given phase as a
	 * function of the phase specific real parameter vector and
	 * the state- and control vector at the evaluation time. */


	// Local variables
	double	p,f,g,h,k,L,time, mass; // renames for the state variables
	double	pdot, fdot, gdot, hdot, kdot, Ldot, mdot;  // time derivatives of state variables

	double	*mee_coords;
	double	*r_vec, *v_vec;
	double	*cart_inert_coords;
	double	**Q_transp;			//contains the rotation matrix for the rotating radial frame in which 
								//disturbing accelerations are expressed (x_rot = Q_transp * x_cart_inert)

	int		center, target;		//these are used for all calls to compute_ephemeris()
	double	julian_date;
	double	thrust;				// [N]
	double	ce;					// [m/s] exhaust velocity 
	double	power;				// [W]
	double	solar_power;		// [W]
	double	thruster_power;		// [W]
	double	payload_power;		// [W]
	double	comms_power;		// [W]
	double	thrust_mag;			// [-] Normalised thrust magnitude

	// thirdbody perturbations
	double	*pertearth;			// cartesian perturbing acceleration due to Earth (inertial)
	double	*pertsun;			// cartesian perturbing acceleration due to Sun (inertial)
	double	*pertmoon;			// cartesian perturbing acceleration due to Moon
	double	*pertjupiter;		// cartesian perturbing acceleration due to jupiter
	double	*pertvenus;			// cartesian perturbing acceleration due to venus
	double	*pertmars;			// cartesian perturbing acceleration due to mars
	
	double	*thirdbody_acc_cart, *thirdbody_acc_rotrad; // third body perturbing accelerations in cartesian and in rotating reference frame

	// gravity harmonic perturbations
	double	*pertharm_cart;		// inertial cartesian perturbing acceleration due to gravity harmonic
	double	*pertharm_rotrad;	// ... in rotating reference frame

	// solar pressure perturbations
	double	*pertsolp_cart;		//inertial cartesian perturbing acceleration due to solar pressure
	double	*pertsolp_rotrad;	// ... in rotating reference frame

	double	**inert2fixed;		//contains the rotation matrix for the transformation from inertial
								//to fixed coordinates (x_fixed = inert2fixed * x_cart_inert)

	double	ephemerides[6];
	double	*body2sun;
	double	*craft2sun;

	// auxiliary values for right hand sides
	double	omega,s2;

	// accelerations
	double	acc_r, acc_theta, acc_h;
	double	thrust_acc_r, thrust_acc_theta, thrust_acc_h;

	double	g_mue;		// gravitional constant depending on phase
	double	deltaL;
	const double LIBR_CENTER = 0.0;

	int		i, j;		//counter variables

#ifdef DEBUG_MODE
printf("Right Hand Side\n");_flushall();
#endif

	// Initializations

	*error = 0;

	Q_transp			= dmatrix(1,3,1,3);
	inert2fixed			= dmatrix(1,3,1,3);
	mee_coords			= dvector(1,6);
	cart_inert_coords	= dvector(1,6);
	r_vec				= dvector(1,3);
	v_vec				= dvector(1,3);
	
	pertearth			= dvector(1,3);
	pertsun				= dvector(1,3);
	pertmoon			= dvector(1,3);
	pertjupiter			= dvector(1,3);
	pertvenus			= dvector(1,3);
	pertmars			= dvector(1,3);

	thirdbody_acc_cart	= dvector(1,3);
	thirdbody_acc_rotrad = dvector(1,3);

	pertharm_cart		= dvector(1,3);
	pertharm_rotrad		= dvector(1,3);

	pertsolp_cart		= dvector(1,3);
	pertsolp_rotrad		= dvector(1,3);	

	craft2sun			= dvector(1,3);
	body2sun			= dvector(1,3);
	
	for (i=1; i<=3; i++){
		for (j=1; j<=3; j++){
			Q_transp[i][j]		= 0.0;
			inert2fixed[i][j]	= 0.0;
		}

		pertearth[i]			= 0.0;
		pertsun[i]				= 0.0;
		pertmoon[i]				= 0.0;
		pertjupiter[i]			= 0.0;
		pertvenus[i]			= 0.0;
		pertmars[i]				= 0.0;

		thirdbody_acc_cart[i]	= 0.0;
		thirdbody_acc_rotrad[i] = 0.0;

		pertharm_cart[i]		= 0.0;
		pertharm_rotrad[i]		= 0.0;
		
		pertsolp_cart[i]		= 0.0;
		pertsolp_rotrad[i]		= 0.0;

		craft2sun[i]			= 0.0;
		body2sun[i]				= 0.0;

		r_vec[i]				= 0.0;
		v_vec[i]				= 0.0;
	}

	acc_r		= 0.0;
	acc_theta	= 0.0;
	acc_h		= 0.0;

	power = 0.0;

	// Retrieve state
	for (i=1; i<=6; i++) {
		mee_coords[i] = x[i-1];
		cart_inert_coords[i] = 0.0;	
	}

	// Rename states
	p = mee_coords[1];
	f = mee_coords[2];
	g = mee_coords[3];
	h = mee_coords[4];
	k = mee_coords[5];
	L = mee_coords[6];
	time = x[6];
	mass = x[7];

	// Retrieve optimised parameters
	deltaL = rpar[2];				// Get optimised phase length
	thrust_mag = u[3];				// Get optimised normalised thrust magnitude

//	if (x[8] < 100)
//		thrust_mag = 0.0;
//	else 
		thrust_mag = u[3];				// Get optimised normalised thrust magnitude

#ifdef DEBUG_MODE
	printf("---PHASE %i---\n", *phase);_flushall();
#endif

	if (*phase < 3)
		center = EARTH;
	else
		center = MOON;
	g_mue = gmue[center];

	switch(*phase)
	{
		case 1:
			thrust = thrust_mag*thrust_arc;
			ce = ce_arc;
			thruster_power = thrust_mag*power_arc;
			payload_power = 0.0;
//			break;
		case 2:
			thrust = thrust_mag*thrust_ppt;
			ce = ce_ppt;
			thruster_power = thrust_mag*power_ppt;
			payload_power = 0.0;
			break;
		case 3:
			thrust = thrust_mag*thrust_arc;
			ce = ce_arc;
			thruster_power = thrust_mag*power_arc;
			payload_power = 0.0;
			break;
		case 4:
			thrust = thrust_mag*thrust_ppt;
			ce = ce_ppt;
			thruster_power = thrust_mag*power_ppt;
			payload_power = 0.0;
			break;
		case 5:
			thrust = 0.0;
			ce = 1.0;
			thruster_power = 0.0;
			payload_power = 1000.0;
			break;
		default:
			printf("Unrecognised phase number in RHS function\n");
			break;
	} // switch
	comms_power = 0.0;

	// auxilliary values
	omega	= 1.0 + f*cos(L) + g*sin(L);
	s2		= 1.0 + h*h + k*k;

	mee2cart(mee_coords, g_mue, cart_inert_coords);		// obtain inertial cartesian coords
	ir_itheta_ih_transp(cart_inert_coords, Q_transp);	// obtain rotating radial coordinate frame Q_transp as [ir itheta ih]'
	julian_date = EPOCH + rpar[0] + time;				// = start time + current mission time

	// split cartesian coords into radius and velocity-vector
	for (i=1; i<=3; i++){						
		r_vec[i] = cart_inert_coords[i];
		v_vec[i] = cart_inert_coords[i+3];
	}


	/******** evaluate accelerations *******/
	// --- thrust ---
	thrust_acc_r		= thrust/mass * u[0];		
	thrust_acc_theta	= thrust/mass * u[1];
	thrust_acc_h		= thrust/mass * u[2];

	switch (PERT_MODE)
	{
	case 3:
		// --- minor 3rd body perturbations ---
		target = JUPITER;
		third_body_pert_cart(julian_date, center, target, r_vec, pertjupiter);
				
		target = VENUS;
		third_body_pert_cart(julian_date, center, target, r_vec, pertvenus);
				
		target = MARS;
		third_body_pert_cart(julian_date, center, target, r_vec, pertmars);

	case 2:
		// --- perturbations due to oblateness ---
		compute_rotation(julian_date, center, inert2fixed);

		if (center == EARTH) {	
			AccelHarmonic (r_vec, inert2fixed, GM_JGM3, R_JGM3, CS_JGM3, 20, 20, pertharm_cart);
		} else {
			AccelHarmonic (r_vec, inert2fixed, GM_LP100J, R_LP100J, CS_LP100J, 20, 20, pertharm_cart);
		}
		dmat_times_vec(Q_transp, pertharm_cart, pertharm_rotrad, 3, 3); //convert to rotating radial frame

		// --- perturbations due to solar pressure ---
		solar_rad_pressure_cart(julian_date, mass, cart_inert_coords, center, pertsolp_cart);
		dmat_times_vec(Q_transp, pertsolp_cart, pertsolp_rotrad, 3, 3);


	case 1:
		// --- major 3rd body perturbations ---
		if (center == EARTH) {	
			target = MOON;
			third_body_pert_cart(julian_date, center, target, r_vec, pertmoon);
		} else {	
			target = EARTH;
			third_body_pert_cart(julian_date, center, target, r_vec, pertearth);
		}
		target = SUN;
		third_body_pert_cart(julian_date, center, target, cart_inert_coords, pertsun);

		// sum 3rd body perturbations
		for (i=1; i<=3; i++){
			thirdbody_acc_cart[i] = pertearth[i] + pertsun[i] + pertmoon[i] + pertjupiter[i] + pertvenus[i] + pertmars[i];
		}
		dmat_times_vec(Q_transp, thirdbody_acc_cart, thirdbody_acc_rotrad, 3, 3); //convert to rotating radial frame
	} 
	
	//total
	acc_r		= thrust_acc_r		+ thirdbody_acc_rotrad[1] + pertharm_rotrad[1] + pertsolp_rotrad[1];
	acc_theta	= thrust_acc_theta	+ thirdbody_acc_rotrad[2] + pertharm_rotrad[2] + pertsolp_rotrad[2];
	acc_h		= thrust_acc_h		+ thirdbody_acc_rotrad[3] + pertharm_rotrad[3] + pertsolp_rotrad[3];

	
	/******** power calculations *******/

	target = SUN;
	compute_ephemeris(julian_date, target, center, ephemerides);
	for(i=1;i<=3;i++) {
		body2sun[i] = ephemerides[i-1]*1000.0;
		craft2sun[i] = body2sun[i] - cart_inert_coords[i];
	}
	
//	if(center == EARTH)
//		solar_power = SOLARINTENSITY*aref*area_efficiency*conversion_efficiency*shadowfunc(body2sun,cart_inert_coords,R_JGM3)*sin(sunangle(craft2sun,u));
//	else
//		solar_power = SOLARINTENSITY*aref*area_efficiency*conversion_efficiency*shadowfunc(body2sun,cart_inert_coords,R_LP100J)*sin(sunangle(craft2sun,u));
//	power = solar_power - thruster_power - payload_power - comms_power;


	// right hand sides
	// -----------------------------------------------------------
	pdot = sqrt(p/g_mue)* 2.0*p/omega*acc_theta; 																				// [km/s]
	fdot = sqrt(p/g_mue) *(  sin(L)*acc_r + ((omega+1.0)*cos(L) + f)*acc_theta/omega - (h*sin(L) - k*cos(L))*g/omega*acc_h );	// [1/s]
	gdot = sqrt(p/g_mue) *( -cos(L)*acc_r + ((omega+1.0)*sin(L) + g)*acc_theta/omega + (h*sin(L) - k*cos(L))*f/omega*acc_h );	// [1/s]
	hdot = sqrt(p/g_mue) * s2*cos(L)/(2.0*omega)*acc_h;  																		// [1/s]
	kdot = sqrt(p/g_mue) * s2*sin(L)/(2.0*omega)*acc_h;  																		// [1/s]
	Ldot = sqrt(p/g_mue) *( h*sin(L) - k*cos(L) )*acc_h/omega + sqrt(p*g_mue)*(omega/p)*(omega/p);								// [1/s]
	mdot = -thrust / ce;  																										// [kg/s]	

if(Ldot < 0.0)
	printf("\n\nLDOT NEGATIVE!!!!!\n\n");_flushall();
if(p < 0.0)
	printf("\n\nP NEGATIVE!!!!!\n\n");_flushall();

	// normalized L shall be the independent variable, so dx/dLn = dx/dt * dt/dL * dL/dLn
	dx[0] = pdot * deltaL / Ldot;				// [km]
	dx[1] = fdot * deltaL / Ldot;				// [-]
	dx[2] = gdot * deltaL / Ldot;
	dx[3] = hdot * deltaL / Ldot;
	dx[4] = kdot * deltaL / Ldot;
	dx[5] = deltaL;								// L_phase_f - L_phase_i
	dx[6] = deltaL / Ldot / 86400.0;			// [days]
	dx[7] = mdot * deltaL / Ldot;				// [kg]
	dx[8] = 0.0;//power * deltaL / Ldot;				// [J]


	//clean-up
	free_dvector(mee_coords,1,6);
	free_dvector(cart_inert_coords,1,6);

	free_dvector(pertearth,1,3);
	free_dvector(pertsun,1,3);
	free_dvector(pertmoon,1,3);
	free_dvector(pertjupiter,1,3);
	free_dvector(pertvenus,1,3);
	free_dvector(pertmars,1,3);

	free_dvector(thirdbody_acc_cart,1,3);
	free_dvector(thirdbody_acc_rotrad,1,3);

	free_dvector(pertharm_cart,1,3);
	free_dvector(pertharm_rotrad,1,3);

	free_dvector(pertsolp_cart,1,3);
	free_dvector(pertsolp_rotrad,1,3);

	free_dmatrix(Q_transp,1,3,1,3);
	free_dmatrix(inert2fixed,1,3,1,3);

	free_dvector(body2sun,1,3);
	free_dvector(craft2sun,1,3);
}