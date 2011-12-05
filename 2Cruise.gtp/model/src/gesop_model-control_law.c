/* -------------------------------------------------------------------------------
 * Gesop_Model.Control_Law
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
#include "mee2cart.h"
#include "constants.h"
#include "disturbances.h"

/* --- Initial guesses for optimisable parameters --- */
/* These parameters remain constant throughout each phase */
 
const double start_time = 0.0; 	// [days]; days since January 1st, 2014, 0h (=midnight)

const double deltaL[1] = {
	4627.0, 
}; 	// phase duration (independent parameter is true anomaly rather than time)
// Needed to normalise the timescale within each phase


/* --- Initial guesses for initial state --- */

const double Ini_p = 27588412.28;			// [m] a(1-e²)
const double Ini_f = 0.229611658;			// e*cos(w+O)
const double Ini_g = -0.004196298;			// e*sin(w+O)
const double Ini_h = 0.230358303;			// tan(i/2)*cos(O)
const double Ini_k = -0.06110852;			// tan(i/2)*sin(O)
const double Ini_L = 6.33048373;			// w+O+v
const double Ini_t = 594.2725;				// [days]
const double Ini_m = 197.034844421;			// [kg]
const double Ini_b = 2.8*3.6e6;				// [J] = [kWh]*[J/kWh]

/* optimisable parameter bounds */
const double start_time_lb = 0.0;		// [days]
const double start_time_ub = 600.0;		// [days]

/* control bounds */
const double LB_Control = - 1.05;		// Allow some slack for optimisation; control constraint will prevent vector exceeding 1.
const double UB_Control =   1.05;

/* state bounds */
const double p_lb =  0.0;		// [m] For parabolic orbit. Should not be less than 6.4e6 pv. still in elliptical orbit.
const double p_ub =  2.0e9;//4.0e8;		// For circular orbit at lunar radius. Will not be greater unless hyperbolic.

const double f_lb = -1.0;		// [-]
const double f_ub =  1.0;		// ATTENTION : For e < 1.0 (f and g could be even bigger, as with e -> infinite, f,g -> infinite)

const double g_lb = -1.0;		// [-]
const double g_ub =  1.0;		// ATTENTION : For e < 1.0

const double h_lb = -1.0;		// [-]
const double h_ub =  1.0;		// ATTENTION : For i < 90.0 deg (h and k could be even bigger, as with i -> 180 deg, h,k -> infinite)

const double k_lb = -1.0;		// [-]
const double k_ub =  1.0;		// ATTENTION : For i < 90.0 deg

const double L_lb = -20.0;		// [rad]
const double L_ub = 20000.0;

const double t_lb = 0.0;		// [days]
const double t_ub = 2206.0;		// 2458864.5 (jpleph limit) - 2456658.5 (EPOCH)

const double m_lb = 50.0;		// [kg]
const double m_ub = 250.0;

const double b_lb = 0.0;		// [J]
const double b_ub = 2.8*3.6e6;	// [2.8 kWh]







void __cdecl Initialize_Control_Law(
	double	*ini_time,		    // Initial time (t0 of phase number one)
	double	*ini_time_lb,	    // Lower bound of initial time
	double	*ini_time_ub,	    // Upper bound of initial time
	int		*error)				// Error flag
{
/* Initialize control law and return start time and its lower/upper bounds. */

#ifdef DEBUG_MODE
printf("Initialize Control Law\n");_flushall();
#endif

	*error = 0;

	// Since "time" is Ln
	*ini_time = 0.0;
	*ini_time_lb = 0.0;
	*ini_time_ub = 0.0;
}







void __cdecl Initialize_Control_Law_Phase(
	const int		*phase,			// Current phase number
	const int		*dimx,			// Dimension of state vector X
	const int		*dimip,			// Dimension of integer parameter vector
	const int		*dimrp,		    // Dimension of real parameter vector
 	const double	*ini_time,	    // Initial time of phase
	int				*ipar,			// Integer parameter vector of phase
	int				*ipar_lb,		// Integer parameter vector of phase
	int				*ipar_ub,		// Integer parameter vector of phase
	double			*rpar,		    // Real parameter vector of phase
	double			*rpar_lb,		// Lower bounds for real parameter vector
 	double			*rpar_ub,		// Upper bounds for real parameter vector
	double			*ini_state,	    // Initial state vector
	double			*ini_state_lb,	// Lower bounds for initial state vector
	double			*ini_state_ub,	// Upper bounds for initial state vector
	double			*fin_time,		// Final time (tf of phase number one)
	double			*fin_time_lb,	// Lower bound for final time
	double			*fin_time_ub,	// Upper bound for final time
	int				*error)			// Error flag
{
/* Initialize the control law for the specified phase and
 * return the parameters, inital states and final time for the
 * phase as well as their lower and upper bounds.
 *
 * For all phases other than phase one, the real parameters and
 * initial states as computed by the "Phase_Connect" routine are
 * an input to the routine, so only the bounds have to be set sensibly
 * (this only applies for items for which Connected flag is set to "true"). */

	int		i;

#ifdef DEBUG_MODE
printf("Initialize Control Law Phase\n");_flushall();
#endif

	*error = 0;

	/*	integer parameters */
	for (i=0; i<*dimip; i++) {
		ipar[i] = 0;
		ipar_lb[i] = 0;
		ipar_ub[i] = 0;
	}

	/* real parameters */
	rpar[0]		= start_time;
	rpar_lb[0]	= start_time_lb;
	rpar_ub[0]	= start_time_ub;

	rpar[2]		= deltaL[*phase-1];
	rpar_lb[2]	= rpar[2]*0.1;
	rpar_ub[2]	= rpar[2]*10.0;


	/* Initial guess for state */
	if (*phase == 1) {	
		ini_state[0] = Ini_p;
		ini_state[1] = Ini_f;
		ini_state[2] = Ini_g;
		ini_state[3] = Ini_h;
		ini_state[4] = Ini_k;
		ini_state[5] = Ini_L;  
		ini_state[6] = Ini_t;
		ini_state[7] = Ini_m;
		ini_state[8] = Ini_b;
	}

	ini_state_lb[0] = p_lb;
	ini_state_ub[0] = p_ub;

	ini_state_lb[1] = f_lb;
	ini_state_ub[1] = f_ub;

	ini_state_lb[2] = g_lb;
	ini_state_ub[2] = g_ub;

	ini_state_lb[3] = h_lb;
	ini_state_ub[3] = h_ub;

	ini_state_lb[4] = k_lb;
	ini_state_ub[4] = k_ub;

	ini_state_lb[5] = L_lb;
	ini_state_ub[5] = L_ub;

	ini_state_lb[6] = t_lb;
	ini_state_ub[6] = t_ub;
	
	ini_state_lb[7] = m_lb;
	ini_state_ub[7] = m_ub;

	ini_state_lb[8] = b_lb;
	ini_state_ub[8] = b_ub;

	/*	final time (normalized L)*/
	*fin_time		= (double) *phase;
	*fin_time_lb	= *fin_time;
	*fin_time_ub	= *fin_time;
}







void __cdecl Finalize_Control_Law(int *error)
{
#ifdef DEBUG_MODE
printf("Finalize Control Law\n");_flushall();
#endif

	/*	clean any internal state of the package */
	*error = 0;
}

void __cdecl Get_State_Bounds(
	const int				*phase,		// Current phase number
	const int				*dimx,	    // Dimension of state vector X
	const int				*dimip,	    // Dimension of integer parameter vector
	const int				*dimrp,	    // Dimension of real parameter vector
	const Phase_Info_Type	*fazinf,	// Record with additional phase info
	const double			*t,		    // Time t of evaluation
	const double			x[],		// State vector X at time t
	const int				ipar[],	    // Integer parameter vector of phase
	const double			rpar[],	    // Real parameter vector of phase
	double					*x_lb,		// Lower bounds for the state vector X
	double					*x_ub,		// Upper bounds for the state vector X
	int						*error)	    // Error flag
{
	/*	return the upper and lower bounds for the states
	 *	(output vectors will have the same range as State). */

#ifdef DEBUG_MODE
printf("Get State Bounds\n");_flushall();
#endif

	*error = 0;

	x_lb[0] = p_lb;
	x_ub[0] = p_ub;

	x_lb[1] = f_lb;
	x_ub[1] = f_ub;

	x_lb[2] = g_lb;
	x_ub[2] = g_ub;

    x_lb[3] = h_lb;
	x_ub[3] = h_ub;

	x_lb[4] = k_lb;
	x_ub[4] = k_ub;

	x_lb[5] = L_lb;
	x_ub[5] = L_ub;

	x_lb[6] = t_lb;
	x_ub[6] = t_ub;

	x_lb[7] = m_lb;
	x_ub[7] = m_ub;

	x_lb[8] = b_lb;
	x_ub[8] = b_ub;
}








void __cdecl Get_Controls(
	const int				*phase,		// Current phase number
	const int				*dimx,	    // Dimension of state vector X
	const int				*dimu,	    // Dimension of control vector U
	const int				*dimip,		// Dimension of integer parameter vector
	const int				*dimrp,	    // Dimension of real parameter vector
	const Phase_Info_Type	*fazinf,	// Record with additional phase info
	const double			*t,		    // Time t of evaluation
	const double			x[],	    // State vector X at time t
	const int				ipar[],	    // Integer parameter vector of phase
	const double			rpar[],     // Real parameter vector of phase
	double					*u,		    // Control vector U at time t
	double					*u_lb,		// Lower bounds for the control vector U
	double					*u_ub,		// Upper bounds for the control vector U
	int						*error)	    // Error flag
{
	/*	compute the control and its bounds as a function
	 *	of time, states and parameters. */

	double	*mee_coords_LCI;
	double	*cart_inert_coords_LCI;
	double	*v_vec_norm;
	double	*v_vec_norm_rotrad;
	double	**Q_transp;
	double	v_norm;
	int		i,j;

#ifdef DEBUG_MODE
printf("Get Controls\n");_flushall();
#endif

	*error = 0;

	mee_coords_LCI			= dvector(1,6);
	cart_inert_coords_LCI	= dvector(1,6);
	v_vec_norm				= dvector(1,3);
	v_vec_norm_rotrad		= dvector(1,3);
	Q_transp				= dmatrix(1,3,1,3);

	for (i=1; i<=3; i++){
		for (j=1; j<=3; j++){
			Q_transp[i][j]			= 0.0;
		}
		v_vec_norm[i]				= 0.0;
		v_vec_norm_rotrad[i]		= 0.0;
	}

	for (i=1; i<=6; i++) {
		mee_coords_LCI[i]			= x[i-1];
		cart_inert_coords_LCI[i]	= 0.0;
	}


	// Calculate velocity vector
	mee2cart(mee_coords_LCI, gmue[MOON], cart_inert_coords_LCI);

	v_norm = dvec_abs(&cart_inert_coords_LCI[3],3);
	for (i=1; i<=3; i++) {
		v_vec_norm[i] = cart_inert_coords_LCI[i+3] / v_norm;
	}

	ir_itheta_ih_transp(cart_inert_coords_LCI, Q_transp); //obtain Q_transp as [ir itheta ih]'
	dmat_times_vec(Q_transp, v_vec_norm, v_vec_norm_rotrad, 3, 3); //convert to rotating radial frame


	/* Initial control profile thrust tangentially to radius. */
	// Future work: set this initial guess separately for each phase
	u[0] =  v_vec_norm_rotrad[1];//0.0;
	u[1] =  v_vec_norm_rotrad[2];//1.0;
	u[2] =  v_vec_norm_rotrad[3];//0.0;
	
	u[3] =	1.0;//(cos(x[5])+1)/2;


	/* Bounds on control parameters */
	u_lb[0]	= LB_Control;
	u_ub[0]	= UB_Control;
	
	u_lb[1]	= LB_Control;
	u_ub[1]	= UB_Control;
		
	u_lb[2]	= LB_Control;
	u_ub[2]	= UB_Control;

	u_lb[3] = 0.0;
	u_ub[3] = 1.0;

	free_dvector(mee_coords_LCI,1,6);
	free_dvector(cart_inert_coords_LCI,1,6);
	free_dvector(v_vec_norm,1,3);
	free_dvector(v_vec_norm_rotrad,1,3);
	free_dmatrix(Q_transp,1,3,1,3);
}









void __cdecl Phase_Is_Over(
	const int				*phase,		// Current phase number
	const int				*dimx,	    // Dimension of state vector X
	const int				*dimip,	    // Dimension of integer parameter vector
	const int				*dimrp,	    // Dimension of real parameter vector
	const Phase_Info_Type	*fazinf,	// Record with additional phase info
	const double			*t,		    // Time t of evaluation
	const double			x[],		// State vector X at time t
	const int				ipar[],	    // Integer parameter vector of phase
	const double			rpar[],     // Real parameter vector of phase
	int						*isover,	// Flag for end of phase
	int						*error)	    // Error flag
{
	/*	check if the final conditions for a phase are met.
	 *	There must be no more than one transition from 'false' to 'true'
	 *	along the trajectory given by the control law.
	 *	The condition is only evaluated in the time interval between the
	 *	minimum and the nominal end time of the phase. */
	
#ifdef DEBUG_MODE
printf("Phase is Over\n");_flushall();
#endif

	*error = 0;
}