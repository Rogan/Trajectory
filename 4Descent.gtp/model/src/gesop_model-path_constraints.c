/* -------------------------------------------------------------------------------
 * Gesop_Model.Path_Constraints
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
#include "kepler.h"
#include "vectors.h"
#include "nrutil.h"
#include "kepler-mee.h"
#include "mee2cart.h"
#include "problem.h"
#include "constants.h"

void __cdecl Path_Constraints(
	const int				*phase,		// Current phase number
	const int				*dimx,		// Dimension of state vector X1
    const int				*dimu,		// Dimension of control vector U1
	const int				*dimip,		// Dimension of integer parameter vector
	const int				*dimrp,		// Dimension of real parameter vector
	const Phase_Info_Type	*fazinf, 	// Record with additional phase info
	const double			*t,			// Time t of evaluation
	const double			x[],		// State vector X at time t
	const double			u[],		// Control vector U at time t
	const double			udot[],		// Control derivative vector at t
	const int				ipar[],		// Integer parameter vector of phase
	const double			rpar[],		// Real parameter vector of phase
 	const int				*dimpc,		// Dimension of path constraint vector
	const Decision_Type		evalc[],	// Which constraints to evaluate
	double          		*pcon,		// Vector of path constraints at t
	int             		*error		// Error flag
	)
{
	/* computes the path constraints at the specified evaluation
	 * time of the given phase as a function of the phase specific
	 * real parameter vector and the state, control and control
	 * derivative vector at the evaluation time. */

	double	*mee_coords_ECI;			// contains the first six elements of the state, ie the equ. elements
	double	*kep_coords_ECI;			// contains the first six elements of the state in earth centred keplerian elements
	double	*cart_inert_coords_ECI;		// contains radius and vel. vector of satellite in earth centred inertial coords
	double	*cart_inert_coords_LCI;		// contains radius and vel. vector of satellite in moon centred inertial coords
	double	r_ECI, v_ECI;				// contains the magnitude of the radius and vel. vector of satellite in earth centred inertial coords
	double	r_LCI, v_LCI;				// contains the magnitude of the radius and vel. vector of satellite in moon centred inertial coords
	double	ephemerides_moon[6];		// contains radius and vel. vector of the Moon in earth centred inertial coords
	double	julian_date;
	int		target, center;
	int		i;

#ifdef DEBUG_MODE
printf("Path Constraints\n");_flushall();
#endif

	*error			= 0;

	// Initialise storage
	mee_coords_ECI			= dvector(1,6);				
	kep_coords_ECI			= dvector(1,6);				
	cart_inert_coords_ECI	= dvector(1,6);
	cart_inert_coords_LCI	= dvector(1,6);
			
	for (i=1; i<=6; i++) {
		mee_coords_ECI[i]			= x[i-1];
		kep_coords_ECI[i]			= 0.0;
		cart_inert_coords_ECI[i]	= 0.0;
		cart_inert_coords_LCI[i]	= 0.0;
	}

	// Calculate parameters
	mee2kep(mee_coords_ECI, kep_coords_ECI);
	mee2cart(mee_coords_ECI, gmue[EARTH], cart_inert_coords_ECI); 

	r_ECI = dvec_abs(cart_inert_coords_ECI,3);
	v_ECI = dvec_abs(&cart_inert_coords_ECI[3],3);

	julian_date = EPOCH + rpar[0] + x[6];	// = start time + current mission time [days] = current time
	target = MOON;
	center = EARTH;
	compute_ephemeris(julian_date, target, center, ephemerides_moon);

	for (i=1; i<=6; i++) {
		cart_inert_coords_LCI[i] = cart_inert_coords_ECI[i] - ephemerides_moon[i-1]*1.0e3;
	}
			
	r_LCI = dvec_abs(cart_inert_coords_LCI,3);
	v_LCI = dvec_abs(&cart_inert_coords_LCI[3],3);	

	
	/* spacecraft thrust is a unit vector */
	if (evalc[0] == yes) {
		pcon[0] = 1.0 - sqrt(pow(u[0],2) + pow(u[1],2) + pow(u[2],2));
	}
	
	/* spacecraft should not go below 100km altitude */
	if (evalc[1] == yes) {
		pcon[1] = r_ECI*1.0e-3/(R_JGM3 + 100.0) - 1.0;
	}
	
	/* spacecraft should not escape Earth orbit */
	if (evalc[2] == yes) {
		pcon[2] = 0.9 - kep_coords_ECI[2];
	}

	/* spacecraft should not go within 50km altitude of the Moon */
	if (evalc[3] == yes) {
		pcon[3] = r_LCI*1.0e-3/(R_LP100J + 50.0e3) - 1.0;
	}

	/* spacecraft inclination should not approach 90° */
	if (evalc[4] == yes) {
		pcon[4] = 1.0 - kep_coords_ECI[3]*R2D/90;
	}

	/* spacecraft mass should be greater than zero(!) */
	if (evalc[5] == yes) {
		pcon[4] = x[7];
	}	

	// Free memory
	free_dvector(mee_coords_ECI,1,6);
	free_dvector(kep_coords_ECI,1,6);
	free_dvector(cart_inert_coords_ECI,1,6);
	free_dvector(cart_inert_coords_LCI,1,6);
}