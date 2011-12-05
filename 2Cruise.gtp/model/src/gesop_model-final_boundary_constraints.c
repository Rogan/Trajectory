/* -------------------------------------------------------------------------------
 * Gesop_Model.Final_Boundary_Constraints
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
#include "matrices.h"
#include "nrutil.h"
#include "kepler-mee.h"
#include "mee2cart.h"
#include "problem.h"
#include "constants.h"

void __cdecl Final_Boundary_Constraints(
	const int				*phase,		// Current phase number
 	const int				*dimx,		// Dimension of state vector X1
    const int				*dimu,		// Dimension of control vector U1
	const int				*dimip,		// Dimension of integer parameter vector
	const int				*dimrp,		// Dimension of real parameter vector
	const Phase_Info_Type	*fazinf,	// Record with additional phase info
	const double			*t,			// Time t of evaluation
	const double			x[],		// State vector X at time t
	const double			u[],		// Control vector U at time t
	const int				ipar[],		// Integer parameter vector of phase
	const double			rpar[],		// Real parameter vector of phase
 	const int				*dimbc,		// Dimension of boundary constraint vector
	const Decision_Type		evalc[],	// Which constraints to evaluate
	double          		*bcon,		// Vector of boundary constraints
	int             		*error		// Error flag
	)
{
	/* computes the boundary constraints associated with the final
	 * time of the specified phase (tf of the given phase) as a
	 * function of the phase final time, the state- and control
	 * vector at time tf and the real parameter vector of the phase. */

	double	*mee_coords_ECI;			// contains the first six elements of the state, ie the equ. elements
	double	*kep_coords_ECI;			// contains the first six elements of the state in earth centred keplerian elements
	double	*cart_inert_coords_ECI;		// contains radius and vel. vector of satellite in earth centred inertial coords
	double	*cart_inert_coords_LCI;		// contains radius and vel. vector of satellite in moon centred inertial coords
	double	r_ECI, v_ECI;				// magnitude of the radius and vel. vector of satellite in earth centred inertial coords
	double	r_LCI, v_LCI;				// magnitude of the radius and vel. vector of satellite in moon centred inertial coords
	double	r_per;						// periapsis of osculating orbit
	double	ephemerides_moon[6];		// radius and vel. vector of the Moon in earth centred inertial coords
	double	epsilon_LCI, epsilon_ECI;	// orbital energy
	double	julian_date;
	int		target, center;
	int		i;

#ifdef DEBUG_MODE
printf("Final Boundary Constraints\n");_flushall();
#endif
	
	*error = 0;	

	// Initialise storage
	mee_coords_ECI			= dvector(1,6);				
	kep_coords_ECI			= dvector(1,6);	/*			
	cart_inert_coords_ECI	= dvector(1,6);
	cart_inert_coords_LCI	= dvector(1,6);*/
			
	for (i=1; i<=6; i++) {
		mee_coords_ECI[i]			= x[i-1];
		kep_coords_ECI[i]			= 0.0;/*
		cart_inert_coords_ECI[i]	= 0.0;
		cart_inert_coords_LCI[i]	= 0.0;*/
	}

	// Calculate parameters
	mee2kep(mee_coords_ECI, kep_coords_ECI);
	//mee2cart(mee_coords_ECI, gmue[EARTH], cart_inert_coords_ECI); 
	//r_ECI = dvec_abs(cart_inert_coords_ECI,3);
	//v_ECI = dvec_abs(&cart_inert_coords_ECI[3],3);
	//epsilon_ECI = v_ECI*v_ECI/2 - gmue[EARTH]/r_ECI;

	//julian_date = EPOCH + rpar[0] + x[6];	// = start time + current mission time [days] = current time
	//target = MOON;
	//center = EARTH;
	//compute_ephemeris(julian_date, target, center, ephemerides_moon);

	//for (i=1; i<=6; i++) {
	//	cart_inert_coords_LCI[i] = cart_inert_coords_ECI[i] - ephemerides_moon[i-1]*1.0e3;
	//}

	//r_per = kep_coords_ECI[1]*(1-kep_coords_ECI[2]);	
	//r_LCI = dvec_abs(cart_inert_coords_LCI,3);
	//v_LCI = dvec_abs(&cart_inert_coords_LCI[3],3);
	//epsilon_LCI = v_LCI*v_LCI/2 - gmue[MOON]/r_LCI;

	// calculate only values required for given phase

	// divide all parameters by required value to scale constraints to Order(1)

	// Continuity constraints acquired from Propagate phase --------------------------------------------

	/* a */
	if (evalc[0] == yes) {
		bcon[0] = kep_coords_ECI[1]/268406334 - 1.0;
	}

	/* e */
	if (evalc[1] == yes) {
		bcon[1] = kep_coords_ECI[2]/0.0828 - 1.0;
	}

	/* i */
	if (evalc[2] == yes) {
		bcon[2] = kep_coords_ECI[3]*R2D/20.382 - 1.0;
	}

	/* o */
	if (evalc[3] == yes) {
		bcon[3] = kep_coords_ECI[4]*R2D/208.745 - 1.0;
	}

	/* O */
	if (evalc[4] == yes) {
		bcon[4] = kep_coords_ECI[5]*R2D/329.162 - 1.0;
	}

	/* v */
	if (evalc[5] == yes) {
		bcon[5] = kep_coords_ECI[6]*R2D/52.667 - 1.0;
	}

	/* t */
	if (evalc[6] == yes) {
		bcon[6] = x[6]/1652.9534 - 1.0;
	}

	/* m */
	if (evalc[7] == yes) {
		bcon[7] = x[7]/181.98 - 1.0;
	}


	// Free memory
	free_dvector(mee_coords_ECI,1,6);
	free_dvector(kep_coords_ECI,1,6);/*
	free_dvector(cart_inert_coords_ECI,1,6);
	free_dvector(cart_inert_coords_LCI,1,6);*/
}