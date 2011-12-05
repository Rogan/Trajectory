/* -------------------------------------------------------------------------------
 * Gesop_Model.Phase_Connect
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
#include "vectors.h"
#include "constants.h"
#include "kepler.h"
#include "mee2cart.h"
#include "problem.h"
#include "cart2mee.h"

void __cdecl Phase_Connect(
	const int				*phase,		// Current phase number
	const int				*dimx1,		// Dimension of state vector X1
    const int				*dimu1,		// Dimension of control vector U1
	const int				*dimip1,	// Dimension of integer parameter vector 1
	const int				*dimrp1,	// Dimension of real parameter vector 1
    const int				*dimx2,		// Dimension of state vector X2
    const int				*dimu2,		// Dimension of control vector U2
	const int				*dimip2,	// Dimension of integer parameter vector 2
	const int				*dimrp2,	// Dimension of real parameter vector 2
	const Phase_Info_Type	*fazinf,	// Record with addition. phase info
	const double			*t,			// Phase separation time t
	const double			x_1[],		// Current phase state vector at t
	const double			u_1[],		// Current phase control vector at t
	const int				ipar_1[],	// Current phase integer parameters
	const double			rpar_1[],	// Current phase real parameters
	const Decision_Type		xconec[],	// Connection flags for State_2
	const Decision_Type		uconec[],	// Connection flags for Control_2
	const Decision_Type		iconec[],	// Connection flags for Int_Param_2
	const Decision_Type		rconec[],	// Connection flags for Real_Param_2
	double					*x_2,		// Next phase state vector at t
	double					*u_2,		// Next phase control vector at t
	int						*ipar_2,	// Next phase int parameter vector
	double					*rpar_2,	// Next phase real parameter vector
	int						*error		// Error flag
	)
{
	/* computes the phase connect conditions for the state, control
	 * and real parameter vector between the given phase and the
	 * following phase (phase+1) at the phase separation time. */

	/* local variables */
	int		i;								// counter variables
	int		center, target;					// these are used for all calls to compute_ephemeris()
	double	julian_date;
	double	ephemerides_moon[6];
	double	*mee_coords_ECI;				// contains the first six elements of the state, ie the equ. elements
	double	*mee_coords_LCI;				// contains mod. eq. elements of satellite in LCI coords
	double	*kep_coords_LCI;				// contains kepl. elements of satellite in LCI coords
	double	*cart_inert_coords_ECI;			// contains the first six elements of the state in inertial cartesian coords.
	double	*cart_inert_coords_LCI;			// contains radius and vel. vector of satellite in LCI coords



#ifdef DEBUG_MODE
printf("Phase Connect\n");_flushall();
#endif

	*error = 0;

	
	/*---Control parameters stay the same in all phases---*/
	for (i=0; i<*dimu2; i++) {
 		u_2[i] = u_1[i];
 	}


	/*---Real optimisation parameters stay the same in all phases---*/
	for (i=0; i<*dimrp2; i++) {
 		rpar_2[i] = rpar_1[i];
 	}

	/*---States stay the same in all phases---*/
	for (i=0; i<=*dimx2; i++) {
		x_2[i] = x_1[i];
	}

	/*---Convert mee_eci states to mee_lci states---*/
	if (*phase == 2) {
		mee_coords_ECI			= dvector(1,6);
		mee_coords_LCI			= dvector(1,6);
		kep_coords_LCI			= dvector(1,6);
		cart_inert_coords_ECI	= dvector(1,6);
		cart_inert_coords_LCI	= dvector(1,6);
		
		for (i=1; i<=6; i++) {
			mee_coords_ECI[i]			= x_1[i-1];
			mee_coords_LCI[i]			= 0.0;
			kep_coords_LCI[i]			= 0.0;
			cart_inert_coords_ECI[i]	= 0.0;
			cart_inert_coords_LCI[i]	= 0.0;
			ephemerides_moon[i-1]		= 0.0;
		}
		mee2cart(mee_coords_ECI, gmue[EARTH], cart_inert_coords_ECI); 
		
		julian_date = rpar_1[0] + EPOCH + x_1[6];  // = start time + current mission time [days] = current time
		target = MOON;
		center = EARTH;
		compute_ephemeris(julian_date, target, center, ephemerides_moon);
		
		// transform from ECI to LCI
		for (i=1; i<=6; i++) {
			cart_inert_coords_LCI[i] = cart_inert_coords_ECI[i] - ephemerides_moon[i-1]*1.0e3;
		}
		cart2mee(cart_inert_coords_LCI, gmue[MOON], mee_coords_LCI);
		
		for (i=1; i<=6; i++) {
			x_2[i-1] = mee_coords_LCI[i];
		}
		
		// clean up
		free_dvector(mee_coords_ECI,1,6);
		free_dvector(mee_coords_LCI,1,6);
		free_dvector(kep_coords_LCI,1,6);
		free_dvector(cart_inert_coords_ECI,1,6);
		free_dvector(cart_inert_coords_LCI,1,6);
	}
}