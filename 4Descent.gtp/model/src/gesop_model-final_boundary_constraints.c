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
#include <math.h>
#include "gesop_model.h"
#include "kepler.h"
#include "nrutil.h"
#include "vectors.h"
#include "kepler-mee.h"
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

	
	double	*mee_coords_LCI;			// contains the first six elements of the state in the equ. elements
	double	*kep_coords_LCI;			// contains the first six elements of the state in keplerian elements
	int		i;


#ifdef DEBUG_MODE
printf("Final Boundary Constraints\n");_flushall();
#endif
	
	*error = 0;

	// Initialise storage
	mee_coords_LCI			= dvector(1,6);				
	kep_coords_LCI			= dvector(1,6);				
				
	for (i=1; i<=6; i++) {
		mee_coords_LCI[i]			= x[i-1];
		kep_coords_LCI[i]			= 0.0;
	}

	// Calculate parameters
	mee2kep(mee_coords_LCI, kep_coords_LCI);

	/* Semimajor axis = 1,838km */
	if (evalc[0] == yes) {
		bcon[0] = 1.0 - kep_coords_LCI[1]/1.8371e6;
	}

	/* Eccentricity < 0.1 */
	if (evalc[1] == yes) {
		bcon[1] = 0.1 - kep_coords_LCI[2];
	}

	/* Inclination = 70° */
	if (evalc[2] == yes) {
		bcon[2] = kep_coords_LCI[3]*R2D/70 - 1.0;
	}

	/* Arg periapsis = 42.296° */
	if (evalc[3] == yes) {
		bcon[3] =  kep_coords_LCI[4]*R2D/42.296 - 1.0;
	}

	/* LLAN (RAAN) = 285° */
	if (evalc[4] == yes) {
		bcon[4] =  kep_coords_LCI[5]*R2D/285.0 - 1.0;
	}

	/* Anomaly = 0.0° */
	if (evalc[5] == yes) {
		bcon[5] =  x[7]/PI;
	}


	// Free memory
	free_dvector(mee_coords_LCI,1,6);
	free_dvector(kep_coords_LCI,1,6);

}