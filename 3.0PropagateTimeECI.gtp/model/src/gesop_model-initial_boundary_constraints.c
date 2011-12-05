/* -------------------------------------------------------------------------------
 * Gesop_Model.Initial_Boundary_Constraints
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

void __cdecl Initial_Boundary_Constraints(
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
	const int				*dimbc,		// dimension of boundary constraint vector
	const Decision_Type		evalc[],	// Which constraints to evaluate
	double          		*bcon,		// Vector of boundary constraints
 	int             		*error		// Error flag
	)
{
	/* computes the boundary constraints associated with the initial
	 * time of the optimal control problem (t0 of *phase)
	 * as a function of the initial time, the state- and control
	 * vector at time t0 and the real parameter vector of *phase. */
	
	double	*mee_coords_LCI;			// contains the first six elements of the state, ie the equ. elements
	double	*kep_coords_LCI;			// contains the first six elements of the state in earth centred keplerian elements
	int		i;

#ifdef DEBUG_MODE
printf("Initial Boundary Constraints\n");_flushall();
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


		/* p */
		if (evalc[0] == yes) {
			bcon[0] = mee_coords_LCI[1]/3.32392054e7 - 1.0;
		}

		/* f */
		if (evalc[1] == yes) {
			bcon[1] = mee_coords_LCI[2] + 0.288;
		}

		/* g */
		if (evalc[2] == yes) {
			bcon[2] = mee_coords_LCI[3] - 0.7688;
		}

		/* h */
		if (evalc[3] == yes) {
			bcon[3] = mee_coords_LCI[4] + 0.007;
		}

		/* k */
		if (evalc[4] == yes) {
			bcon[4] = mee_coords_LCI[5] - 0.002;
		}

		/* L */
		if (evalc[5] == yes) {
			bcon[5] = mee_coords_LCI[6]/701.3 - 1.0;
		}

		/* t */
		if (evalc[6] == yes) {
			bcon[6] =  x[6]/274.395 - 1.0;
		}

		/* m */
		if (evalc[7] == yes) {
			bcon[7] =  x[7]/181.98 - 1.0;
		}


		// Free memory
		free_dvector(mee_coords_LCI,1,6);
		free_dvector(kep_coords_LCI,1,6);

}