/* -------------------------------------------------------------------------------
 * Gesop_Model.Terminal_Cost
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

void __cdecl Terminal_Cost(
	const int				*phase,		// Current phase number
	const int				*dimx,		// Dimension of state vector X
	const int				*dimip,		// Dimension of integer parameter vector
	const int				*dimrp,		// Dimension of real parameter vector
	const Phase_Info_Type	*fazinf,	// Record with additional phase info
	const double			*t,			// Time t of evaluation
	const double			x[],		// State vector X at time t
	const int				ipar[],		// Integer parameter vector of phase
	const double			rpar[],		// Real parameter vector of phase
	double          		*cost,		// Value of the cost functional
	int             		*error)		// Error flag
{
/* computes the cost associated with the final time of the
* specified phase (tf of the given phase) as a function
* of the time tf, the state vector at time tf and the real
* parameter vector of that phase.  */

#ifdef DEBUG_MODE
printf("Terminal Cost\n");_flushall();
#endif

	*error = 0;

	/* Maximise final mass (minimise fuel used) */
	*cost = -x[7]/200.0;

}