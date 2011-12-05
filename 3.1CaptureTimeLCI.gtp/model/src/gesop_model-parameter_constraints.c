/* -------------------------------------------------------------------------------
 * Gesop_Model.Parameter_Constraints
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

void __cdecl Parameter_Constraints(
	const int             *phase,		// Current phase number
	const int             *dimip,		// Dimension of integer parameter vector
	const int             *dimrp,		// Dimension of real parameter vector
	const Phase_Info_Type *fazinf,		// Record with additional phase info
	const int             ipar[],		// Integer parameter vector of phase
	const double          rpar[],		// Real parameter vector of phase
 	const int             *dimprc,		// Dimension of parameter constraint vector
	const Decision_Type   evalc[],		// Which constraints to evaluate
	double          	  *prcon,		// Vector of parameter constraints
	int             	  *error		// Error flag
	)
{
	/* computes the nonlinear parameter constraints of the specified
	 * phase as a function of the phase_info and the real parameter
	 * vector of that phase. */

	int i;

#ifdef DEBUG_MODE
printf("Parameter Constraints\n");_flushall();
#endif

	*error	= 0;

	/* default initialization  */
	for (i=0; i<*dimprc; i++) {
 		prcon[i] = 0.0;
 	}

	/* deltaL should be positive */
	if(evalc[0] == yes)
		prcon[0] = rpar[2];
}