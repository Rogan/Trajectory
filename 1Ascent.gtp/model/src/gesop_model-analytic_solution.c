/* -------------------------------------------------------------------------------
 * Gesop_Model.Analytic_Solution
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

void __cdecl Analytic_Solution(
	const int             *phase,		// Current phase number
	const int             *dimx,		// dimension of state vector X
	const int             *dimip,		// dimension of integer parameter vector
	const int             *dimrp,		// dimension of real parameter vector
	const Phase_Info_Type *fazinf,		// Record with additional phase info
	const double          *t_1,			// Initial time t1
	const double          x_1[],		// State vector X at time t1
	const int             ipar[],		// Integer parameter vector of phase
	const double          rpar[],		// Real parameter vector of phase
	const double       	  *t_2,			// Solution time t2
	double          	  *x_2,			// State vector X at time t2
	int             	  *error)		// error flag
{
	int i;

#ifdef DEBUG_MODE
printf("Analytic Solution\n");_flushall();
#endif

	*error = 0;

	for (i=0; i<*dimx; i++) {
 		x_2[i] = 0.0;
 	}
}