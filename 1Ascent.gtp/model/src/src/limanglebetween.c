/* This file contains the functions 

limanglebetween()

This function is used to reduce an angle by multiples of 2*pi until it is between 
phi_inf and phi_sup. */

#include <stdio.h>
#include <math.h>
#include "limanglebetween.h"
#include "constants.h"

/***********************************************************************************/

void limanglebetween(
	double *angle,
	const double phi_inf,
	const double phi_sup
	)
{

	/* This function changes 'angle' to fit between 'phi_inf' and 'phi_sup' by adding or
		subtracting multiples of 2*pi. If it is not possible to fit 'angle' between these
		borders, an error message will be issued. phi_sup must be greater than phi_inf.

	Input:
		angle	[rad] angle that will be changed to fulfil phi_inf <= angle <= phi_sup
		phi_inf [rad] upper border for angle
		phi_sup [rad] lower border for angle	phi_sup>=phi_inf!

	Output:
		angle	[rad] 

    */

	if (phi_inf > phi_sup) {
		perror("phi_sup must be greater than phi_inf in limanglebetween()\n");
		return;
	}


	//instructions


	
	if (*angle < phi_inf)
		*angle = *angle + 2*PI* ceil((phi_inf-*angle)/(2*PI));
		
	// Now angle is > phi_inf
	
	
	if (*angle > phi_sup) 
		*angle = *angle - 2*PI * ceil((*angle-phi_sup)/(2*PI));
		
	// Now angle is < phi2
	
	
	// if at this point again angle < phi_inf, then
	// the relation phi_inf <= angle <= phi_sup cannot be fulfilled by
	// adding multiples of 2*pi.
	
	if (*angle < phi_inf)
		perror("angle cannot be brought within [phi_inf,phi_sup] in limanglebetween()\n");


}
