/*	cart2kep.c
	This file contains the function 
	
		cart2kep()

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cart2kep.h"
#include "vectors.h"
#include "nrutil.h"
#include "constants.h"
#include "limanglebetween.h"

/*****************************************************************************/

void cart2kep(
	double *cart,
	const double mue,
	double *out_kep
	)
{

	/*	This function converts cartesian coordinates to keplerian elements

		It is assumed that cart and out_kep are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		cart = dvector(1,6), out_kep = dvector(1,6)

		Error: 
			if angular momentum is zero. 
			if radius is zero.
			if mue is zero.
		Warning:
			if escape velocity is reached (true anomaly calculation breaks).
			if angular momentum is exclusively polar (equatorial orbit, no RAAN; set to zero).
			if velocity in radial direction AND radius equal to semimajor axis (linear orbit).
			if "p" is greater than semimajor axis, eccentricity is undefined. Assuming "p" is semilatus 
				rectum, this means a hyperbolic orbit which should be caught by escape velocity.
			if 1-ecc*ecc equals zero (escape velocity).
			if sin(E) = 0 AND cos(E) = ecc (escape velocity).
			if cart[3] = 0 AND -cart[1]*h_vec_norm[2]+cart[2]*h_vec_norm[1] = 0 ???


	input:
	cart			cartesian elements in m and m/s, being
		cart[1]		x
		cart[2]		y
		cart[3]		z
		cart[4]		x_dot
		cart[5]		y_dot
		cart[6]		z_dot

	mue				gravitational parameter of central body (m³/s²)

	ouput:
	out_kep				keplerian elements in m and rad
		out_kep[1]		semimajor axis
		out_kep[2]		eccentricity
		out_kep[3]		inclination
		out_kep[4]		argument of periapsis
		out_kep[5]		longitude of ascending node
		out_kep[6]		true anomaly 
	*/


	double	*r_vec, *v_vec, *h_vec, *h_vec_norm;
	double	r,v,hmag,p,n,E,u,numerator,denominator;
	double	semi,ecc,inc,omega,raan,nue;
	int		i;

	r_vec		= dvector(1,3);
	v_vec		= dvector(1,3);
	h_vec		= dvector(1,3);
	h_vec_norm	= dvector(1,3);

	for(i=1;i<=3;i++) {
		r_vec[i]		= cart[i];
		v_vec[i]		= cart[i+3];
		h_vec[i]		= 0.0;
		h_vec_norm[i]	= 0.0;
	}

	r = dvec_abs(r_vec,3);
	v = dvec_abs(v_vec,3);
	
	dcross_prod(r_vec, v_vec, h_vec);
	hmag = dvec_abs(h_vec,3);

	/* --- ERROR HANDLING --- */

	if(r == 0) {
		perror("cart2kep: Radius zero.\n");
		exit(1);
	}
	if(mue == 0) {
		perror("cart2kep: Gravitational constant mu is zero.\n");
		exit(1);
	}
	if(hmag == 0) {
		perror("cart2kep: Angular momentum zero.\n");
		exit(1);
	}
	if(v >= sqrt(2*mue/r)) {
		perror("cart2kep: Escape velocity.\n");
	}

	/* --- Okay, back to the good stuff --- */

	for (i = 1; i <=3; i++) {
		h_vec_norm[i] = h_vec[i]/hmag;
	}

	semi = 1.0/(2.0/r - v*v/mue);
	p = hmag*hmag/mue;

	if(p > semi) {
		perror("cart2kep: Eccentricity undefined.\n");
	}

	ecc = sqrt(1.0-p/semi);
	inc = atan2(dvec_abs(h_vec_norm,2),h_vec_norm[3]);
	
	if(h_vec_norm[3] != hmag)
		raan = atan2(h_vec_norm[1],(-h_vec_norm[2]));
	else
		raan = 0.0;

	n = sqrt(mue/pow(semi,3));
	E = atan2(ddot_prod(r_vec,v_vec)/(semi*semi*n),1.0-r/semi);	// numerator can't be zero as long as orbit isn't linear
	nue = atan2(sqrt(1.0-ecc*ecc)*sin(E),cos(E)-ecc);			// can't both be zero as long as eccentricity < 1
	
	numerator = cart[3];
	denominator = -cart[1]*h_vec_norm[2]+cart[2]*h_vec_norm[1];
	
	if((numerator == 0.0) && (denominator == 0.0))	{
		perror("cart2kep: I have no idea when this error may occur.\n");
	}
	u = atan2(numerator, denominator);

	omega = u - nue;	

	limanglebetween(&omega, 0.0, 2.0*PI);
	limanglebetween(&raan, 0.0, 2.0*PI);
	limanglebetween(&nue, 0.0, 2.0*PI);

	out_kep[1] = semi;
	out_kep[2] = ecc;
	out_kep[3] = inc;
	out_kep[4] = omega;
	out_kep[5] = raan;
	out_kep[6] = nue;

	free_dvector(r_vec,1,3);
	free_dvector(v_vec,1,3);
	free_dvector(h_vec,1,3);
	free_dvector(h_vec_norm,1,3);

}