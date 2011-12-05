/*	mee2cart.c
	This file contains the functions
	
	mee2cart()

*/

#include <stdio.h>
#include <math.h>
#include "mee2cart.h"

/*****************************************************************************/

void mee2cart(
	double *mee,
	const double mue,
	double *out_cart
	)
{

	/*	This function converts modified equinoctial elements to cartesian coordinates

		It is assumed that mee and out_cart are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		mee = dvector(1,6), out_cart = dvector(1,6)

		If (1 + f*cos(L) + g*sin(L)) equals zero, the position is undefined; presumably this is impossible?
		If p <= zero, the velocity is undefined.


	input:
	mee[6]			modified equinoctial elements being
		mee[1]		p
		mee[2]		f
		mee[3]		g
		mee[4]		h
		mee[5]		k
		mee[6]		L 

	mue				gravitational parameter of central body /(m³/s²)

	ouput:
	out_cart[6]		cartesian elements in m and m/s, being
		out_cart[1]	x
		out_cart[2]	y
		out_cart[3]	z
		out_cart[4]	x_dot
		out_cart[5]	y_dot
		out_cart[6]	z_dot */



	//local variables
	double p,f,g,h,k,L;
	double xi,r,alpha_2,s_2;
	

	//renames
	p = mee[1];
	f = mee[2];
	g = mee[3];
	h = mee[4];
	k = mee[5];
	L = mee[6];

	//auxiliary scalars
	xi = 1 + f*cos(L) + g*sin(L);
	if (xi == 0.0) {
		perror("mee2cart: radius undefined\n");_flushall();
	} else {
		r = p/xi;
	}
	alpha_2 = pow(h,2.0) - pow(k,2.0);
	s_2		= 1 + pow(h,2.0) + pow(k,2.0);

	out_cart[1] = r/s_2 * (cos(L) + alpha_2*cos(L) + 2*h*k*sin(L));
	out_cart[2] = r/s_2 * (sin(L) - alpha_2*sin(L) + 2*h*k*cos(L));
	out_cart[3] = 2.0*r/s_2 * (h*sin(L) - k*cos(L));

	if(p <= 0.0) {
		perror("mee2cart: velocity undefined\n");_flushall();
	} else {
		out_cart[4] = -1.0/s_2 * sqrt(mue/p) * ( sin(L) + alpha_2*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + alpha_2*g);
		out_cart[5] = -1.0/s_2 * sqrt(mue/p) * (-cos(L) + alpha_2*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + alpha_2*f);
		out_cart[6] =  2.0/s_2 * sqrt(mue/p) * (h*cos(L) + k*sin(L) + f*h + g*k);
	}
}
