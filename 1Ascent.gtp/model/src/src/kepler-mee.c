/* This file contains the functions 

kep2mee() 
	and 
mee2kep()

These functions are used to convert modified equinoctial elements (mee) to keplerian
elements and vice-versa. */

#include <stdio.h>
#include <math.h>
#include "kepler-mee.h"
#include "constants.h"

/***********************************************************************************/

void kep2mee(
	double *kep,
	double *out_mee
	)
{
	/*	converts keplerian elements to modified equinoctial elements

		It is assumed that kep and out_mee are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		kep = dvector(1,6), out_mee = dvector(1,6)

		No errors can be generated

	input:
	kep[6]			keplerian elements being
		kep[1]		semimajor axis
		kep[2]		eccentricity
		kep[3]		inclination
		kep[4]		argument of periapsis
		kep[5]		longitude of ascending node
		kep[6]		true anomaly

	output:
	out_mee[6]		modified equinoctial elements being
		out_mee[1]	p
		out_mee[2]	f
		out_mee[3]	g
		out_mee[4]	h
		out_mee[5]	k
		out_mee[6]	L */


	//local variables
	double semi, ecc, inc, arg, raan, truean;
	

	//renames
	semi	= kep[1];
	ecc		= kep[2];
	inc		= kep[3];
	arg		= kep[4];
	raan	= kep[5];
	truean	= kep[6];

	out_mee[1] = semi * ( 1.0 - pow(ecc,2.0) );
	out_mee[2] = ecc * cos(arg + raan);
	out_mee[3] = ecc * sin(arg + raan);
	out_mee[4] = tan(inc/2.0) * cos(raan);
	out_mee[5] = tan(inc/2.0) * sin(raan);
	out_mee[6] = arg + raan + truean;

		
}




/***********************************************************************************/

void mee2kep(
	double *mee,		// 6 modified equinoctial elements as input
	double *out_kep			// 6 Kepler elements as output
	)
{
	/*	converts modified equinoctial elements to keplerian elements
		It is assumed that mee and out_kep are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		mee = dvector(1,6), out_kep = dvector(1,6)

		If the inclination is zero, the raan is set to 0.0
		If the eccentricity is zero, the arg of periapsis is set to 0.0
		If the eccentricity is one, the semimajor axis is calculated as undefined
	
	  input:
	  mee[6]		modified equinoctial elements being:
		p
		f
		g
		h
		k
		L

	   
	  output:
      out_kep[6]    keplerian elements being:
	      semi      semimajor axis
		  ecc	    eccentricity
		  inc	    inclination in rad
		  arg	    argument of periapsis in rad
		  omega	    longitude of ascending node in rad
		  truean	true anomaly in rad*/

	   
	
	// local variables
	double p,f,g,h,k,L;
	double semi, ecc, inc, arg, raan, truean;
	extern void limanglebetween(double *, const double, const double);
	double temp1;
	
	
	//instructions


	// renames
	p = mee[1];
	f = mee[2];
	g = mee[3];
	h = mee[4];
	k = mee[5];
	L = mee[6];
	
	temp1 = pow(f,2) + pow(g,2);

	if(temp1 == 1) {
		perror("mee2kep: semimajor axis undefined\n");
	}
	
	//conversion	
	semi = p/(1-temp1);
	
	ecc = sqrt(temp1); 

	inc = 2*atan( sqrt( pow(h,2) + pow(k,2) ) );
	limanglebetween(&inc, 0.0, 2*PI);

	if((f == 0) && (g == 0)) // zero eccentricity; no periapsis
		arg = 0.0;
	else	
		arg = atan2(g,f) - atan2(k,h);
	limanglebetween(&arg, 0.0, 2*PI);

	if((k == 0) && (h == 0)) // zero inclination; no ascending node
		raan = 0.0;
	else 
		raan = atan2(k,h);
	limanglebetween(&raan, 0.0, 2*PI);

	truean = L - atan2(g,f);
	limanglebetween(&truean, 0.0, 2*PI);


	//renames
	out_kep[1] = semi;
	out_kep[2] = ecc;
	out_kep[3] = inc;
	out_kep[4] = arg;
	out_kep[5] = raan;
	out_kep[6] = truean;

}

/***********************************************************************************/