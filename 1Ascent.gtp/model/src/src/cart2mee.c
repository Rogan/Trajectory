/*	This file contains the function 
	
	cart2mee()

*/

#include <math.h>
#include "cart2mee.h"
#include "vectors.h"
#include "nrutil.h"
#include "constants.h"
#include "limanglebetween.h"

/*****************************************************************************/

void cart2mee(
	double *cart,
	const double mue,
	double *out_mee
	)
{	/*	This function converts cartesian coordinates to modified equinoctial elements

		It is assumed that cart and out_kep are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		cart = dvector(1,6), out_kep = dvector(1,6)

	input:
	cart[6]			cartesian elements in m and m/s, being
		cart[1]		x
		cart[2]		y
		cart[3]		z
		cart[4]		x_dot
		cart[5]		y_dot
		cart[6]		z_dot

	mue				gravitational parameter of central body /(m³/s²)

	ouput:
	mee[6]			modified equinoctial elements being
		mee[1]		p
		mee[2]		f
		mee[3]		g
		mee[4]		h
		mee[5]		k
		mee[6]		L */

	
	double	radius,vel,hmag,rdotv,rzerod;
	double	ssqrd, sinl, cosl;
	double	p,f,g,h,k,L;
	double	*hvec, *hvec_norm, *uhat, *vhat, *eccen;
	double	*fhat, *ghat;
	double	*rvec, *vvec;
	int		i;

	hvec		= dvector(1,3);
	hvec_norm	= dvector(1,3);
	uhat		= dvector(1,3);
	vhat		= dvector(1,3);
	eccen		= dvector(1,3);
	fhat		= dvector(1,3);
	ghat		= dvector(1,3);
	rvec		= dvector(1,3);
	vvec		= dvector(1,3);
	
	for (i=1; i<=3; i++){
		hvec[i]			= 0.0;
		hvec_norm[i]	= 0.0;
		uhat[i]			= 0.0;
		vhat[i]			= 0.0;
		eccen[i]		= 0.0;
		fhat[i]			= 0.0;
		ghat[i]			= 0.0;
		rvec[i]			= 0.0;
		vvec[i]			= 0.0;
	}

	// split cartesian coords in radius and velocity-vector	
	for (i=1; i<=3; i++){
		rvec[i] = cart[i];
		vvec[i] = cart[i+3];
	}

	// compute angular momentum and eccentricity vectors
	dcross_prod(rvec, vvec, hvec);
	dcross_prod(vvec, hvec, eccen);

	// compute hmag and pmee
	radius	= dvec_abs(rvec,3);	
	vel		= dvec_abs(vvec,3);	
	hmag	= dvec_abs(hvec,3);	
	p		= hmag*hmag/mue;

	rdotv	= ddot_prod(rvec, vvec);
	rzerod	= rdotv/radius;



	for (i=1;i<=3;i++) {
        uhat[i] = rvec[i]/radius;
		vhat[i] = (radius*vvec[i] - rzerod*rvec[i])/hmag;
		
		eccen[i] = eccen[i]/mue - uhat[i];
	}

	// compute hvec norm vector
	for (i=1;i<=3;i++) {
		hvec_norm[i] = hvec[i]/hmag;
	}	
	
	// compute xkmee and hmee
	k = hvec_norm[1]/(1.0 + hvec_norm[3]);
	h = -hvec_norm[2]/(1.0 + hvec_norm[3]);
	
	// construct unit vectors in the equinoctial frame
	fhat[1] = 1.0 - k*k+ h*h;
    fhat[2] = 2.0*k*h;
	fhat[3] = -2.0*k;
	
	ghat[1] = fhat[2];
	ghat[2] = 1.0 + k*k - h*h;
	ghat[3] = 2.0*h;
	
	ssqrd = 1.0 + k*k + h*h;
	
	// normalize
	
	for (i=1;i<=3;i++) {
		fhat[i] = fhat[i]/ssqrd;
		ghat[i] = ghat[i]/ssqrd;
	}
	
	// compute fmee and gmee
	f = eccen[1]*fhat[1] + eccen[2]*fhat[2] + eccen[3]*fhat[3];
	g = eccen[1]*ghat[1] + eccen[2]*ghat[2] + eccen[3]*ghat[3];
	
	// compute true longitude
	cosl = uhat[1] + vhat[2];
	sinl = uhat[2] - vhat[1];
	L = atan2(sinl,cosl);
	limanglebetween(&L, 0.0, 2.0*PI);


	out_mee[1] = p;
	out_mee[2] = f;
	out_mee[3] = g;
	out_mee[4] = h;
	out_mee[5] = k;
	out_mee[6] = L;

	free_dvector(hvec,1,3);
	free_dvector(hvec_norm,1,3);
	free_dvector(uhat,1,3);
	free_dvector(vhat,1,3);
	free_dvector(eccen,1,3);
	free_dvector(fhat,1,3);
	free_dvector(ghat,1,3);
	free_dvector(rvec,1,3);
	free_dvector(vvec,1,3);

}
