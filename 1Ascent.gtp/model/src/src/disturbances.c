/*	This file contains the functions

		ir_itheta_ih_transp()
		Fqk()
		third_body_pert_cart()
		AccelEarthObl()
		AccelHarmonic()
		solar_rad_pressure_cart()
		shadowfunc()
		sunangle()
		AccelEarthObl2()
*/

#include <stdio.h>
#include <math.h>
#include "disturbances.h"
#include "vectors.h"
#include "matrices.h"
#include "kepler.h"
#include "nrutil.h"
#include "constants.h"
#include "problem.h"

/**************************************************************************************/

void ir_itheta_ih_transp(
	double	*cart_coords,
	double	**out_rotmatrix
	)
{
	/*	This function computes the rotation matrix Q' = [ir itheta ih]' that represents
		the principal axes of the rotating radial frame in which the disturbing accelerations
		of the modified equinoctial elements are expressed.

		It is assumed that cart_coords is a double vector previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of
		cart_coords = dvector(1,6)
		It is assumed that out_rotmatrix is a double matrix previously allocated with dmatrix()
		from nrutil.c (see: Numerical Recipes in C), in form of
		out_rotmatrix = dmatrix(1,3,1,3)

		Input:
			cart_coords[1] = x			in m
			cart_coords[2] = y
			cart_coords[3] = z
			cart_coords[4] = x_dot		in m/s
			cart_coords[5] = y_dot
			cart_coords[6] = z_dot


		Output:
			out_rotmatrix is a 3 x 3 matrix:

			out_rotmatrix	= [ir itheta ih]
							= [ r/abs(r) , (r x v) x r /((abs(r x v) * abs(r)) , r x v / abs(r x v)  ]

			where r = [ x y z ]' and v = [ x_dot y_dot z_dot]'

			ir is pointing in radial direction, itheta is normal to the radial direction in the orbit plane
			in velocity increasing direction, ih is normal to the orbit plane.

	*/

	//local variables
	double	*r;			//contains the radius vector
	double	*v;			//contains the velocity vector
	double	*rxv;		//contains the crossproduct of r and v
	double	norm_r;		//contains the norm of r
	double	norm_rxv;	//contains the norm of r x v


	double	*ir;
	double	*itheta;
	double	*ih;

	int		i; //counter variable



	//initialization
	r = dvector(1,3);
	v = dvector(1,3);
	rxv = dvector(1,3);
	ir = dvector(1,3);
	itheta = dvector(1,3);
	ih = dvector(1,3);

	for (i=1; i<=3; i++){
		r[i] = cart_coords[i];
		v[i] = cart_coords[i+3];
		rxv[i] = 0.0;
		ir[i] = 0.0;
		itheta[i] = 0.0;
		ih[i] = 0.0;
	}



	//instructions

	norm_r = dvec_abs(r,3); //obtain norm of r
	dvec_times_scalar(r, 3, 1.0/norm_r, ir); //compute ir as r/abs(r)

	dcross_prod(r,v,rxv); //compute r x v
	norm_rxv = dvec_abs(rxv,3); //compute norm of r x v

	dvec_times_scalar(rxv,3, 1.0/norm_rxv, ih); //compute ih as r x v /abs(r x v)

	dcross_prod(rxv,r,itheta); //compute (r x v) x r
	dvec_times_scalar(itheta, 3 ,1.0/(norm_r*norm_rxv), itheta); //divide itheta by abs(r x v)*abs(r)



	out_rotmatrix[1][1] = ir[1];
	out_rotmatrix[1][2] = ir[2];
	out_rotmatrix[1][3] = ir[3];

	out_rotmatrix[2][1] = itheta[1];
	out_rotmatrix[2][2] = itheta[2];
	out_rotmatrix[2][3] = itheta[3];

	out_rotmatrix[3][1] = ih[1];
	out_rotmatrix[3][2] = ih[2];
	out_rotmatrix[3][3] = ih[3];




	//clean up
	free_dvector(r,1,3);
	free_dvector(v,1,3);
	free_dvector(rxv,1,3);
	free_dvector(ir,1,3);
	free_dvector(itheta,1,3);
	free_dvector(ih,1,3);

}

/**************************************************************************************/

void Fqk(
	double	*r_vec,
	double	*s_vec,
	double	*out_fqk
	)
{
	/*	This function computes the scalar F(q_k), used for the computation of 3rd body
		perturbations. It is first defined in Battin, R.H., "An Introduction to the Mathematics and Methods of Astrodynamics"
		AIAA Education Series, New York, 1987.
		Here it was taken directly from Betts, J.T. and Erb, S.O., "Optimal Low Thrust
		Trajectories to the Moon", SIAM J. Applied Dynamical Systems, Vol.2, No.2, pp.144-170

		It is assumed that r_vec and s_vec are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of
		r_vec = dvector(1,3), s_vec = dvector(1,3)


		Input:
			r_vec[1] = x_sat		//cartesian coordinates of the satellite with respect to
			r_vec[2] = y_sat		//primary body
			r_vec[3] = z_sat

			s_vec[1] = x_sec		//cartesian coordinates of the secondary body with respect to
			s_vec[2] = y_sec		//primary body
			s_vec[3] = z_sec

		Output:
			out_fqk					scalar quantity, fqk = qk*( (3+3qk+qk²) / (1+(1+qk)^(3/2)) )
									and qk = r_vec'*(r_vec-2*s_vec) / (s_vec'*s_vec)

	*/


	//local variables
	double qk; //see def. of Output
	double num, den; //numerator and denominator of fraction

	int i; //counter variable


	//initializations
	qk = 0.0;
	num = 0.0;
	den = 0.0;



	//instructions

	for (i=1; i<=3; i++){
		num += r_vec[i] * (r_vec[i]-2*s_vec[i]);
		den += pow(s_vec[i],2.0);
	}

	qk = num/den;

	*out_fqk = qk*( (3.0 + 3.0*qk + pow(qk,2.0)) / (1.0 + pow(1.0+qk,1.5)) );

}

/**************************************************************************************/

void third_body_pert_cart(
		double	julian_date,
		int		center,
		int		target,
		double	*r_vec,
		double	*out_pert_vec
		)
{
	/*	This function computes the perturbing acceleration vector due to the influence
		of planets other than the primary (central) body.
		The algorithm is based on Betts, J.T. "Optimal Interplanetary Orbit Transfers by
		Direct Transcription", The Journal of the Astronautical Sciences, Vol. 12, No.3,
		September 1994, pp.247-268 and corresponds to that introduced by Battin, R.H.
		It uses the JPL Ephemerides (see kepler.h) to determine the vector s_vec from the primary
		to the other body, then computes F(qk) by means of the function Fqk() coded above
		and returns -µk/dk³*(r_vec + F(qk)*s_vec), where µk is the gravitational parameter
		of the third body and dk is the scalar distance from the satellite to that body.

		It is assumed that r_vec and out_pert_vec are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of
		r_vec = dvector(1,3), out_pert_vec = dvector(1,3)

		Input:
			julian_date		julian date of the evaluation
			center			center body (primary)
			target			third body that causes the perturbation (for numerical values,
							cf. "kepler.h" or the documentation to the ephemerides)

			r_vec[1] = x_sat	cartesian coordinates of the satellite with respect to
			r_vec[2] = y_sat	primary body (inertial) [m]
			r_vec[3] = z_sat

		Output:
			out_pert_vec[1]		vector of perturbing acceleration in cartesian coordinates
			out_pert_vec[2]		(inertial) [m/s²]
			out_pert_vec[3]
	*/


	//local variables
	double	ephemerides[6];		// contains the output of compute_ephemeris()
	double	*central2third;		// cartesian vector from central to third body
	double	*third2sat;			// cartesian vector from third body to satellite
	double	norm_third2sat;		// scalar distance from third body to satellite
	double	Fq_k;				// scalar value of F(qk)

	int		i;					//counter variable



	//initializations
	central2third = dvector(1,3);
	third2sat = dvector(1,3);

	for (i=0; i<6; i++){
		ephemerides[i] = 0.0;
	}

	for (i=1; i<=3; i++){
		central2third[i] = 0.0;
		third2sat[i] = 0.0;
	}

	norm_third2sat = 0.0;
	Fq_k = 0.0;

	
	//instructions
	compute_ephemeris(julian_date, target, center, ephemerides); //third body with respect to central body (JPL Ephemerides)

	for (i=1; i<=3; i++){
		central2third[i] = ephemerides[i-1]*1.0e3;	// copy the position of third body from primary [m]
		third2sat[i] = r_vec[i] - central2third[i];	// compute the position of satellite from third body [m]
	}

	norm_third2sat = dvec_abs(third2sat,3);			//compute scalar distance from satellite to third body

	Fqk(r_vec, central2third, &Fq_k);				//compute F(qk) with external function

	for (i=1; i<=3; i++){
		out_pert_vec[i] = -gmue[target]/pow(norm_third2sat,3.0)*(r_vec[i] + Fq_k*central2third[i]);
		//gmue[target] is gravitational parameter of target body defined in "kepler.h"
	}

	//clean up
	free_dvector(central2third,1,3);
	free_dvector(third2sat,1,3);

}

/**************************************************************************************/

void AccelEarthObl(		
		double	*r_vec,
		double	*out_pert_vec
		)
{
	/*	calculate acceleration in ECI due to earth zonal gravity effects of J2 and J3. Longitude is 
		irrelevant since only zonal harmonics are used, hence no coordinate transformation is required.

		It is assumed that r_vec and out_pert_vec are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of
		r_vec = dvector(1,3), out_pert_vec = dvector(1,3)

		Input:
			r_vec[1] = x_sat	cartesian coordinates of the satellite with respect to
			r_vec[2] = y_sat	primary body (inertial) [m]
			r_vec[3] = z_sat

		Output:
			out_pert_vec[1]		vector of perturbing acceleration in cartesian coordinates
			out_pert_vec[2]		(inertial) [m/s²]
			out_pert_vec[3]
	*/

	//local variables
	double	r;			// distance from earth to satellite
	double	*ez;		// unity vector in z-direction in ECI (which is north direction)
	double	*iN;		// unity vector in north direction in geodetic frame expressed in ECI
	double	*ir;		// unity vector in radius direction in ECI
	double	*gN;		// acceleration vector in north direction in geodetic frame expressed in ECI
	double	gN_mag;		// magnitude of gN vector
	double	*gr;		// acceleration vector in radius direction in ECI
	double	gr_mag;		// magnitude of gr vector
	double	*temp; 
	double	temp_mag;
	double	proj;		// vector ez projected onto ir;
	double	sin_gl;		// sine of geocentric latitude
	double	cos_gl;		// cosine of geocentric latitude
	double	RHO;		// ratio of satellite position to Earth radius
	int		i;			// counter variable

	ez = dvector(1,3);
	iN = dvector(1,3);
	ir = dvector(1,3);
	temp = dvector(1,3);
	gN = dvector(1,3);
	gr = dvector(1,3);

	for (i=1; i<=3; i++){
		ez[i] = 0.0;
		iN[i] = 0.0;
		ir[i] = 0.0;
		temp[i] = 0.0;
		gN[i] = 0.0;
		gr[i] = 0.0;
	}

	r = dvec_abs(r_vec,3);

	for (i=1; i<=3; i++) {
		ir[i] = r_vec[i]/r;
	}

	ez[1] = 0.0;
	ez[2] = 0.0;
	ez[3] = 1.0;

	proj = ddot_prod(ez, ir);
	for (i=1; i<=3; i++) {
		temp[i] = ez[i] - proj*ir[i];
	}

	temp_mag = dvec_abs(temp,3);

	for (i=1; i<=3; i++) {
		iN[i] = temp[i]/temp_mag;
	}

	sin_gl = r_vec[3]/r;
	cos_gl = sqrt(r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2])/r;
	RHO = R_JGM3*1000.0/r;

	gN_mag = -gmue[EARTH]/(r*r)*(3.0*CS_JGM3[2][0]*RHO*RHO*sin_gl*cos_gl-1.5*CS_JGM3[3][0]*pow(RHO,3.0)*(5.0*sin_gl*sin_gl-1.0)*cos_gl);
	gr_mag = gmue[EARTH]/(r*r)*((1.5*CS_JGM3[2][0]*RHO*RHO*(1.0-3.0*sin_gl*sin_gl)+2.0*CS_JGM3[3][0]*pow(RHO,3.0)*(3.0-5.0*sin_gl*sin_gl))*sin_gl);

	for (i=1; i<=3; i++) {
		out_pert_vec[i] = gN_mag*iN[i]-gr_mag*ir[i];
	}


	free_dvector(ez,1,3);
	free_dvector(iN,1,3);
	free_dvector(ir,1,3);
	free_dvector(temp,1,3);
	free_dvector(gN,1,3);
	free_dvector(gr,1,3);
}
/**************************************************************************************/


//------------------------------------------------------------------------------
//
// AccelHarmonic
//
// Purpose:
//
//   Computes the acceleration due to the harmonic gravity field of the
//   central body
//
// Input:
//   r           Satellite position vector in the inertial system [m]
//   E           Transformation matrix to body-fixed system
//   GM          Gravitational coefficient [km³/s²]
//   R_ref       Reference radius [km]
//   CS          Spherical harmonics coefficients (un-normalized)
//   n_max       Maximum degree
//   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
// Output:
//   a_inert     3x1-Vector containing acceleration (a=d²r/dt²) in inertial coordinates [m/s²]
//
//	adapted from O. Montenbruck, E. Gill, "Satellite Orbits", Springer, 2000, enclosed CD-ROM
//------------------------------------------------------------------------------

void AccelHarmonic (
	double			*r,				// satellite position (inertial cartesian)
	double			**E,			// conversion to body fixed coords
	double			GM,				// primary body gravitational constant
	double			R_ref,			// primary body reference radius
	const double	CS[21][21],		// gravitational coefficients
	int				n_max,			// degree
	int				m_max,			// order
	double			*a_inert		// output
	)
{

	// Local variables

	int		n,m;					// Loop counters
	double	r_sqr, rho, Fac;		// Auxiliary quantities
	double	x0,y0,z0;				// Normalized coordinates
	double	ax,ay,az;				// Acceleration vector
	double	C,S;					// Gravitational coefficients
	double	*r_bf;					// Body-fixed position
	double	*a_bf;					// Body-fixed acceleration
	double	**E_transp;				// E transposed (body-fixed to inertial)

	double	**V;					// Harmonic functions
	double	**W;					// work array (0..n_max+1,0..n_max+1)



	//Initializations
	r_bf = dvector(1,3);
	a_bf = dvector(1,3);
	E_transp = dmatrix(1,3,1,3);

	// V and W intermediate coefficients required to one index higher than accelerations (p.68)
	V = dmatrix(0,n_max+1,0,n_max+1);
	W = dmatrix(0,n_max+1,0,n_max+1);


	for (n=1; n<=3; n++){
		for (m=0; m<= 3; m++){
			E_transp[n][m] = 0.0;
		}
		r_bf[n] = 0.0;
		a_bf[n] = 0.0;
	}

	for (n=0; n<= n_max+1; n++){
		for (m=0; m<= n_max+1; m++){
			V[n][m] = 0.0;
			W[n][m] = 0.0;
		}
	}


	// convert position to body-fixed coords
	dmat_times_vec(E,r,a_bf,3,3);			// r_bf = E * r;
	dvec_times_scalar(a_bf,3,1.0e-3,r_bf);	// convert to [km]


	// Auxiliary quantities
	r_sqr =  pow(r_bf[1],2.0) + pow(r_bf[2],2.0) + pow(r_bf[3],2.0);    // Square of distance [km²]
	rho   =  R_ref*R_ref / r_sqr;

	// Normalized, body-fixed coords of satellite
	x0 = R_ref * r_bf[1] / r_sqr;
	y0 = R_ref * r_bf[2] / r_sqr;
	z0 = R_ref * r_bf[3] / r_sqr;

	// Evaluate harmonic functions p.66 eq (3.27)
	//   V_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * cos(m*lambda)
	// and
	//   W_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * sin(m*lambda)
	// up to degree and order n_max+1
	//

	// Calculate zonal terms V[n,0]; set W[n,0]=0.0
	// p.67 (3.31)
	V[0][0] = R_ref / sqrt(r_sqr);
	W[0][0] = 0.0;

	// p.67 (3.30) Need this here because of [n-2] index in equation
	V[1][0] = z0 * V[0][0];
	W[1][0] = 0.0;

	for (n=2; n<=n_max+1; n++) {
		V[n][0] = ( (2*n-1) * z0 * V[n-1][0] - (n-1) * rho * V[n-2][0] ) / n;
		W[n][0] = 0.0;
	};

	// Calculate tesseral and sectorial terms V[m,m] .. V[n_max+1,m]
	for (m=1; m<=m_max+1; m++) {

		// Calculate diagonal terms p.66 (3.29)
		V[m][m] = (2*m-1) * ( x0*V[m-1][m-1] - y0*W[m-1][m-1] );
		W[m][m] = (2*m-1) * ( x0*W[m-1][m-1] + y0*V[m-1][m-1] );

		// Calculate first tesseral (3.30)
		if (m<=n_max) {
			V[m+1][m] = (2*m+1)*z0*V[m][m] - 2*m*rho*V[m-1][m];
			W[m+1][m] = (2*m+1)*z0*W[m][m] - 2*m*rho*W[m-1][m];
		};

		// Fill remainder of the tesserals for given m (3.30)
		for (n=m+2; n<=n_max+1; n++) {
			V[n][m] = ( (2*n-1)*z0*V[n-1][m] - (n+m-1)*rho*V[n-2][m] ) / (n-m);
			W[n][m] = ( (2*n-1)*z0*W[n-1][m] - (n+m-1)*rho*W[n-2][m] ) / (n-m);
		};
	};

	// Calculate accelerations ax,ay,az (p.68)
	ax = ay = az = 0.0;

	for (m=0; m<=m_max; m++)
		for (n=m; n<=n_max ; n++)
			if (m==0) {
				C = CS[n][0];	// = C_n,0 ... see definition of CS in const.c
				ax -=       C * V[n+1][1];
				ay -=       C * W[n+1][1];
				az -= (n+1)*C * V[n+1][0];
			}
		else {
			C = CS[n][m];		// = C_n,m ... see definition of CS in const.c
			S = CS[m-1][n];		// = S_n,m ... see definition of CS in const.c
			Fac = 0.5 * (n-m+1) * (n-m+2);
			ax +=   + 0.5 * ( - C * V[n+1][m+1] - S * W[n+1][m+1] )
				    + Fac * ( + C * V[n+1][m-1] + S * W[n+1][m-1] );
			ay +=   + 0.5 * ( - C * W[n+1][m+1] + S * V[n+1][m+1] )
				    + Fac * ( - C * W[n+1][m-1] + S * V[n+1][m-1] );
			az += (n-m+1) * ( - C * V[n+1][m]   - S * W[n+1][m]   );
		};

	// Body-fixed acceleration
	a_bf[1] = ax;
	a_bf[2] = ay;
	a_bf[3] = az;

	dvec_times_scalar(a_bf,3,GM/(R_ref*R_ref), a_bf);	// a_bf = (GM/(R_ref*R_ref)) * Vector(ax,ay,az);

	// Convert acceleration back to space-fixed coords
	dmat_transp(E,E_transp,3,3);
	dmat_times_vec(E_transp, a_bf, a_inert, 3,3);  // return  Transp(E)*a_bf;


	// Clean-up
	free_dvector(r_bf,1,3);
	free_dvector(a_bf,1,3);

	free_dmatrix(E_transp,1,3,1,3);
	free_dmatrix(V,0,n_max+1,0,n_max+1);
	free_dmatrix(W,0,n_max+1,0,n_max+1);


}

/**************************************************************************************/

void solar_rad_pressure_cart(
	double  julian_date,
	double  mass,
	double  *r_vec,
	int		center,
	double	*out_pert_vec
	)
{
	/*	This function computes the perturbing acceleration due to solar radiation
		pressure in inertial cartesian coordinates.
		out_pert_vec = - shadowfunc * P_sun * (1+eps) * A_sat/m_sat * r_sat2sun/abs(r_sat2sun)

		It is assumed that r_vec and out_pert_vec are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of
		r_vec = dvector(1,3), out_pert_vec = dvector(1,3)

		Input:
			julian_date			julian date of the evaluation

			mass				mass of satellite at julian_date [kg]

			r_vec[1] = x_sat	cartesian coordinates of the satellite with respect to
			r_vec[2] = y_sat	primary body (inertial) [m]
			r_vec[3] = z_sat

			center				primary body (as per kepler.h, myTransfer.h)

		Output:
			out_pert_vec[1]		vector of perturbing acceleration in cartesian coordinates
			out_pert_vec[2]		(inertial) [m/s²]
			out_pert_vec[3]

		Rogan edit: made generic by adding centre input parameter. Formerly assumed moon-centred.
	*/


	//local variables
	double	ephemerides[6];		//contains the state of the sun rel. to the "center" body
	int		target;				//target for compute_ephemeris
	double	*sat2sun;			//vector from satellite to the sun
	double	*body2sun;			//vector from "center" body to the sun
	double	norm_sat2sun;		//scalar distance from satellite to the sun
	double	sol_press;			//solar pressure at position of satellite
	double	norm_out_pertvec;	//scalar value of perturbing acceleration divided by distance from satellite to sun

	int		i;					//counter variable


	//initialization
	sat2sun = dvector(1,3);
	body2sun = dvector(1,3);


	for (i=0; i<6; i++){
		ephemerides[i] = 0.0;
	}

	for (i=1; i<=3; i++){
		sat2sun[i] = 0.0;
		body2sun[i] = 0.0;
	}

	norm_sat2sun = 0.0;
	norm_out_pertvec = 0.0;
	sol_press = 0.0;


	//compute vector from "center" body to sun and from satellite to sun
	target = SUN;
	compute_ephemeris(julian_date, target, center, ephemerides); //sun with respect to primary body (JPL Ephemerides)

	for (i=1; i<=3; i++){
		body2sun[i] = ephemerides[i-1]*1.0e3;	// copy the position of sun from "center" body [m]
		sat2sun[i] = body2sun[i] - r_vec[i];	// calculate the position of sun from satellite [m]
	}
	norm_sat2sun = dvec_abs(sat2sun,3);  

	//compute perturbing acceleration
	sol_press = SOLARLUMINOSITY/(4 * PI * pow(norm_sat2sun,2.0) * SPEEDOFLIGHT);
	norm_out_pertvec = - sol_press * (1+reflectivity) * aref/mass /norm_sat2sun;// * shadowfunc(body2sun, r_vec);
	dvec_times_scalar(sat2sun, 3, norm_out_pertvec, out_pert_vec);

	//clean-up
	free_dvector(sat2sun,1,3);
	free_dvector(body2sun,1,3);
}



/**************************************************************************************/

double shadowfunc(
	double	*body2sun,
	double	*r_vec,
	double	R_body
	)
{
	/*	This function returns the shadow function for a satellite whose view of the sun
		is occulted by the central body.

		It is assumed that sat2sun, body2sun and r_vec are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of
		sat2sun = dvector(1,3), body2sun = dvector(1,3), r_vec = dvector(1,3)

		The used algorithms are taken from O. Montenbruck and E. Gill, "Satellite Orbits"


		Input:
			
			body2sun[1] = x_sun-x_body	cartesian coordinates of the sun with respect to
			body2sun[2] = y_sun-y_body	the central body (inertial) (units unimportant, just direction)
			body2sun[3] = z_sun-z_body

			r_vec[1] = x_sat			cartesian coordinates of the satellite with respect to
			r_vec[2] = y_sat			the central body (inertial)
			r_vec[3] = z_sat

			R_body						average radius of the central body [m]

		Output:
			shadowfunc = 0				if satellite is in umbra
			shadowfunc = 1				if the satellite is in sunlight
			0 < shadowfunc < 1			if the satellite is in penumbra.
	*/


	//local variables
	

	double	norm_body2sun;	// scalar distance from body to sun
	double	norm_r_vec;		// scalar distance from body to satellite 
	double	s0;				//scalar distance between body and fundamental plane
	double	l,l1,l2;		//radii of cones   to sat / of penumbra / of umbra
	double	f1,f2;			//halfangle of shadowcones penumbra / umbra
	double	c1,c2;			//distances from fundamental plane to vertices of shadowcones penumbra / umbra
	double	temp;

	//instructions
	

	// (conical shadow model)
	norm_body2sun = dvec_abs(body2sun,3);
	norm_r_vec = dvec_abs(r_vec,3);

	s0 = -ddot_prod(r_vec,body2sun)/norm_body2sun; // projection of r_vec onto body2sun

	l = sqrt(pow(norm_r_vec,2.0) - pow(s0,2.0) );

	temp = (R_SUN+R_body)/norm_body2sun;
	if(temp > 1 || temp < -1)
		printf("shadowfunc: mathematical violation.\n");
	f1 = asin(temp);

	temp = (R_SUN-R_body)/norm_body2sun;
	if(temp > 1 || temp < -1)
		printf("shadowfunc: mathematical violation.\n");
	f2 = asin(temp);

	c1 = s0 + R_body/sin(f1);
	c2 = s0 - R_body/sin(f2);

	l1 = c1*tan(f1);
	l2 = fabs(c2*tan(f2));	// in our case, the fundamental plane will always between the two
							// vertices of the umbral cones. Note that this implies c2, l2 <0
							// revise this section, if considerably larger orbits are to be
							// examined that involve annular eclipses

	if ( (l>l1) || (l2>l1) )
		return(1.0); //no eclipse
	else {
		if (l<l2)
			return(0.0); //full eclipse (umbra)
		else
			return((l-l2)/(l1-l2)); //penumbra
	}
}


/**************************************************************************************/

double sunangle(
		double *craft2sun,
		const double *u
		)
{
	/*	This function returns the angle between the satellite thrust vector and the sun.

		It is assumed that craft2sun is a double vector previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of
		craft2sun = dvector(1,3)

		Input:
				craft2sun[1] = x_sun-x_craft	cartesian coordinates of the sun with respect to
				craft2sun[2] = y_sun-y_craft	the spacecraft (inertial)
				craft2sun[3] = z_sun-z_craft

				u[0] = radial thrust component
				u[1] = tangential thrust component
				u[2] = out of plane thrust component

		Output:
				angle between thrust vector and sun [rad]

	*/

	double	dot_product;
	double	norm_craft2sun;

	norm_craft2sun = dvec_abs(craft2sun,3);
	
	dot_product = u[0]*craft2sun[1] + u[1]*craft2sun[2] + u[2]*craft2sun[3];

	return acos(dot_product/norm_craft2sun);
}


/**************************************************************************************/


void AccelEarthObl2(		
		double	*r_vec,
		double	*out_pert_vec
		)
{
	/*	calculate acceleration in ECI due to earth zonal gravity effect of J2. Longitude is 
		irrelevant since only zonal harmonics are used, hence no coordinate transformation is required.

		It is assumed that r_vec and out_pert_vec are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of
		r_vec = dvector(1,3), out_pert_vec = dvector(1,3)

		Input:
			r_vec[1] = x_sat	cartesian coordinates of the satellite with respect to
			r_vec[2] = y_sat	primary body (inertial) [m]
			r_vec[3] = z_sat

		Output:
			out_pert_vec[1]		vector of perturbing acceleration in cartesian coordinates
			out_pert_vec[2]		(inertial) [m/s²]
			out_pert_vec[3]
	*/

	//local variables
	double	temp;
	double	norm_r_vec;
	int		i;

	norm_r_vec = dvec_abs(r_vec,3);

	temp = -1.5*gmue[EARTH]*R_JGM3*R_JGM3*CS_JGM3[2][0]*1e6/pow(norm_r_vec,7);

	// Formula from Montenbruck and Gill, "Satellite Orbits" p.173
	for (i=0; i<=3; i++) {
//		out_pert_vec[i] = -1.5*gmue[EARTH]*R_JGM3*R_JGM3*J2/pow(r,7)*((5.0*r_vec[3]*r_vec[3]-r*r)*r_vec[i]-2.0*r_vec[3]*r*r*ez[i-1]); // where ez = [0 0 1];
		out_pert_vec[i] = temp*(5.0*r_vec[3]*r_vec[3]-norm_r_vec*norm_r_vec)*r_vec[i];
	}
	out_pert_vec[3] -= temp*2.0*r_vec[3]*norm_r_vec*norm_r_vec;

}


/**************************************************************************************/