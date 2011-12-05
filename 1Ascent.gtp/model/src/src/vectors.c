/* This file contains the functions

	dcross_prod()
	ddot_prod()
	dvec_add()
	dvec_sub()
	dvec_abs()
	dvec_mult_w_scalar()
	dvec_angle()
	print_vector()

*/

#include <stdio.h>
#include <math.h>
#include "vectors.h"

/**************************************************************************************/

void dcross_prod(
	double *a, 
	double *b, 
	double *result
	)
{
	/*	This function performs the cross product result = a x b.
		It is assumed that a, b and result are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		a = dvector(1,3), b = dvector(1,3), result = dvector(1,3).
		
		NB. ASSUMES 1x3 VECTORS! HOW STUPID IS THAT?
	*/


	result[1] = a[2]*b[3] - b[2]*a[3];
	result[2] = a[3]*b[1] - b[3]*a[1];
	result[3] = a[1]*b[2] - b[1]*a[2];

}

/**************************************************************************************/

double ddot_prod(
	double *a, 
	double *b
	)
{
	/*	This function performs the dot product result = a * b.
		It is assumed that a and b are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		a = dvector(1,3), b = dvector(1,3), result = dvector(1,3).

		NB. ASSUMES 1x3 VECTORS! HOW STUPID IS THAT?

	*/

	return a[1]*b[1] + a[2]*b[2] + a[3]*b[3];

}

/**************************************************************************************/

void dvec_add(
	double *a,
	double *b,
	int dim,
	double *result
	)
{
	/*	This function adds the vectors a and b in result: result = a + b
		It is assumed that a, b and result are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		a = dvector(1,dim), b = dvector(1,dim), result = dvector(1,dim).
	*/


	//local variables
	int i; //counter variable


	//instructions

	for (i=1; i<=dim; i++) {
		result[i] = a[i] + b[i];
	}

}

/**************************************************************************************/

void dvec_sub(
	double *a,
	double *b,
	int dim,
	double *result
	)
{
	/*	This function subtracts the vectors a and b in result: result = a - b
		It is assumed that a, b and result are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		a = dvector(1,dim), b = dvector(1,dim), result = dvector(1,dim).
	*/


	//local variables
	int i; //counter variable


	//instructions

	for (i=1; i<=dim; i++) {
		result[i] = a[i] - b[i];
	}

}

/**************************************************************************************/

double dvec_abs(
	double *vec,
	int dim
	)
{
	/*	This function returns the euclidian 2-norm of the vector 'vec' in out_abs_val
		It is assumed that vec is a double vector previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		vec = dvector(1,dim).
		out_abs_val = sqrt( vec[1]² + vec[2]² +...+ vec[dim]² )
		
	*/


	//local variables
	double out_abs_val; //output variable
	int i; //counter variable

	//instructions
	out_abs_val = 0.0;
	for(i=1; i<=dim; i++){
		out_abs_val += pow(vec[i],2.0);
	}
	out_abs_val = sqrt(out_abs_val);

	return out_abs_val;

}

/**************************************************************************************/

void dvec_times_scalar(
	double *in_vec,
	int dim,
	double scalar,
	double *result
	)
{
	/*	This function multiplies the double vector 'in_vec' with the scalar 'scalar' and
		returns the result in 'result'.

		It is assumed that in_vec and result are double vectors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		in_vec = dvector(1,dim), result = dvector(1,dim).

	*/

	//local variables 
	int i; //counter variable


	for (i=1; i<= dim; i++){
		result[i] = in_vec[i]*scalar;
	}

}

/**************************************************************************************/

double dvec_angle(
	double *a,
	double *b
	)
{
	/* This function returns the angle between two cartesian vectors.

		It is assumed that a and b are double cevtors previously allocated with dvector()
		from nrutil.c (see: Numerical Recipes in C), in form of
		a = dvector(1,3), b = dvector(1,3);
	*/

	double	norm_a;
	double	norm_b;

	norm_a = dvec_abs(a,3);
	norm_b = dvec_abs(b,3);

	return acos(ddot_prod(a,b)/norm_a/norm_b);
}


/**************************************************************************************/


void print_vector(
	double *vec,
	int u
	)
{
	/* This function is for debugging purposes. It prints out a vector with dimension
	   u.
	*/

	int i; //counter variable

	printf("[%f", vec[1]);
	for(i=2; i<=u; i++){
		printf(" %f", vec[i]);
	}
	printf("]\n");_flushall();
}


/**************************************************************************************/