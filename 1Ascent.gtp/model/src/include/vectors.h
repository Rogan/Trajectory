#ifndef _VECTORS_H_
#define _VECTORS_H_

void dcross_prod(
	double *a, 
	double *b, 
	double *result
	);

double ddot_prod(
	double *a, 
	double *b
	);

void dvec_add(
	double *a,
	double *b,
	int dim,
	double *result
	);

void dvec_sub(
	double *a,
	double *b,
	int dim,
	double *result
	);

double dvec_abs(
	double *vec,
	int dim
	);

void dvec_times_scalar(
	double *in_vec,
	int dim,
	double scalar,
	double *result
	);

double dvec_angle(
	double *a,
	double *b
	);

void print_vector(
	double *vec,
	int u
	);

#endif //_VECTORS_H