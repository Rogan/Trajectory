/* matrices.h contains the prototypes

	dmat_prod()
	dmat_times_vec()
	dmat_transp()
	print_matrix()
	rot_x()
	rot_y()
	rot_z()
	
 */


/**************************************************************************************/

#ifndef _MATRICES_H_
#define _MATRICES_H_

void dmat_prod (
	double **A,
	double **B,
	double **C,
	int u,
	int v,
	int w
	);

void dmat_times_vec (
	double **A,
	double *b,
	double *c,
	int u,
	int v
	);

void dmat_transp(
	double **in_mat,
	double **out_transp_mat,
	int u,
	int v
	);

void print_matrix(
	double **mat,
	int u,
	int v
	);

void rot_x(
	double **rotmatrix,
	double angle
	);


void rot_y(
	double **rotmatrix,
	double angle
	);

void rot_z(
	double **rotmatrix,
	double angle
	);


#endif //_MATRICES_H_




