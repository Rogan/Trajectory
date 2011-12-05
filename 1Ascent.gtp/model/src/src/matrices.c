/* matrices.c contains the functions

	dmat_prod()
	dmat_times_vec()
	dmat_transp()
	print_matrix()
	rot_x()
	rot_y()
	rot_z()

 */

#include <stdio.h>
#include <math.h>
#include "vectors.h"

/**************************************************************************************/

void dmat_prod (
	double **A,
	double **B,
	double **C,
	int u,
	int v,
	int w
	)
{
	/*	This routine computes C = A*B, with dim(A) = u x v, dim(B) = v x w and dim(C) = u x w.
		It is assumed that A, B, and C are double matrices previously allocated with dmatrix()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		A = dmatrix(1,u,1,v), B = dmatrix(1,v,1,w), C = dmatrix(1,u,1,w)

		CAUTION: NO ERROR TESTING IS PERFORMED. IF THESE MATRICES ARE NOT GIVEN IN THIS WAY,
		NO WARRANTY FOR THE VALIDITY OF THE RESULT NOR FOR THE REMAINS OF YOUR MEMORY...
	*/


	//local variables

	int i, j, k; //counter variables 

	//initialization

	for (i=1; i<= u; i++){
		for (j=1; j<=w; j++){
			C[i][j] = 0.0;
		}
	}




	// **********************************
	// * Perform simple matrix multiply *
	// **********************************
    for(i=1;i<=u;i++) {			// 1 to rows(A)
		for(j=1;j<=w;j++) {		// 1 to cols(B)
			for(k=1;k<=v;k++) { // 1 to cols(A)
			
				C[i][j] = C[i][j] + A[i][k]*B[k][j];
			}
        }
    }


}

/**************************************************************************************/

void dmat_times_vec (
	double **A,
	double *b,
	double *c,
	int u,
	int v
	)
{
	/*	This routine computes c = A*b, with dim(A) = u x v, dim(b) = v x 1 and dim(c) = u x 1.
		It is assumed that A is a double matrix previously allocated with dmatrix()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		A = dmatrix(1,u,1,v)
		and b and c are double vectors previously allocated with dvector() in form of
		b = dvector(1,v), c = dvector(1,u).

		CAUTION: NO ERROR TESTING IS PERFORMED. IF THESE MATRICES ARE NOT GIVEN IN THIS WAY,
		NO WARRANTY FOR THE VALIDITY OF THE RESULT NOR FOR THE REMAINS OF YOUR MEMORY...
	*/


	//local variables

	int i, j; //counter variables 


	// initialization
	for (i=1; i<=u; i++){
		c[i] = 0.0;
	}

	
	// **********************************
	// * Perform simple matrix multiply *
	// **********************************
    for(i=1;i<=u;i++) {			// 1 to rows(A)
		for(j=1;j<=v;j++) {		// 1 to cols(A)
			
			c[i] = c[i] + A[i][j]*b[j];
			
        }
    }


}

/**************************************************************************************/

void dmat_transp(
	double **in_mat,
	double **out_transp_mat,
	int u,
	int v
	)
{
	/*	This function transposes matrix 'in_mat' and saves the result in matrix
		'out_transp_mat'.
		
		It is assumed that 'in_mat' and 'out_transp_mat' are matrices previously 
		allocated with dmatrix() from nrutil.c (see: Numerical Recipes in C), in form of 
		in_mat = dmatrix(1,u,1,v), out_transp_mat = dmatrix(1,v,1,u)
	*/

	//local variables
	int i,j; //counter variables

	for(i=1; i<=u; i++){
		for(j=1; j<=v; j++){
			out_transp_mat[j][i] = in_mat[i][j];
		}
	}

}
/**************************************************************************************/

void print_matrix(
	double **mat,
	int u,
	int v
	)
{
	/* This function is for debugging purpose. It prints out a matrix with dimensions
	   u times v.
	*/

	int i,j; //counter variables

	for(i=1; i<=u; i++){
		for(j=1; j<=v; j++){
		printf("%f ", mat[j][i]);
		}
		printf("\n");
	}
}
/**************************************************************************************/

void rot_x(
	double **rotmatrix,
	double angle
	)
{
	/*	This function returns the 3 x 3 - rotation matrix of 'angle' around the x-axis
		It is assumed that rotmatrix is a double matrix previously allocated with dmatrix()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		rotmatrix = dmatrix(1,3,1,3)
		'angle' must be given in rad.
		The result will be of the form
			  
						[	1  0           0			]
		rotmatrix =		[	0  cos(angle)  sin(angle)	]
						[	0  -sin(angle) cos(angle)	]
	*/

	rotmatrix[1][1] = 1.0;
	rotmatrix[1][2] = 0.0;
	rotmatrix[1][3] = 0.0;
	
	rotmatrix[2][1] = 0.0;
	rotmatrix[2][2] = cos(angle);
	rotmatrix[2][3] = sin(angle);
	
	rotmatrix[3][1] = 0.0;
	rotmatrix[3][2] = -sin(angle);
	rotmatrix[3][3] = cos(angle);


}

/**************************************************************************************/


void rot_y(
	double **rotmatrix,
	double angle
	)
{
	/*	This function returns the 3 x 3 - rotation matrix of 'angle' around the y-axis
		It is assumed that rotmatrix is a double matrix previously allocated with dmatrix()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		rotmatrix = dmatrix(1,3,1,3)
		'angle' must be given in rad.
		The result will be of the form

						[cos(angle)     0    	-sin(angle)	] 	
		rotmatrix =		[0              1  	    0			]
						[sin(angle)     0  	    cos(angle)	]	  
						
	*/

	rotmatrix[1][1] = cos(angle);
	rotmatrix[1][2] = 0.0;
	rotmatrix[1][3] = -sin(angle);
	
	rotmatrix[2][1] = 0.0;
	rotmatrix[2][2] = 1.0;
	rotmatrix[2][3] = 0.0;
	
	rotmatrix[3][1] = sin(angle);
	rotmatrix[3][2] = 0.0;
	rotmatrix[3][3] = cos(angle);


}

/**************************************************************************************/


void rot_z(
	double **rotmatrix,
	double angle
	)
{
	/*	This function returns the 3 x 3 - rotation matrix of 'angle' around the z-axis
		It is assumed that rotmatrix is a double matrix previously allocated with dmatrix()
		from nrutil.c (see: Numerical Recipes in C), in form of 
		rotmatrix = dmatrix(1,3,1,3)
		'angle' must be given in rad.
		The result will be of the form
	
						[cos(angle)     sin(angle)      0	]
		rotmatrix =		[-sin(angle)    cos(angle)      0	]
						[0				0  	            1	]
	  
						
	*/

	rotmatrix[1][1] = cos(angle);
	rotmatrix[1][2] = sin(angle);
	rotmatrix[1][3] = 0.0;
	
	rotmatrix[2][1] = -sin(angle);
	rotmatrix[2][2] = cos(angle);
	rotmatrix[2][3] = 0.0;
	
	rotmatrix[3][1] = 0.0;
	rotmatrix[3][2] = 0.0;
	rotmatrix[3][3] = 1.0;


}

/**************************************************************************************/


