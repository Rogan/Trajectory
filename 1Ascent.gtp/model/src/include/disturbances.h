#ifndef _DISTURBANCES_H_
#define _DISTURBANCES_H_

/* disturbances.c */
void	ir_itheta_ih_transp(double *,double **	);
void	Fqk(double *,double *,double *	);
void	third_body_pert_cart(double ,int ,int ,double *,double * );
void	AccelEarthObl(double *r_vec, double *out_pert_vec);
void	AccelHarmonic (double *, double **, double ,double ,const double [21][21],int ,int ,double *);
void	solar_rad_pressure_cart(double ,double ,double *, int, double *);
double	shadowfunc(double *, double *, double);
double	sunangle(double *, const double *);
void	AccelEarthObl2(double *r_vec, double *out_pert_vec);

#endif