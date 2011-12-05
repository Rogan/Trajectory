#ifndef _PROBLEM_H_
#define _PROBLEM_H_

//#define DEBUG_MODE
#define PERT_MODE 2				// 1=restricted 3 body, 2=sun, 3=jupiter,venus,mars,oblateness,solarpressure
#define EPOCH 2456658.50000		// 00:00.00 1-Jan-2014 http://aa.usno.navy.mil/data/docs/JulianDate.php
#define TERMINATION_MODE 1		//1=pECI, 2=aECI, 3=epsilonECI, 4=rLCI

#define NO_CONNECT 0
#define SOFT_CONNECT 1
#define HARD_CONNECT 2

#define PIECEWISE_CONSTANT 0
#define PIECEWISE_LINEAR 1
#define PIECEWISE_CUBIC_SPLINE 2

/* problem.c */
const double aref;
const double emmisivity;
const double reflectivity;
const double area_efficiency;
const double conversion_efficiency;
const double thrust_arc;
const double thrust_ppt;
const double ce_arc;
const double ce_ppt;
const double power_arc;
const double power_ppt;

#endif