/*******************************************************************************
 Global constants, defined for the problem:
*******************************************************************************/

#include "problem.h"

const int moonmod = 1;

// PPT data from Matthias Lau, arcjet data from Birk Wollenhaupt
const double thrust_arc = 102.5e-3;		/* Arcjet thrust [N] */
const double thrust_ppt = 4.88e-3;		/* PPT thrust [N] */
const double ce_arc = 4.768e3;			/* Arcjet exhaust velocity [m/s] */
const double ce_ppt = 27.0e3;			/* PPT exhaust velocity [m/s] */
const double power_arc = 801.0;			/* Arcjet power consumption [W] */
const double power_ppt = 204.0;			/* PPT power consumption [W] */

const double aref = 5.4;				/* Surface area of spacecraft [m²] after compensating for 45° central sections */
const double emissivity = 0.85;			/* Emissivity of spacecraft [-] */
const double reflectivity = 0.3;		/* Reflectivity of spacecraft [-] */

// Solar panel data from Oliver Zeile
const double conversion_efficiency = 0.27;			/* Energy conversion efficiency of solar panels [-] */
const double area_efficiency = 0.8;					/* Surface coverage of solar panels [-] */
