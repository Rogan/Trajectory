#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#define PI 3.1415926535897932384626433832795	/* should be in math.h, but doesn't work */

/* C.F. Yoder, Astronomic and Geodetic Properties of Earth and the Solar System */
#define R_SUN 6.960e+8							/* radius of the Sun [m] */
#define SPEEDOFLIGHT 299792458.0				/* speed of light [m/s] */
#define SOLARLUMINOSITY 3.846e26				/* solar luminosity [W] (equivalent to total output power of the Sun) */
#define SOLARINTENSITY 1367.56					/* [W/m²] solarluminosity/(4*pi*a_sun^2) */

#define A_SUN 1.49598261e11						/* Semimajor axis of the Earth's orbit around the Sun */
#define GMST2000 1.74476716333061				/* Greenwich Mean Sidereal Time at J2000.0 [rad] */
#define J2000 2451545.0							/* Julian date at J2000 epoch */
#define OMEGA_EARTH	7.2921158553e-5				/* Rotational speed of Earth [rad/s] (p.60 Montenbruck & Gill) */

#define D2R 1.7453292519943295769236907684886e-2	/* Conversion from degrees to radians */
#define R2D 5.7295779513082320876798154814105e+1	/* Conversion from radians to degrees */

const double AU;
const double gmue[12];
const int naif_id[12];

const double R_LP100J;
const double GM_LP100J;
const double CS_LP100J[21][21];

const double R_JGM3;
const double GM_JGM3;
const double CS_JGM3[21][21];   

#endif