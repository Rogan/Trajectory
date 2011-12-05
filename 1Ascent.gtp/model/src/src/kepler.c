/*	This file contains the functions

	init_ephemeris()
	compute_ephemeris()
	free_ephemeris()
	compute_rotation()

*/

#define ABCORR	"NONE"
#define FRAME	"J2000"
#define SPK		"C:/Users/shimmin/Documents/bw1trajectory/STK/3.1CaptureTimeLCI.gtp/model/src/ephemeris/de421.bsp"							// ephemeris
#define LSK		"C:/Users/shimmin/Documents/bw1trajectory/STK/3.1CaptureTimeLCI.gtp/model/src/ephemeris/naif0009.tls"						// time (leap seconds)
#define ECF		"C:/Users/shimmin/Documents/bw1trajectory/STK/3.1CaptureTimeLCI.gtp/model/src/ephemeris/earth_assoc_itrf93.tf"				// Earth fixed reference frame
#define LCF		"C:/Users/shimmin/Documents/bw1trajectory/STK/3.1CaptureTimeLCI.gtp/model/src/ephemeris/moon_080317.tf"						// Moon fixed reference frame
#define LBPC	"C:/Users/shimmin/Documents/bw1trajectory/STK/3.1CaptureTimeLCI.gtp/model/src/ephemeris/moon_pa_de421_1900-2050.bpc"		// Moon orientation
#define EBPC	"C:/Users/shimmin/Documents/bw1trajectory/STK/3.1CaptureTimeLCI.gtp/model/src/ephemeris/earth_070425_370426_predict.bpc"	// Earth orientation
#define TPC		"C:/Users/shimmin/Documents/bw1trajectory/STK/3.1CaptureTimeLCI.gtp/model/src/ephemeris/pck00009.tpc"						// planetary size and shape
#define	LCFA	"C:/Users/shimmin/Documents/bw1trajectory/STK/3.1CaptureTimeLCI.gtp/model/src/ephemeris/moon_assoc_me.tf"					// makes MOON_ME default fixed frame for the Moon

#include <stdio.h>
#include "kepler.h"
#include "constants.h"
#include "../cspice/include/SpiceUsr.h"

/* specifically, uses
		furnsh_c.c
		spkez_c.c
		unitim_c.c
		pxform_c.c
		unload_c.c
	these need to be included in the compiler directive

	Edit: easier to include \src\cspice\lib\*.lib in the linker directive
*/

/**************************************************************************************/


void init_ephemeris(int *error)
{
	*error = 0;

	furnsh_c ( SPK );	// Ephemeris
	furnsh_c ( LSK );	// Time conversion (light seconds)
	furnsh_c ( LCF );	// Lunar centric frame
	furnsh_c ( ECF );	// Earth centric frame
	furnsh_c ( EBPC );	// Earth orientation
	furnsh_c ( LBPC );	// Moon orientation
	furnsh_c ( LCFA );	// Set default fixed lunar frame

}




void free_ephemeris()
{
	unload_c ( SPK );
	unload_c ( LSK );
	unload_c ( LCF );
	unload_c ( ECF );
	unload_c ( EBPC );
	unload_c ( LBPC );
	unload_c ( LCFA );
}



void compute_ephemeris(
	double	julian_date,
	int		TARG,
	int		OBS,
	double	starg[6]
	)
{
	double lt = 0.0;
	double ET;

	ET = unitim_c(julian_date, "JED", "ET");

	spkez_c(naif_id[TARG], ET, FRAME, ABCORR, naif_id[OBS], starg, &lt);

}


void compute_rotation(
	double	julian_date,
	int		center,
	double	**rotmatrix
	)
{
	double	ET;
	double	temp[3][3];
	int		i,j;

	ET = unitim_c(julian_date, "JED", "ET");

	if(center == EARTH)
		pxform_c ( FRAME, "ITRF93", ET, temp );
	else // center == MOON
		pxform_c ( FRAME, "MOON_ME", ET, temp );

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			rotmatrix[i+1][j+1] = temp[i][j];
}
