#ifndef _GESOP_KEPLER
#define _GESOP_KEPLER

/* Definitions */
// better as an enum
#define	MERCURY 1
#define VENUS 2
#define EARTH 3
#define MARS 4
#define JUPITER 5
#define SATURN 6
#define URANUS 7
#define NEPTUNE 8
#define PLUTO 9
#define MOON 10
#define SUN 11


/* Functions */

void init_ephemeris();


void free_ephemeris();


void compute_ephemeris(
	double	julian_date,
	int		TARG,
	int		OBS,
	double	starg[6]
);


void compute_rotation(
	double	julian_date,
	int		TARG,
	double	**rotmatrix
	);


#endif