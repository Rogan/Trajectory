/* -------------------------------------------------------------------------------
 * Gesop_Model.Auxiliary_Evaluations
 * -------------------------------------------------------------------------------
 *  Copyright
 *                        University of Stuttgart
 *                   Institute of Space Systems (IRS)
 *                          Pfaffenwaldring 31
 *                       70569 Stuttgart, Germany
 *
 *                        FAX: (+49)-711-685-63596
 *				     E-Mail: roeser@irs.uni-stuttgart.de
 * ---------------------------------------------------------------------------
 * Copyright (c) 2010, University of Stuttgart, IRS
 * All rights reserved.
 * ---------------------------------------------------------------------------
 * THE UNIVERSITY OF STUTTGART, IRS AND ASTOS GmbH DISCLAIM ALL WARRANTIES 
 * WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF 
 * MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THEY BE LIABLE FOR ANY 
 * SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER 
 * RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF 
 * CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN 
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 * ---------------------------------------------------------------------------
 * PROJECT:   BW-1
 * AUTHORS:   D. Fischer, C. Moellman, R. Shimmin
 * ---------------------------------------------------------------------------
 */

#include <stdio.h>
#include "gesop_model.h"
#include "kepler.h"
#include "vectors.h"
#include "matrices.h"
#include "nrutil.h"
#include "problem.h"
#include "cart2kep.h"
#include "mee2cart.h"
#include "disturbances.h"
#include "kepler-mee.h"
#include "limanglebetween.h"
#include "constants.h"

/*--------------
 *   Constants -
 *--------------*/
static DLL_Group_Item Algebraic_Function_Names[90] =
	{
	{"AUX_ECI_x",				-1, "X in ECI",											-1, "Meter",			-1,"",-1	},	//0
	{"AUX_ECI_y",				-1, "Y in ECI",											-1, "Meter",			-1,"",-1	},	//1
	{"AUX_ECI_z",				-1, "Z in ECI",											-1, "Meter",			-1,"",-1	},	//2
	{"AUX_ECI_x_dot",			-1, "X dot in ECI",										-1, "Meter/Second",		-1,"",-1	},	//3
	{"AUX_ECI_y_dot",			-1, "Y dot in ECI",										-1, "Meter/Second",		-1,"",-1	},	//4
	{"AUX_ECI_z_dot",			-1, "Z dot in ECI",										-1, "Meter/Second",		-1,"",-1	},	//5
	{"AUX_dist_to_Earth",		-1, "Satellite's distance to Earth",					-1, "Meter",			-1,"",-1	},	//6
	{"AUX_dist_to_Moon",		-1, "Satellite's distance to Moon",						-1, "Meter",			-1,"",-1	},	//7
	{"AUX_v_e",					-1, "Absolute velocity wrt. Earth",						-1, "Meter/Second",		-1,"",-1	},	//8
	{"AUX_v_m",					-1, "Absolute velocity wrt. Moon",						-1, "Meter/Second",		-1,"",-1	},	//9
	{"AUX_epsilon_e",			-1, "Orbital energy wrt. Earth",						-1, "Meter**2/Second**2",-1,"",-1	},	//10
	{"AUX_epsilon_m",			-1, "Orbital energy wrt. Moon",							-1, "Meter**2/Second**2",-1,"",-1	},	//11
	{"AUX_gamma",				-1, "Angle between velocity and thrust vector",			-1, "Radian",			-1,"",-1	},	//12
	{"AUX_ECI_kep_a",			-1, "Semimajor axis",									-1, "Meter",			-1,"",-1	},	//13
	{"AUX_ECI_kep_e",			-1, "Eccentricity",										-1, "None",				-1,"",-1	},	//14
	{"AUX_ECI_kep_i",			-1, "Inclination",										-1, "Radian",			-1,"",-1	},	//15
	{"AUX_ECI_kep_omega",		-1, "Argument of periapsis",							-1, "Radian",			-1,"",-1	},	//16
	{"AUX_ECI_kep_RAAN",		-1, "Right ascension of the ascending node",			-1, "Radian",			-1,"",-1	},	//17
	{"AUX_ECI_kep_nu",			-1, "True anomaly",										-1, "Radian",			-1,"",-1	},	//18
	{"AUX_SEL_x",				-1, "X in SEL",											-1, "Meter",			-1,"",-1	},	//19
	{"AUX_SEL_y",				-1, "Y in SEL",											-1, "Meter",			-1,"",-1	},	//20
	{"AUX_SEL_z",				-1, "Z in SEL",											-1, "Meter",			-1,"",-1	},	//21
	{"AUX_SEL_x_dot",			-1, "X dot in SEL",										-1, "Meter/Second",		-1,"",-1	},	//22
	{"AUX_SEL_y_dot",			-1, "Y dot in SEL",										-1, "Meter/Second",		-1,"",-1	},	//23
	{"AUX_SEL_z_dot",			-1, "Z dot in SEL",										-1, "Meter/Second",		-1,"",-1	},	//24
	{"AUX_LCI_x",				-1, "X in LCI",											-1, "Meter",			-1,"",-1	},	//25
	{"AUX_LCI_y",				-1, "Y in LCI",											-1, "Meter",			-1,"",-1	},	//26
	{"AUX_LCI_z",				-1, "Z in LCI",											-1, "Meter",			-1,"",-1	},	//27
	{"AUX_LCI_x_dot",			-1, "X dot in LCI",										-1, "Meter/Second",		-1,"",-1	},	//28
	{"AUX_LCI_y_dot",			-1, "Y dot in LCI",										-1, "Meter/Second",		-1,"",-1	},	//29
	{"AUX_LCI_z_dot",			-1, "Z dot in LCI",										-1, "Meter/Second",		-1,"",-1	},	//30
	{"AUX_LCI_kep_a",			-1, "Semimajor axis in LCI",							-1, "Meter",			-1,"",-1	},	//31
	{"AUX_LCI_kep_e",			-1, "Eccentricity in LCI",								-1, "None",				-1,"",-1	},	//32
	{"AUX_SEL_kep_i",			-1, "Inclination in SEL",								-1, "Radian",			-1,"",-1	},	//33
	{"AUX_SEL_kep_omega",		-1, "Argument of periapsis in SEL",						-1, "Radian",			-1,"",-1	},	//34
	{"AUX_SEL_kep_RAAN",		-1, "RAAN in SEL",										-1, "Radian",			-1,"",-1	},	//35
	{"AUX_SEL_kep_nu",			-1, "true anomaly in SEL",								-1, "Radian",			-1,"",-1	},	//36
	{"AUX_LCI_kep_i",			-1, "Inclination in LCI",								-1, "Radian",			-1,"",-1	},	//37
	{"AUX_LCI_kep_omega",		-1, "Argument of periapsis in LCI",						-1, "Radian",			-1,"",-1	},	//38
	{"AUX_LCI_kep_RAAN",		-1, "RAAN in LCI",										-1, "Radian",			-1,"",-1	},	//39
	{"AUX_Moon_x",				-1, "Moon X",											-1, "Meter",			-1,"",-1	},	//40
	{"AUX_Moon_y",				-1, "Moon Y",											-1, "Meter",			-1,"",-1	},	//41
	{"AUX_Moon_z",				-1, "Moon Z",											-1, "Meter",			-1,"",-1	},	//42
	{"AUX_dist_Earth_Moon",		-1, "Distance from Earth to Moon",						-1, "Meter",			-1,"",-1	},	//43
	{"AUX_Sun_x",				-1, "Sun X",											-1, "Meter",			-1,"",-1	},	//44
	{"AUX_Sun_y",				-1, "Sun Y",											-1, "Meter",			-1,"",-1	},	//45
	{"AUX_Sun_z",				-1, "Sun Z",											-1, "Meter",			-1,"",-1	},	//46
	{"AUX_dist_Earth_Sun",		-1, "Distance from Earth to Sun",						-1, "Meter",			-1,"",-1	},	//47
	{"AUX_Earth_x",				-1, "Earth X in SEL",									-1, "Meter",			-1,"",-1	},	//48
	{"AUX_Earth_y",				-1, "Earth Y in SEL",									-1, "Meter",			-1,"",-1	},	//49
	{"AUX_Earth_z",				-1, "Earth Z in SEL",									-1, "Meter",			-1,"",-1	},	//50

	{"AUX_L1_x",				-1, "L1 X in sel. coords",								-1, "Meter",			-1,"",-1	},	//51 obsolete
	{"AUX_L1_y",				-1, "L1 Y in sel. coords",								-1, "Meter",			-1,"",-1	},	//52 obsolete
	{"AUX_L1_z",				-1, "L1 Z in sel. coords",								-1, "Meter",			-1,"",-1	},	//53 obsolete
	{"AUX_dist_to_L1",			-1, "distance from satellite to L1",					-1, "Meter",			-1,"",-1	},	//54 obsolete

	{"AUX_jacobi_integral",		-1, "value of jacobi integal",							-1, "None",				-1,"",-1	},	//55 obsolete
	{"AUX_i_MOP",				-1, "inclination in moon orbit plane frame",			-1, "Radian",			-1,"",-1	},	//56 obsolete
	{"AUX_omega_MOP",			-1, "argument of periapsis in moon orbit plane frame",	-1, "Radian",			-1,"",-1	},	//57 obsolete
	{"AUX_cap_omega_MOP",		-1, "longitude of asc. node in moon orbit plane frame",	-1, "Radian",			-1,"",-1	},	//58 obsolete

	{"AUX_phase",				-1, "Phase between satellite and moon",					-1, "Radian",			-1,"",-1	},	//59 
	{"AUX_thrust",				-1,	"Thrust level",										-1, "Newton",			-1,"",-1	},	//60
	{"AUX_acc_thrust",			-1,	"Thrust acc.",										-1, "Meter/Second**2",	-1,"",-1	},	//61
	{"AUX_acc_thirdbody_earth",	-1,	"TB Earth acc.",									-1, "Meter/Second**2",	-1,"",-1	},	//62
	{"AUX_acc_thirdbody_moon",	-1,	"TB Moon acc.",										-1, "Meter/Second**2",	-1,"",-1	},	//63
	{"AUX_acc_thirdbody_sun",	-1,	"TB Sun acc.",										-1, "Meter/Second**2",	-1,"",-1	},	//64
	{"AUX_acc_thirdbody_jupiter",	-1,	"TB Jupiter acc.",								-1, "Meter/Second**2",	-1,"",-1	},	//65
	{"AUX_acc_thirdbody_venus",	-1,	"TB Venus acc.",									-1, "Meter/Second**2",	-1,"",-1	},	//66
	{"AUX_acc_thirdbody_mars",	-1,	"TB Mars acc.",										-1, "Meter/Second**2",	-1,"",-1	},	//67
	{"AUX_acc_grav_harm",		-1,	"Acc. due to gravity harmonics",					-1, "Meter/Second**2",	-1,"",-1	},	//68
	{"AUX_acc_solp",			-1,	"Acc. due to solar rad. pressure",					-1, "Meter/Second**2",	-1,"",-1	},	//69
	{"AUX_Moon_a_ECI",			-1,	"Lunar semimajor axis in ECI",						-1, "Meter",			-1,"",-1	},	//70
	{"AUX_Moon_e_ECI",			-1,	"Lunar eccentricity in ECI",						-1, "None",				-1,"",-1	},	//71
	{"AUX_Moon_i_ECI",			-1,	"Lunar inclination in ECI",							-1, "Radian",			-1,"",-1	},	//72
	{"AUX_Moon_o_ECI",			-1,	"Lunar longitude of ascending node in ECI",			-1, "Radian",			-1,"",-1	},	//73
	{"AUX_Moon_O_ECI",			-1,	"Lunar argument of periapsis in ECI",				-1, "Radian",			-1,"",-1	},	//74
	{"AUX_Moon_nu_ECI",			-1,	"Lunar mean anomaly in ECI",						-1, "Radian",			-1,"",-1	},	//75
	{"AUX_delta_V",				-1,	"Delta V",											-1, "Meter/Second",		-1,"",-1	},	//76
	{"AUX_periapsis",			-1,	"Periapsis distance in ECI",						-1, "Meter",			-1,"",-1	},	//77
	{"FBC_periapsis",			-1,	"Terminal condition: Periapsis distance in ECI",	-1, "None",				-1,"",-1	},	//78 obsolete
	{"IBC_a_max",				-1,	"Initial condition: Semimajor axis maximum",		-1, "None",				-1,"",-1	},	//79 obsolete
	{"IBC_a_min",				-1,	"Initial condition: Semimajor axis minimum",		-1, "None",				-1,"",-1	},	//80 obsolete
	{"IBC_e",					-1,	"Initial condition: Eccentricity",					-1, "None",				-1,"",-1	},	//81 obsolete
	{"IBC_i",					-1,	"Initial condition: Inclination",					-1, "None",				-1,"",-1	},	//82 obsolete
	{"IBC_mass",				-1,	"Initial condition: Mass",							-1, "None",				-1,"",-1	},	//83 obsolete
	{"AUX_energy_used",			-1,	"Energy used by thrusters",							-1, "Joule",			-1,"",-1	},	//84
	{"AUX_energy_generated",	-1,	"Energy generated by solar panels",					-1, "Joule",			-1,"",-1	},	//85
	{"AUX_oblateness_higher",	-1,	"Higher order Earth oblateness",					-1, "Meter/Second**2",	-1,"",-1	},	//86
	{"AUX_energy_integral",		-1,	"Integral of instantaneous power",					-1, "Joule",			-1,"",-1	},	//87
	{"AUX_energy_net",			-1,	"Difference of thruster and PV integral energies",	-1, "Joule",			-1,"",-1	},	//88
	{"AUX_oblateness_J2",		-1,	"Earth oblateness due to J2 harmonic",				-1, "Meter/Second**2",	-1,"",-1	},	//89
};

/*---------------
 *   Operations -
 *---------------*/
void __cdecl Algebraic_Functions_Info (
	const int		phase,				// Current phase number
	int				*item_size,			// Number of algebraic output functions
	DLL_Group_Item	**item_vec,			// Vector of function descriptions
	int				*error)				// Error flag
	/* SYNOPSIS:
	 * Set header vector (possibly specific to each phase).
	 * Take care that names are unique within each phase, taking into
	 * account also the names given in the TOPS file! */
{
	int i;

#ifdef DEBUG_MODE
	printf("Algebraic Functions Info\n");_flushall();
#endif

	*error = 0;
	//*item_size = 0; 
	*item_size = sizeof(Algebraic_Function_Names)/sizeof(DLL_Group_Item);

	for (i=0; i<*item_size; i++) {
		Algebraic_Function_Names[i].Name_len = strlen(Algebraic_Function_Names[i].Name);
		Algebraic_Function_Names[i].Desc_len = strlen(Algebraic_Function_Names[i].Desc);
		Algebraic_Function_Names[i].Unit_len = strlen(Algebraic_Function_Names[i].Unit);
		Algebraic_Function_Names[i].Group_len = strlen(Algebraic_Function_Names[i].Group);
	}
	*item_vec = Algebraic_Function_Names;
}

void __cdecl Algebraic_Functions (
	const int				*phase,		// Current phase number
	const int				*dimx,	    // Dimension of state vector X at time t
	const int				*dimu,	    // Dimension of control vector U at time t
	const int				*dimip,	    // Dimension of integer parameter vector of phase
	const int				*dimrp,	    // Dimension of real parameter vector of phase
	const Phase_Info_Type	*fazinf,	// Record with additional phase info
	const double			*t,			// Time t of evaluation
	const double			x[],		// State vector X at time t
	const double			u[],		// Control vector U at time t
	const double			udot[],		// Control derivative vector at t
	const int				ipar[],		// Integer parameter vector of phase
	const double			rpar[],		// Real parameter vector of phase
	const int				*dimof,		// Number of algebraic output functions
	double					*outfun,	// Values of algebraic output functions
	int						*error)		// Flag for error handling

	/* computes the algebraic "output" functions at the specified
	 * evaluation time of the given phase, as a function of the
	 * phase specific real parameter vector and the state, control
	 * and control derivative vector at the evaluation time. */

{	
	/* local variables */
	int		i,j;			//counter variables

	double	p, f, g, h, k, L;

	double	*mee_coords_ECI;				// contains the first six elements of the state, ie the equ. elements
	double	*mee_coords_LCI;				// contains the first six elements of the state in LCI frame
	double	*kep_coords_ECI;				// contains the first six elements of the state in keplerian coords.
	double	*kep_coords_LCI;				// contains kepl. elements of satellite in LCI frame
	double	*kep_coords_SEL;				// contains kepl. elements in selenographic frame
	double	*cart_inert_coords_ECI;			// contains the first six elements in cartesian ECI coords
	double	*cart_inert_coords_LCI;			// contains the first six elements in cartesian LCI coords
	double	*cart_inert_coords_SEL;			// contains the first six elements in cartesian selenographic coords
											// the sel. frame is not inertial any more, as it is rotated!
	double	*r_vec_ECI, *v_vec_ECI;			// contains radius vector and velocity vector in inertial cartesian coords ECI in m and m/s
	double	*r_vec_LCI, *v_vec_LCI;			// contains radius vector and velocity vector in inertial cartesian coords LCI in m and m/s
	double	*r_vec, *v_vec;					// contains radius vector and velocity vector of satellite (points to either ECI or LCI vectors depending on phase)
	double	*v_vec_norm;					// unit vector in velocity direction
	double	*v_vec_norm_rotrad;				// unit vector in velocity direction converted to rotating radial (body fixed) frame
	double	r_ECI;							// magnitude of radius vector in ECI
	double	v_ECI;							// magnitude of velocity vector in ECI
	double	r_LCI;							// magnitude of radius vector in LCI
	double	v_LCI;							// magnitude of velocity vector in LCI

	int		center, target;					// these are used for all calls to compute_ephemeris()
	double	julian_date;
	double	r_norm;							// [m] magnitude of position vector
	double	v_norm;							// [m/s] magnitude of velocity vector

	double	ephemerides_moon[6];
	double	ephemerides_sun[6];
	double	*earth_coords_LCI;				// position vector of earth in LCI [m]
	double	*earth_coords_SEL;				// position vector of earth in sel. frame [m]
	double	*moon_cart_coords_ECI;
	double	*moon_coords_kep_ECI;
	double	*sun_cart_coords_ECI;
	double	**Q_transp;						// contains the rotation matrix for the rotating radial frame in which 
											// disturbing accelerations are expressed (x_rot = Q_transp * x_cart_inert)
	double	**inert2selen;					// contains the rotation matrix from inertial to selenographic coords
	double	**inert2geo;					// contains the rotation matrix from inertial to geographic coords

	double	gamma, dot_prod;

	double	phase_moon;						// phase of sat. wrt. Moon

	double	thrust;							// [N]
	double	ce;								// [m/s] exhaust velocity
	double	g_mue;							// [m²/s³] gravitational constant of central body
	double	thrust_mag;						// [-] normalised

	// thirdbody perturbations
	double	*pertearth;						// cartesian perturbing acceleration due to Earth (inertial)
	double	*pertsun;						// cartesian perturbing acceleration due to Sun (inertial)
	double	*pertmoon;						// cartesian perturbing acceleration due to Moon
	double	*pertjupiter;					// cartesian perturbing acceleration due to jupiter
	double	*pertvenus;						// cartesian perturbing acceleration due to venus
	double	*pertmars;						// cartesian perturbing acceleration due to mars
	
	double	*thirdbody_acc_cart, *thirdbody_acc_rotrad; // third body perturbing accelerations in cartesian and in rotating reference frame	

	// gravity harmonic perturbations
	double	*pertharm_cart;					// inertial cartesian perturbing acceleration due to J2 and J3 gravity harmonic
	double	*pertharm_rotrad;				// ... in rotating reference frame

	double	*pertharm_cart_CS;				// inertial cartesian perturbing acceleration due to full degree/order 20 harmonic
	double	*pertharm_rotrad_CS;			// ... in rotating reference frame

	double	*pertharm_cart2;				// inertial cartesian perturbing acceleration due to J2 gravity harmonic
	double	*pertharm_rotrad2;				// ... in rotating reference frame

	// solar pressure perturbations
	double	*pertsolp_cart;					//inertial cartesian perturbing acceleration due to solar pressure
	double	*pertsolp_rotrad;				// ... in rotating reference frame

	double	thruster_power;
	double	solar_power;

	double	*craft2sun;						// intertial cartesian position vector of sun relative to craft

	// power calculations
	double	battery_level;					// [J]
	double	net_power;						// [W]

	// time integrals
	static double	deltaV = 0;				// [m/s]
	static double	old_t = 0;				// [d]
	static double	thruster_energy = 0;	// [J] 
	static double	solar_energy = 0;		// [J]
	static double	net_energy = 2.8*3.6e6;	// [J]

#ifdef DEBUG_MODE
printf("Algebraic Functions\n");_flushall();
#endif


	*error = 0;

	Q_transp				= dmatrix(1,3,1,3);
	inert2selen				= dmatrix(1,3,1,3);
	inert2geo				= dmatrix(1,3,1,3);

	mee_coords_ECI			= dvector(1,6);
	mee_coords_LCI			= dvector(1,6);
	kep_coords_ECI			= dvector(1,6);	
	kep_coords_LCI			= dvector(1,6);
	kep_coords_SEL			= dvector(1,6);
	cart_inert_coords_ECI	= dvector(1,6);
	cart_inert_coords_LCI	= dvector(1,6);
	cart_inert_coords_SEL	= dvector(1,6);
	earth_coords_SEL		= dvector(1,3);
	earth_coords_LCI		= dvector(1,3);
	moon_cart_coords_ECI	= dvector(1,6);
	moon_coords_kep_ECI		= dvector(1,6);
	sun_cart_coords_ECI		= dvector(1,6);
	r_vec_ECI				= dvector(1,3);
	v_vec_ECI				= dvector(1,3);
	r_vec_LCI				= dvector(1,3);
	v_vec_LCI				= dvector(1,3);
	v_vec_norm				= dvector(1,3);
	v_vec_norm_rotrad		= dvector(1,3);

	craft2sun				= dvector(1,3);
	
	pertearth				= dvector(1,3);
	pertsun					= dvector(1,3);
	pertmoon				= dvector(1,3);
	pertjupiter				= dvector(1,3);
	pertvenus				= dvector(1,3);
	pertmars				= dvector(1,3);

	thirdbody_acc_cart		= dvector(1,3);
	thirdbody_acc_rotrad	= dvector(1,3);

	pertharm_cart			= dvector(1,3);
	pertharm_rotrad			= dvector(1,3);
	
	pertharm_cart_CS		= dvector(1,3);
	pertharm_rotrad_CS		= dvector(1,3);
	
	pertharm_cart2			= dvector(1,3);
	pertharm_rotrad2		= dvector(1,3);

	pertsolp_cart			= dvector(1,3);
	pertsolp_rotrad			= dvector(1,3);	


	for (i=1; i<=3; i++){
		for (j=1; j<=3; j++){
			Q_transp[i][j]			= 0.0;
			inert2selen[i][j]		= 0.0;
			inert2geo[i][j]			= 0.0;
		}
		r_vec_ECI[i]				= 0.0;
		v_vec_ECI[i]				= 0.0;
		r_vec_LCI[i]				= 0.0;
		v_vec_LCI[i]				= 0.0;
		v_vec_norm[i]				= 0.0;
		v_vec_norm_rotrad[i]		= 0.0;

		earth_coords_SEL[i]			= 0.0;
		earth_coords_LCI[i]			= 0.0;
		
		pertearth[i]				= 0.0;
		pertsun[i]					= 0.0;
		pertmoon[i]					= 0.0;
		pertjupiter[i]				= 0.0;
		pertvenus[i]				= 0.0;
		pertmars[i]					= 0.0;

		thirdbody_acc_cart[i]		= 0.0;
		thirdbody_acc_rotrad[i]		= 0.0;

		pertharm_cart[i]			= 0.0;
		pertharm_rotrad[i]			= 0.0;

		pertharm_cart2[i]			= 0.0;
		pertharm_rotrad2[i]			= 0.0;

		pertharm_cart_CS[i]			= 0.0;
		pertharm_rotrad_CS[i]		= 0.0;
		
		pertsolp_cart[i]			= 0.0;
		pertsolp_rotrad[i]			= 0.0;
	}
	
	for (i=1; i<=6; i++) {
		mee_coords_ECI[i]			= 0.0;
		mee_coords_LCI[i]			= 0.0;
		kep_coords_ECI[i]			= 0.0;
		kep_coords_LCI[i]			= 0.0;
		kep_coords_SEL[i]			= 0.0;
		cart_inert_coords_ECI[i]	= 0.0;
		cart_inert_coords_LCI[i]	= 0.0;
		cart_inert_coords_SEL[i]	= 0.0;
		ephemerides_moon[i-1]		= 0.0;
		ephemerides_sun[i-1]		= 0.0;
		moon_cart_coords_ECI[i]		= 0.0;
		moon_coords_kep_ECI[i]		= 0.0;
		sun_cart_coords_ECI[i]		= 0.0;
	}

	/* default initialization  */
	for (i=0; i<*dimof; i++) {
		outfun[i] = 0.0;
	}
	
		
	p = x[0];
	f = x[1];
	g = x[2];
	h = x[3];
	k = x[4];
	L = x[5];
	
	thrust_mag = u[3];
	

	/************ compute ephemerides (independent of spacecraft position) ************/

	julian_date = EPOCH + rpar[0] + x[6];	// = reference time + start time + elapsed mission time [days] = current time
	center = EARTH;
	target = MOON;
	compute_ephemeris(julian_date, target, center, ephemerides_moon);
	target = SUN;
	compute_ephemeris(julian_date, target, center, ephemerides_sun);
	
	for (i=1;i<=6;i++) {
		moon_cart_coords_ECI[i] = ephemerides_moon[i-1]*1.0e3;
		sun_cart_coords_ECI[i] = ephemerides_sun[i-1]*1.0e3;
	}
	
	// compute kepl. elements of moon orbit
	cart2kep(moon_cart_coords_ECI, gmue[EARTH], moon_coords_kep_ECI);
	
	compute_rotation(julian_date, MOON, inert2selen);
	compute_rotation(julian_date, EARTH, inert2geo);

	// compute earth coords in LCI
	for (i=1;i<=3;i++) {
		earth_coords_LCI[i] = -ephemerides_moon[i-1]*1.0e3;
	}
	dmat_times_vec(inert2selen,earth_coords_LCI,earth_coords_SEL,3,3);


	/************ compute spacecraft position in different frames ************/

/*	if(*phase < 3) {	// ECI

		for (i=1; i<=6; i++) {
			mee_coords_ECI[i] = x[i-1];
		}

		mee2cart(mee_coords_ECI, gmue[EARTH], cart_inert_coords_ECI); 

		for (i=1; i<=6; i++) {
			cart_inert_coords_LCI[i] = cart_inert_coords_ECI[i] - moon_cart_coords_ECI[i];
		}

		// compute rotating radial coordinate frame
		ir_itheta_ih_transp(cart_inert_coords_ECI, Q_transp); //obtain Q_transp as [ir itheta ih]'

	} else {	// LCI
*/
		for (i=1; i<=6; i++) {
			mee_coords_LCI[i] = x[i-1];
		}

		mee2cart(mee_coords_LCI, gmue[MOON], cart_inert_coords_LCI); 

		for (i=1; i<=6; i++) {
			cart_inert_coords_ECI[i] = cart_inert_coords_LCI[i] + moon_cart_coords_ECI[i];
		}
		
		// compute rotating radial coordinate frame
		ir_itheta_ih_transp(cart_inert_coords_LCI, Q_transp); //obtain Q_transp as [ir itheta ih]'

//	}

	/************ additional frames and parameters ************/
	mee2kep(mee_coords_LCI, kep_coords_LCI);	// may be a problem if not in lunar orbit?
	mee2kep(mee_coords_ECI, kep_coords_ECI);
	phase_moon = moon_coords_kep_ECI[4] + moon_coords_kep_ECI[5] + moon_coords_kep_ECI[6] - mee_coords_ECI[6];
	limanglebetween(&phase_moon, -PI, PI);


	// split cartesian coords into radius and velocity-vector
	for (i=1; i<=3; i++){						
		r_vec_ECI[i] = cart_inert_coords_ECI[i];
		v_vec_ECI[i] = cart_inert_coords_ECI[i+3];
		r_vec_LCI[i] = cart_inert_coords_LCI[i];
		v_vec_LCI[i] = cart_inert_coords_LCI[i+3];
	}

	// and determine norms
	r_ECI = dvec_abs(r_vec_ECI,3);
	r_LCI = dvec_abs(r_vec_LCI,3);
	v_ECI = dvec_abs(v_vec_ECI,3);
	v_LCI = dvec_abs(v_vec_LCI,3);

	// set phase-based parameters
	if(*phase < 3) {	// ECI
		r_vec = r_vec_ECI;
		v_vec = v_vec_ECI;
		r_norm = r_ECI;
		v_norm = v_ECI;
	} else {	// LCI
		r_vec = r_vec_LCI;
		v_vec = v_vec_LCI;
		r_norm = r_LCI;
		v_norm = v_LCI;
	}

	// --- calculate direction of thrust vector ---
	for (i=1; i<=3; i++) {
		v_vec_norm[i] = v_vec[i] / v_norm;
	}
	
	// convert norm. v vector to rotating radial frame
	dmat_times_vec(Q_transp, v_vec_norm, v_vec_norm_rotrad, 3, 3);
	
	dot_prod = u[0]*v_vec_norm_rotrad[1]+u[1]*v_vec_norm_rotrad[2]+u[2]*v_vec_norm_rotrad[3];

	// bring dot_prod to allowed range for acos calculation (could be out of allowed range due to numerical errors)
	if (dot_prod > 1.0)
		dot_prod = 1.0;
	else if (dot_prod < -1.0)
		dot_prod = -1.0;
	gamma = acos(dot_prod); // see eq. (5.3) Diplomarbeit S. Erb
		


	// compute coords in SEL if in lunar orbit
	if(kep_coords_LCI[2] < 1.0) {
		dmat_times_vec(inert2selen,cart_inert_coords_LCI,cart_inert_coords_SEL,3,3);
		dmat_times_vec(inert2selen,&cart_inert_coords_LCI[3],&cart_inert_coords_SEL[3],3,3);
		cart2kep(cart_inert_coords_SEL, gmue[MOON], kep_coords_SEL);
	}

	/************ replicate RHS ************/

//	if(*phase < 3) {
//		center = EARTH;
//	else
		center = MOON;
	g_mue = gmue[center];

	// thrust magnitude
	switch(*phase)
	{
	case 1:
		thrust = thrust_mag*thrust_arc;
		ce = ce_arc;
		thruster_power = thrust_mag*power_arc;
//		break;
	case 2:
		thrust = thrust_mag*thrust_ppt;
		ce = ce_ppt;
		thruster_power = thrust_mag*power_ppt;
//		break;
	case 3:
		thrust = thrust_mag*thrust_arc;
		ce = ce_arc;
		thruster_power = thrust_mag*power_arc;
		break;
	case 4:
		thrust = thrust_mag*thrust_ppt;
		ce = ce_ppt;
		thruster_power = thrust_mag*power_ppt;
		break;
	case 5:
		thrust = 0.0;
		ce = 1.0;
		thruster_power = 0.0;
		break;
	default:
		printf("Unrecognised phase number in RHS function\n");
		break;
	} // switch



	// --- compute accelerations due to perturbations --- //

	// compute 3rd body perturbations, unscaled	
	if (center == EARTH) {
		target = MOON;
		third_body_pert_cart(julian_date, center, target, r_vec, pertmoon);
	}
	else {	
		target = EARTH;
		third_body_pert_cart(julian_date, center, target, r_vec, pertearth);
	}

	target = SUN;
	third_body_pert_cart(julian_date, center, target, r_vec, pertsun);
	
	target = JUPITER;
	third_body_pert_cart(julian_date, center, target, r_vec, pertjupiter);
	
	target = VENUS;
	third_body_pert_cart(julian_date, center, target, r_vec, pertvenus);
	
	target = MARS;
	third_body_pert_cart(julian_date, center, target, r_vec, pertmars);
	
	for (i=1; i<=3; i++){
		thirdbody_acc_cart[i] = pertearth[i] + pertsun[i] + pertmoon[i] + pertjupiter[i] + pertvenus[i] + pertmars[i];
	}
	dmat_times_vec(Q_transp, thirdbody_acc_cart, thirdbody_acc_rotrad, 3, 3); //convert to rotating radial frame
	
	//perturbations due to oblateness
	if (center == EARTH) {
		AccelEarthObl(r_vec_ECI, pertharm_cart);							// J2 and J3
		AccelEarthObl2(r_vec_ECI, pertharm_cart2);							// just J2
		dmat_times_vec(Q_transp, pertharm_cart, pertharm_rotrad, 3, 3);		//convert to rotating radial frame
		dmat_times_vec(Q_transp, pertharm_cart2, pertharm_rotrad2, 3, 3);	//convert to rotating radial frame

		// Calculate higher order Earth oblateness
		AccelHarmonic(r_vec_ECI, inert2geo, GM_JGM3, R_JGM3, CS_JGM3, 20, 20, pertharm_cart_CS);
		dmat_times_vec(Q_transp, pertharm_cart_CS, pertharm_rotrad_CS, 3, 3); //convert to rotating radial frame
	}
	else {
		AccelHarmonic (r_vec_LCI, inert2selen, GM_LP100J, R_LP100J, CS_LP100J, 20, 20, pertharm_cart);
		dmat_times_vec(Q_transp, pertharm_cart, pertharm_rotrad, 3, 3); //convert to rotating radial frame
	}
	
	
	//perturbations due to solar pressure
	solar_rad_pressure_cart(julian_date, x[7], r_vec, center, pertsolp_cart);
	dmat_times_vec(Q_transp, pertsolp_cart, pertsolp_rotrad, 3, 3);	//convert to rotating radial frame

	//--------------------------------------------------------------------------------------------

	// Calculate delta V
	deltaV += thrust/x[7] * rpar[2] * (*t-old_t);

	// Calculate power integrals
	for(i=1;i<=3;i++) {
		craft2sun[i] = sun_cart_coords_ECI[i] - r_vec_ECI[i];
	}
	if(center == EARTH)
		solar_power = SOLARINTENSITY*aref*area_efficiency*conversion_efficiency*shadowfunc(sun_cart_coords_ECI,r_vec_ECI,R_JGM3)*sin(sunangle(craft2sun,u));
	else
		solar_power = SOLARINTENSITY*aref*area_efficiency*conversion_efficiency*shadowfunc(sun_cart_coords_ECI,r_vec_LCI,R_LP100J)*sin(sunangle(craft2sun,u));
	thruster_energy += thruster_power * rpar[2] * (*t-old_t);
	solar_energy += solar_power * rpar[2] * (*t-old_t);
	battery_level = solar_energy - thruster_energy + 2.8*3.6e6;
	net_power = solar_power - thruster_power;
	net_energy += net_power * rpar[2] * (*t-old_t);
	if(net_energy > 2.8*3.6e6) {
		net_energy = 2.8*3.6e6;
	} else if (net_energy < 0.0) {
		net_energy = 0.0;
	}
	old_t = *t;



	/************ return values ************/

	outfun[0] =		cart_inert_coords_ECI[1];
	outfun[1] =		cart_inert_coords_ECI[2];
	outfun[2] =		cart_inert_coords_ECI[3];
	outfun[3] =		cart_inert_coords_ECI[4];
	outfun[4] =		cart_inert_coords_ECI[5];
	outfun[5] =		cart_inert_coords_ECI[6];

	outfun[6] =		r_ECI;
	outfun[7] =		r_LCI;
	outfun[8] =		v_ECI;
	outfun[9] =		v_LCI;
	outfun[10] =	v_ECI*v_ECI/2.0-gmue[EARTH]/r_ECI;				// epsilon (orbital energy wrt the earth)
	outfun[11] =	v_LCI*v_LCI/2.0-gmue[MOON]/r_LCI;				// epsilon (orbital energy wrt the moon)
	outfun[12] =	gamma;

	outfun[13] =	kep_coords_ECI[1];
	outfun[14] =	kep_coords_ECI[2];
	outfun[15] =	kep_coords_ECI[3];
	outfun[16] =	kep_coords_ECI[4];
	outfun[17] =	kep_coords_ECI[5];
	outfun[18] =	kep_coords_ECI[6];

	outfun[19] =	cart_inert_coords_SEL[1];
	outfun[20] =	cart_inert_coords_SEL[2];
	outfun[21] =	cart_inert_coords_SEL[3];
	outfun[22] =	cart_inert_coords_SEL[4];
	outfun[23] =	cart_inert_coords_SEL[5];
	outfun[24] =	cart_inert_coords_SEL[6];

	outfun[25] =	cart_inert_coords_LCI[1];
	outfun[26] =	cart_inert_coords_LCI[2];
	outfun[27] =	cart_inert_coords_LCI[3];
	outfun[28] =	cart_inert_coords_LCI[4];
	outfun[29] =	cart_inert_coords_LCI[5];
	outfun[30] =	cart_inert_coords_LCI[6];
		
	outfun[31] =	kep_coords_SEL[1];
	outfun[32] =	kep_coords_SEL[2];
	outfun[33] =	kep_coords_SEL[3];
	outfun[34] =	kep_coords_SEL[4];
	outfun[35] =	kep_coords_SEL[5];
	outfun[36] =	kep_coords_SEL[6];

	outfun[37] =	kep_coords_LCI[3];		// a, e and ny in LCI and SEL frame identical!!
	outfun[38] =	kep_coords_LCI[4];
	outfun[39] =	kep_coords_LCI[5];

	outfun[40] =	moon_cart_coords_ECI[1];
	outfun[41] =	moon_cart_coords_ECI[2];
	outfun[42] =	moon_cart_coords_ECI[3];
	outfun[43] =	dvec_abs(moon_cart_coords_ECI,3);

	outfun[44] =	sun_cart_coords_ECI[1];
	outfun[45] =	sun_cart_coords_ECI[2];
	outfun[46] =	sun_cart_coords_ECI[3];
	outfun[47] =	dvec_abs(sun_cart_coords_ECI,3);

	outfun[48] =	earth_coords_SEL[1];
	outfun[49] =	earth_coords_SEL[2];
	outfun[50] =	earth_coords_SEL[3];
	
//	outfun[51] =	L1_sel[1];
//	outfun[52] =	L1_sel[2];
//	outfun[53] =	L1_sel[3];

//	outfun[54] =	dist_sat_L1;
//	outfun[55] =	jacobi;

//	outfun[56] =	kep_coords_MOP[3];
//	outfun[57] =	kep_coords_MOP[4];
//	outfun[58] =	kep_coords_MOP[5];

	outfun[59] =	phase_moon;

	outfun[60] =	thrust;

	outfun[61] =	thrust/x[7];
	outfun[62] =	dvec_abs(pertearth,3);
	outfun[63] =	dvec_abs(pertmoon,3);
	outfun[64] =	dvec_abs(pertsun,3);
	outfun[65] =	dvec_abs(pertjupiter,3);
	outfun[66] =	dvec_abs(pertvenus,3);
	outfun[67] =	dvec_abs(pertmars,3);
	outfun[68] =	dvec_abs(pertharm_rotrad,3);
	outfun[69] =	dvec_abs(pertsolp_rotrad,3);

	outfun[70] =	moon_coords_kep_ECI[1];
	outfun[71] =	moon_coords_kep_ECI[2];
	outfun[72] =	moon_coords_kep_ECI[3];
	outfun[73] =	moon_coords_kep_ECI[4];
	outfun[74] =	moon_coords_kep_ECI[5];
	outfun[75] =	moon_coords_kep_ECI[6];

	outfun[76] =	deltaV;
	outfun[77] =	kep_coords_ECI[1]*(1-kep_coords_ECI[2]);

//	outfun[78] =	outfun[77]/2.2668e7 - 1.0;
//	outfun[79] =	1.0 - kep_coords_ECI[1]/2.4455e7;
//	outfun[80] =	kep_coords_ECI[1]/2.0e7 - 1.0;
//	outfun[81] =	kep_coords_ECI[2] - 0.732;
//	outfun[82] =	kep_coords_ECI[3]*r2d/21.7 - 1.0;
//	outfun[83] =	x[7]/250.0 - 1.0;

	outfun[84] =	thruster_energy;
	outfun[85] =	solar_energy;
	outfun[86] =	dvec_abs(pertharm_rotrad_CS,3);
	outfun[87] =	net_energy;
	outfun[88] =	battery_level;
	outfun[89] =	dvec_abs(pertharm_rotrad2,3);

	// Free memory
	free_dvector(mee_coords_ECI,1,6);
	free_dvector(mee_coords_LCI,1,6);
	free_dvector(kep_coords_ECI,1,6);	
	free_dvector(kep_coords_LCI,1,6);
	free_dvector(kep_coords_SEL,1,6);	
	free_dvector(cart_inert_coords_ECI,1,6);
	free_dvector(cart_inert_coords_LCI,1,6);
	free_dvector(cart_inert_coords_SEL,1,6);
	free_dvector(moon_cart_coords_ECI,1,6);
	free_dvector(moon_coords_kep_ECI,1,6);
	free_dvector(sun_cart_coords_ECI,1,6);
	free_dvector(earth_coords_SEL,1,3);
	free_dvector(earth_coords_LCI,1,3);
	free_dvector(r_vec_ECI,1,3);
	free_dvector(v_vec_ECI,1,3);
	free_dvector(r_vec_LCI,1,3);
	free_dvector(v_vec_LCI,1,3);
	free_dvector(v_vec_norm,1,3);
	free_dvector(v_vec_norm_rotrad,1,3);
	free_dmatrix(Q_transp,1,3,1,3);

	free_dmatrix(inert2selen,1,3,1,3);
	free_dmatrix(inert2geo,1,3,1,3);

	free_dvector(craft2sun,1,3);

	free_dvector(pertearth,1,3);
	free_dvector(pertsun,1,3);
	free_dvector(pertmoon,1,3);
	free_dvector(pertjupiter,1,3);
	free_dvector(pertvenus,1,3);
	free_dvector(pertmars,1,3);

	free_dvector(thirdbody_acc_cart,1,3);
	free_dvector(thirdbody_acc_rotrad,1,3);

	free_dvector(pertharm_cart,1,3);
	free_dvector(pertharm_rotrad,1,3);
	
	free_dvector(pertharm_cart2,1,3);
	free_dvector(pertharm_rotrad2,1,3);
	
	free_dvector(pertharm_cart_CS,1,3);
	free_dvector(pertharm_rotrad_CS,1,3);
	
	free_dvector(pertsolp_cart,1,3);
	free_dvector(pertsolp_rotrad,1,3);
}