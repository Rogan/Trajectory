/*	This file contains no functions, but rather sets some universal constants 
	declared in the header file
*/

#include "constants.h"

/*****************************************************************************/

/* Migrated from kepler.c */

const double AU = 149597870.66;

// These two arrays could probably be more neatly presented as enums or merged as structs

/* JPL DE-405 GRAVITATIONAL CONSTANTS in m�/s� */
const double gmue[12] = {
		    0.0		*1.0e9,
	    22032.080	*1.0e9,		//Mercury
	   324858.599	*1.0e9,		//Venus
	   398600.433	*1.0e9,		//Earth	
	    42828.314	*1.0e9,		//Mars
	126712767.863	*1.0e9,		//Jupiter
	 37940626.063	*1.0e9,		//Saturn
	  5794549.007	*1.0e9,		//Uranus
	  6836534.064	*1.0e9,		//Neptune
	      981.601	*1.0e9,		//Pluto
	     4902.80058	*1.0e9,		//Moon
 132712440017.987	*1.0e9		//Sun
};

/* NAIF object ID codes */
const int naif_id[12] = {
	0,
	1,		//Mercury barycenter = Mercury center
	2,		//Venus barycenter = Venus center
	399,	//Earth	center
	4,		//Mars barycenter ~= Mars center
	5,		//Jupiter barycenter
	6,		//Saturn barycenter
	7,		//Uranus barycenter
	8,		//Neptune barycenter
	9,		//Pluto barycenter
	301,	//Moon center
	1		//Sun
};
	

/*****************************************************************************/

	
// Moon gravity field LP100J (Konopliv)

  // Gravitational coefficients C, S are efficiently stored in a single
  // array CS. The lower triangle matrix CS holds the non-sectorial C
  // coefficients C_n,m (n>m). Sectorial C coefficients C_n,n are the 
  // diagonal elements of CS (n=m) and the upper triangular matrix stores
  // the S_n,m (for m!=0) coefficients in columns, for the same degree n.
  // Mapping of CS to C, S is achieved through 
  // C_n,m = CS(n,m), S_n,m = CS(m-1,n) ((Table runs from index(0,0)))

const double R_LP100J = .1738000000000000E+04;			// Radius Moon [km]; LP100J
const double GM_LP100J =  .4902800476000000E+04;		// Reference gravity coefficient [km�/s�]; LP100J

const double CS_LP100J[21][21] = {
{ 0.000000E+00,	0.000000E+00, 6.183557E-09,	5.903576E-06, 1.549241E-06,
 -3.499878E-06,-2.030386E-06,-8.547835E-08,	7.852770E-07, 1.375669E-08,
 -6.069273E-07,	3.510464E-07, 8.785556E-07,-7.878745E-09,-3.083653E-07,
  1.889323E-07,	4.668953E-07,-2.892686E-08,-5.039355E-08, 6.158279E-08,	
 -2.585099E-07														 },
{ 0.000000E+00,	0.000000E+00, 7.684213E-09,	1.669827E-06,-1.511609E-06,
  1.814513E-07,-2.758493E-07, 2.394904E-07,	1.535459E-07,-1.036544E-07,
 -1.181695E-08,-1.022052E-07, 7.338319E-08,	2.368704E-08, 9.911799E-09,
  2.379832E-08,-4.203777E-09,-3.204361E-08,-1.481146E-08,-8.404291E-09,
  3.061872E-09														 },
{-2.032162E-04,	1.276073E-08, 2.235138E-05,-2.460004E-07,-8.025185E-07,
  2.862258E-07,-7.282225E-08, 3.353339E-08,	9.408261E-09, 1.745566E-08,
  4.760463E-09,	3.369715E-09,-1.485740E-08,-5.023367E-09,-8.002114E-11,
 -1.444098E-09,	7.563582E-10,-4.742682E-10,	9.311577E-10, 1.185745E-09,
 -1.057544E-10														 },
{-8.407298E-06,	2.846102E-05, 4.850405E-06,	1.713208E-06, 8.331840E-08,
  3.068594E-10,-1.526551E-08, 1.814770E-09,-6.585045E-10,-1.191681E-09,
  8.566674E-10,	8.688822E-10,-1.437570E-10,-1.556399E-10,-4.426564E-10,
 -1.745485E-10,	5.545983E-11, 1.979120E-10,	7.375470E-11,-7.054657E-11,
 -3.474563E-11														 },
{ 9.645077E-06,-5.703048E-06,-1.589100E-06,-8.189809E-08,-1.275449E-07,
 -6.797378E-09,-8.311863E-09, 3.801287E-10,	5.079819E-10,-3.696206E-10,
 -1.863135E-11,	1.055713E-10, 2.603273E-11,	2.251636E-11,-3.015668E-12,
 -8.651339E-12,-4.739445E-12,-4.526767E-12,	5.122471E-13, 4.684512E-12,
 -8.725647E-13},
{-7.341367E-07,-8.830035E-07, 7.076999E-07,	1.514322E-08, 2.170378E-08,
  7.698559E-09, 1.686082E-09, 7.241261E-11,-6.042409E-11,-3.773668E-11,
 -1.363160E-11,-5.916106E-12, 1.382827E-12,	1.993054E-12, 1.620916E-12,
 -1.081831E-13,-4.060517E-13,-4.159690E-13,	2.108235E-13, 3.058644E-14,
  6.942523E-14},
{ 1.367853E-05,	1.194777E-06,-5.438400E-07,-6.751643E-08, 1.406861E-09,
  1.135188E-09,-1.098313E-09,-2.982575E-11,	1.695067E-11,-4.926613E-14,
 -9.495022E-13,-1.223667E-12,-1.487791E-13,-1.075632E-13, 8.889327E-14,
  7.804940E-14,-2.737567E-14,-2.618900E-14,-1.529492E-14, 4.422880E-15,
 -3.110448E-15},
{ 2.167857E-05,	5.473658E-06,-6.575417E-08,	8.654178E-09,-1.963683E-09,
 -9.570447E-11,-7.124501E-11,-3.329078E-11,	2.696036E-12,-7.394317E-13,
  3.181235E-13,-4.502723E-15,-4.937402E-14,-1.889780E-14, 4.812295E-15,
  6.746356E-15,	4.300420E-16,-6.388434E-17,-2.753577E-17, 4.217441E-16,
  1.408565E-17},
{ 9.523089E-06,-3.552751E-08, 2.512294E-07,-1.862059E-08, 4.373439E-09,
 -2.116046E-10,-4.433097E-11,-7.873775E-12,-3.193752E-12, 1.919989E-13,
  6.633726E-17,	1.504723E-15, 1.882634E-15,	8.180925E-16,-5.425901E-16,
 -4.136312E-16,-1.450075E-16, 1.213584E-16,-5.575990E-17,-1.871685E-18,
  5.928348E-18},
{-1.550488E-05,	1.280791E-06, 1.401844E-07,-1.545108E-08,-1.688192E-09,
 -1.539007E-10,-2.936081E-11,-7.509049E-12,-4.182886E-13,-7.332206E-14,
 -7.230301E-15,-2.013232E-16,-3.403849E-16,	8.005804E-17, 4.953694E-18,
 -4.488494E-18,	7.918075E-18, 5.222037E-18,-1.938616E-18, 5.598880E-19,
 -3.695024E-19},
{-4.636715E-06,	5.194869E-07, 1.880154E-08,	2.083116E-09,-2.109229E-09,
  4.883792E-11,-6.212629E-13,-3.250805E-12,-4.186040E-13,-8.867469E-14,
  4.171993E-15,-3.468548E-16,-1.464066E-17,	2.600524E-18, 8.716276E-18,	
 -1.624860E-18,	1.120071E-18,-3.949450E-19,-3.522040E-19,-1.320613E-20,
  3.370133E-20},
{-4.604805E-06,	7.269174E-09, 3.901520E-08,	1.500363E-09,-4.074805E-10, 
  7.411216E-12,	1.939377E-12,-1.212847E-13,-1.082595E-13,-1.269484E-14,
 -4.465511E-15,-6.066655E-16, 1.089928E-17,-2.712667E-18,-1.623186E-19,
 -1.178411E-19,	2.208415E-20,-4.531702E-21,	1.326693E-20, 1.212673E-21,
 -3.583738E-21},
{-9.453210E-06,-3.329599E-07,-1.068457E-08,	2.716012E-09, 3.073795E-10,
 -1.985543E-12,	1.969692E-12, 4.674400E-13,	1.793645E-14,-2.894070E-15,
 -1.053973E-15,-3.792094E-17, 4.087941E-18,-9.761491E-19,-8.516684E-20,
  5.309678E-21,	2.086210E-21, 1.385106E-22,-8.019025E-23, 8.041888E-22,
  3.974132E-23},
{ 1.363850E-06,	6.530422E-07,-1.450923E-07,-9.929775E-10, 2.298457E-10,
 -2.467626E-11,	4.910397E-14, 1.897512E-14,-3.006884E-15,-7.922713E-16,
 -7.400408E-17,-9.293750E-18,-1.543956E-18,	8.482520E-19, 1.568215E-20,
  2.078962E-21,-1.098570E-22,-1.000096E-22,-6.992009E-23, 2.654856E-23,
  4.061647E-24},
{ 2.471068E-06,	3.491508E-07, 1.299892E-08,	1.838746E-09,-5.937240E-11,
 -1.305975E-11,-5.857007E-13,-5.104621E-14,	3.028876E-16, 1.139638E-16,
 -7.417176E-18,-1.063501E-17,-1.128790E-18,	3.854369E-21,-5.886268E-21,
  2.970664E-22,-1.147007E-22, 6.634346E-24,-1.413362E-24,-1.997281E-24,	
 -3.410010E-25},
{-2.703449E-07,-4.850783E-07,-6.994233E-09,-2.508492E-09,-1.528978E-10,
 -2.452747E-12,	2.616560E-14, 6.161741E-14,	5.925688E-15, 7.858544E-18,
 -1.098021E-17,-3.774952E-18,-2.034814E-19,	4.199203E-21, 1.061324E-21,
  1.780790E-22,-1.051947E-23,-2.310145E-26,-9.772859E-25,-2.872782E-26,
  4.751995E-26},
{ 2.259089E-06,-5.794624E-08, 5.139173E-08,	7.803275E-11, 4.518107E-11,
  5.525298E-12,	5.141238E-13,-6.880886E-15,-6.303059E-16,-2.122778E-16,
  3.398018E-18,	3.658667E-19,-6.948115E-20,-7.372335E-21,-8.349957E-22,
  1.050372E-25,-1.156946E-23, 5.794636E-25,-5.220583E-26,-1.168617E-26,
 -3.780774E-27},
{-6.183981E-06,	3.097805E-07,-3.100659E-09,-2.342974E-10, 8.887374E-11,
  2.078475E-12,	2.533667E-13,-3.315183E-14,-4.546275E-16, 5.806060E-17,
  7.244710E-18,	3.456877E-19, 3.394437E-22,-2.334577E-22, 3.088407E-22,	
  7.581203E-24,-1.548707E-24,-2.740148E-25,	4.919814E-28,-9.428684E-29,	
  -6.610233E-28},
{-3.476290E-06,	7.228552E-08,-4.601468E-09,	1.166260E-09,-7.104547E-11,	
 -8.607658E-13,-4.005427E-13,-2.005942E-15,	4.998133E-16, 3.115842E-17,
  1.318464E-19,-8.820551E-20, 1.858943E-20,-3.244006E-22,-1.114326E-23,
 -3.849166E-24,	3.863587E-25, 1.110594E-25,	6.816005E-29, 3.060220E-29,
  3.224313E-29},
{-3.120009E-07,-1.410628E-07, 1.047452E-08,-1.223891E-09,-5.809645E-11,
  2.825311E-13,	2.072534E-15, 9.724138E-15,	1.543957E-16,-1.198094E-18,	
 -1.183074E-18,	2.628569E-20, 3.703394E-21,	6.785224E-23,-2.545975E-23,	
  1.026046E-24,	9.165692E-26,-1.317408E-26,	2.689282E-27, 1.348206E-28,
 -1.150088E-30},
{ 3.299472E-06,	3.102846E-08, 1.101222E-08,	3.297104E-10, 4.518263E-11,
  8.093938E-13,-7.185532E-14,-5.207906E-15,	1.796679E-16, 5.653554E-18,
  2.109493E-19,-5.128015E-21,-8.206011E-23,	1.064948E-22, 6.373959E-24,
 -1.222753E-26,	4.715335E-26, 3.289519E-27,	1.488307E-28,-1.462056E-29,	
  1.383488E-30}
};

  // Earth gravity field JGM3 (Joint Gravity Model 3, Tapley et. al., 1996)

const double R_JGM3 = 6378.1363;			// Radius Earth [km]; JGM3
const double GM_JGM3 = 398600.4415;			// Reference gravity coefficient [km�/s�]; JGM3

const double CS_JGM3[21][21] = {
{ 1.000000e+00, 0.000000e+00, 1.543100e-09, 2.680119e-07,-4.494599e-07,     
 -8.066346e-08, 2.116466e-08, 6.936989e-08, 4.019978e-08, 1.423657e-08,     
 -8.128915e-08,-1.646546e-08,-2.378448e-08, 2.172109e-08, 1.443750e-08,     
  4.154186e-09, 1.660440e-08,-1.427822e-08,-1.817656e-08, 7.160542e-11,     
  2.759192e-09                                                       },     
{ 0.000000e+00, 0.000000e+00,-9.038681e-07,-2.114024e-07, 1.481555e-07,     
 -5.232672e-08,-4.650395e-08, 9.282314e-09, 5.381316e-09,-2.228679e-09,     
 -3.057129e-09,-5.097360e-09, 1.416422e-09,-2.545587e-09,-1.089217e-10,     
 -1.045474e-09, 7.856272e-10, 2.522818e-10, 3.427413e-10,-1.008909e-10,     
  3.216826e-10                                                       },     
{-1.082627e-03,-2.414000e-10, 1.574536e-06, 1.972013e-07,-1.201129e-08,     
 -7.100877e-09, 1.843134e-10,-3.061150e-09,-8.723520e-10,-5.633921e-10,     
 -8.989333e-10,-6.863521e-10, 9.154575e-11, 3.005522e-10, 5.182512e-11,     
  3.265044e-11,-4.271981e-11, 1.297841e-11,-4.278803e-12,-1.190759e-12,     
  3.778260e-11                                                       },     
{ 2.532435e-06, 2.192799e-06, 3.090160e-07, 1.005589e-07, 6.525606e-09,     
  3.873005e-10,-1.784491e-09,-2.636182e-10, 9.117736e-11, 1.717309e-11,     
 -4.622483e-11,-2.677798e-11, 9.170517e-13,-2.960682e-12,-3.750977e-12,     
  1.116419e-12, 5.250141e-12, 2.159727e-12, 1.105860e-13,-3.556436e-13,     
 -1.178441e-12                                                       },     
{ 1.619331e-06,-5.087253e-07, 7.841223e-08, 5.921574e-08,-3.982396e-09,     
 -1.648204e-09,-4.329182e-10, 6.397253e-12, 1.612521e-11,-5.550919e-12,     
 -3.122269e-12, 1.982505e-12, 2.033249e-13, 1.214266e-12,-2.217440e-13,     
  8.637823e-14,-1.205563e-14, 2.923804e-14, 1.040715e-13, 9.006136e-14,     
 -1.823414e-14                                                       },     
{ 2.277161e-07,-5.371651e-08, 1.055905e-07,-1.492615e-08,-2.297912e-09,     
  4.304768e-10,-5.527712e-11, 1.053488e-11, 8.627743e-12, 2.940313e-12,     
 -5.515591e-13, 1.346234e-13, 9.335408e-14,-9.061871e-15, 2.365713e-15,     
 -2.505252e-14,-1.590014e-14,-9.295650e-15,-3.743268e-15, 3.176649e-15,     
 -5.637288e-17                                                       },     
{-5.396485e-07,-5.987798e-08, 6.012099e-09, 1.182266e-09,-3.264139e-10,     
 -2.155771e-10, 2.213693e-12, 4.475983e-13, 3.814766e-13,-1.846792e-13,     
 -2.650681e-15,-3.728037e-14, 7.899913e-15,-9.747983e-16,-3.193839e-16,     
  2.856094e-16,-2.590259e-16,-1.190467e-16, 8.666599e-17,-8.340023e-17,     
 -8.899420e-19                                                       },     
{ 3.513684e-07, 2.051487e-07, 3.284490e-08, 3.528541e-09,-5.851195e-10,     
  5.818486e-13,-2.490718e-11, 2.559078e-14, 1.535338e-13,-9.856184e-16,     
 -1.052843e-14, 1.170448e-15, 3.701523e-16,-1.095673e-16,-9.074974e-17,     
  7.742869e-17, 1.086771e-17, 4.812890e-18, 2.015619e-18,-5.594661e-18,     
  1.459810e-18                                                       },     
{ 2.025187e-07, 1.603459e-08, 6.576542e-09,-1.946358e-10,-3.189358e-10,     
 -4.615173e-12,-1.839364e-12, 3.429762e-13,-1.580332e-13, 7.441039e-15,     
 -7.011948e-16, 2.585245e-16, 6.136644e-17, 4.870630e-17, 1.489060e-17,     
  1.015964e-17,-5.700075e-18,-2.391386e-18, 1.794927e-18, 1.965726e-19,     
 -1.128428e-19                                                       },     
{ 1.193687e-07, 9.241927e-08, 1.566874e-09,-1.217275e-09,-7.018561e-12,     
 -1.669737e-12, 8.296725e-13,-2.251973e-13, 6.144394e-14,-3.676763e-15,     
 -9.892610e-17,-1.736649e-17, 9.242424e-18,-4.153238e-18,-6.937464e-20,     
  3.275583e-19, 1.309613e-19, 1.026767e-19,-1.437566e-20,-1.268576e-20,     
 -6.100911e-21                                                       },     
{ 2.480569e-07, 5.175579e-08,-5.562846e-09,-4.195999e-11,-4.967025e-11,     
 -3.074283e-12,-2.597232e-13, 6.909154e-15, 4.635314e-15, 2.330148e-15,     
  4.170802e-16,-1.407856e-17,-2.790078e-19,-6.376262e-20,-1.849098e-19,     
  3.595115e-20,-2.537013e-21, 4.480853e-21, 4.348241e-22, 1.197796e-21,     
 -1.138734e-21                                                       },     
{-2.405652e-07, 9.508428e-09, 9.542030e-10,-1.409608e-10,-1.685257e-11,     
  1.489441e-12,-5.754671e-15, 1.954262e-15,-2.924949e-16,-1.934320e-16,     
 -4.946396e-17, 9.351706e-18,-9.838299e-20, 1.643922e-19,-1.658377e-20,     
  2.905537e-21, 4.983891e-22, 6.393876e-22,-2.294907e-22, 6.437043e-23,     
  6.435154e-23                                                       },     
{ 1.819117e-07,-3.068001e-08, 6.380398e-10, 1.451918e-10,-2.123815e-11,     
  8.279902e-13, 7.883091e-15,-4.131557e-15,-5.708254e-16, 1.012728e-16,     
 -1.840173e-18, 4.978700e-19,-2.108949e-20, 2.503221e-20, 3.298844e-21,     
 -8.660491e-23, 6.651727e-24, 5.110031e-23,-3.635064e-23,-1.311958e-23,     
  1.534228e-24                                                       },     
{ 2.075677e-07,-2.885131e-08, 2.275183e-09,-6.676768e-11,-3.452537e-13,     
  1.074251e-12,-5.281862e-14, 3.421269e-16,-1.113494e-16, 2.658019e-17,     
  4.577888e-18,-5.902637e-19,-5.860603e-20,-2.239852e-20,-6.914977e-23,     
 -6.472496e-23,-2.741331e-23, 2.570941e-24,-1.074458e-24,-4.305386e-25,     
 -2.046569e-25                                                       },     
{-1.174174e-07,-9.997710e-09,-1.347496e-09, 9.391106e-11, 3.104170e-13,     
  3.932888e-13,-1.902110e-14, 2.787457e-15,-2.125248e-16, 1.679922e-17,     
  1.839624e-18, 7.273780e-20, 4.561174e-21, 2.347631e-21,-7.142240e-22,     
 -2.274403e-24,-2.929523e-24, 1.242605e-25,-1.447976e-25,-3.551992e-26,     
 -7.473051e-28                                                       },     
{ 1.762727e-08, 6.108862e-09,-7.164511e-10, 1.128627e-10,-6.013879e-12,     
  1.293499e-13, 2.220625e-14, 2.825477e-15,-1.112172e-16, 3.494173e-18,     
  2.258283e-19,-1.828153e-21,-6.049406e-21,-5.705023e-22, 1.404654e-23,     
 -9.295855e-24, 5.687404e-26, 1.057368e-26, 4.931703e-27,-1.480665e-27,     
  2.400400e-29                                                       },     
{-3.119431e-08, 1.356279e-08,-6.713707e-10,-6.451812e-11, 4.698674e-12,     
 -9.690791e-14, 6.610666e-15,-2.378057e-16,-4.460480e-17,-3.335458e-18,     
 -1.316568e-19, 1.643081e-20, 1.419788e-21, 9.260416e-23,-1.349210e-23,     
 -1.295522e-24,-5.943715e-25,-9.608698e-27, 3.816913e-28,-3.102988e-28,     
 -8.192994e-29                                                       },     
{ 1.071306e-07,-1.262144e-08,-4.767231e-10, 1.175560e-11, 6.946241e-13,     
 -9.316733e-14,-4.427290e-15, 4.858365e-16, 4.814810e-17, 2.752709e-19,     
 -2.449926e-20,-6.393665e-21, 8.842755e-22, 4.178428e-23,-3.177778e-24,     
  1.229862e-25,-8.535124e-26,-1.658684e-26,-1.524672e-28,-2.246909e-29,     
 -5.508346e-31                                                       },     
{ 4.421672e-08, 1.958333e-09, 3.236166e-10,-5.174199e-12, 4.022242e-12,     
  3.088082e-14, 3.197551e-15, 9.009281e-17, 2.534982e-17,-9.526323e-19,     
  1.741250e-20,-1.569624e-21,-4.195542e-22,-6.629972e-24,-6.574751e-25,     
 -2.898577e-25, 7.555273e-27, 3.046776e-28, 3.696154e-29, 1.845778e-30,     
  6.948820e-31                                                       },     
{-2.197334e-08,-3.156695e-09, 7.325272e-10,-1.192913e-11, 9.941288e-13,     
  3.991921e-14,-4.220405e-16, 7.091584e-17, 1.660451e-17, 9.233532e-20,     
 -5.971908e-20, 1.750987e-21,-2.066463e-23,-3.440194e-24,-1.487095e-25,     
 -4.491878e-26,-4.558801e-27, 5.960375e-28, 8.263952e-29,-9.155723e-31,     
 -1.237749e-31                                                       },     
{ 1.203146e-07, 3.688524e-09, 4.328972e-10,-6.303973e-12, 2.869669e-13,     
 -3.011115e-14, 1.539793e-15,-1.390222e-16, 1.766707e-18, 3.471731e-19,     
 -3.447438e-20, 8.760347e-22,-2.271884e-23, 5.960951e-24, 1.682025e-25,     
 -2.520877e-26,-8.774566e-28, 2.651434e-29, 8.352807e-30,-1.878413e-31,     
  4.054696e-32                                                       }     
};