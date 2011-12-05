/* -------------------------------------------------------------------------------
 * Gesop_Model.Description
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
#include "constants.h"
#include "problem.h"

/*--------------
 *   Constants -
 *--------------*/

	const int Number_of_Phases = 1;

	// Max 60 chars
	static char Problem_Description[] = "Low thrust transfer for Lunar Mission BW1";

	static DLL_Item Phase_Description[1]=
	{
//		{"Ascent",	-1, "Ascend beyond VABs using arcjet",			-1},
//		{"Cruise",	-1, "Cruise to lunar SOI using PPTs",			-1},
//		{"Capture",	-1, "Capture into lunar orbit using arcjet",	-1},
//		{"Descent",	-1, "Descend to LLO using PPTs",				-1},
		{"Science",	-1, "Science orbit: stationkeeping only",		-1},
	};

	static DLL_Unit_Item Time_Description  = 
		{"Lnmee",	-1, "Normalized Modified Equinoctial Element L",	-1, "None",	-1};
	
	static DLL_Item Initial_Cost_Item =		{"",	-1, "",	-1};
	
	static DLL_Item Integral_Cost_Item =	{"",	-1, "",	-1};
	
	static DLL_Item Terminal_Cost_Item =	{"",	-1, "",	-1};
	
	static DLL_Connected_Unit_Item States[9] =
	{
		{"pmee",		-1, "Modified Equinoctial Element p",			-1, "Meter",	-1, SOFT_CONNECT},
		{"fmee",		-1, "Modified Equinoctial Element f",			-1, "None",		-1, SOFT_CONNECT},
		{"gmee",		-1, "Modified Equinoctial Element g",			-1, "None",		-1, SOFT_CONNECT},
		{"hmee",		-1, "Modified Equinoctial Element h",			-1, "None",		-1, SOFT_CONNECT},
		{"kmee",		-1, "Modified Equinoctial Element k",			-1, "None",		-1, SOFT_CONNECT},
		{"Lmee",		-1, "Modified Equinoctial Element L",			-1, "Radian",	-1, SOFT_CONNECT},
		{"RealTime",	-1, "Real Time",								-1, "Day",		-1, SOFT_CONNECT},
		{"Mass",		-1, "Mass",										-1, "Kilogram",	-1, SOFT_CONNECT},
		{"Battery",		-1, "Battery charge",							-1,	"Joule",	-1, SOFT_CONNECT},
	};

	static DLL_Control_Item Control_Vector[4] =
	{
		{"ur",			-1, "Thrust component in r direction",				-1, "None",		-1, PIECEWISE_CONSTANT},
		{"utheta",		-1, "Thrust component perp. to r in orbit plane",	-1, "None",		-1, PIECEWISE_CONSTANT},
		{"uh",			-1, "Thrust component perp. to orbit plane",		-1, "None",		-1, PIECEWISE_CONSTANT},
		{"thrust_mag",	-1, "Normalised thrust magnitude",					-1, "None",		-1, PIECEWISE_CONSTANT},
	};

	static DLL_Connected_Unit_Item RPars[3] =
	{
		{"Ini_T",			-1, "Start time (JD)",						-1,	"Day",		-1, SOFT_CONNECT},
		{"Thrust_Level",	-1, "Thrust magnitude",						-1, "None",		-1, NO_CONNECT},		
		{"Delta_L",			-1, "Phase length (anomaly)",				-1, "Radian",	-1, NO_CONNECT},
	};

	static DLL_Connected_Item IPars[1] = 
	{
		{"", -1, "", -1, NO_CONNECT}
	};

	// In a DLL_Constraint_Item, the second last parameter represents whether the constraint is an equality constraint
	// In a DLL_Constraint_Item, the last parameter represents whether the constraint is enforced
	static DLL_Constraint_Item IBCon[5] =
	{
		{"IBC_a_LCI",			-1, "Semimajor axis",					-1, true,	true},
		{"IBC_e_LCI",			-1, "Eccentricity",						-1, false,	true},
		{"IBC_i_LCI",			-1, "Inclination",						-1, true,	true},
		{"IBC_t",				-1, "Time",								-1, true,	false},
		{"IBC_m",				-1, "Mass",								-1, true,	false},
	};
	
	static DLL_Constraint_Item FBCon[1]=
	{
		{"",	-1, "",	-1, false,	true},
	};

	static DLL_Constraint_Item PCon[5] =
	{
		{"PC_u",					-1, "Unit control vector",							-1, true,	true},
		{"PC_r_ECI",				-1, "Over 100km Earth altitude",					-1, false,	true},
		{"PC_e_ECI",				-1, "In Earth orbit",								-1, false,	false},
		{"PC_r_LCI",				-1, "Over 50km lunar altitude",						-1, false,	true},
		{"PC_i_ECI",				-1, "Prograde orbit",								-1, false,	true},
	};

	static DLL_Constraint_Item PrCon[1] = 
	{
		{"PRC_deltaL",				-1, "Time step must be positive",					-1, false,	true},
	};

	

/*--------------
 *  Operations -
 *--------------*/
	void __cdecl No_Phases (
		int		*num_phases,	// Number of phases
		int		*error			// Error flag
		)
	/* SYNOPSIS:
	 * returns the number of phases specified. */
	{
#ifdef DEBUG_MODE
printf("No. Phases\n");_flushall();
#endif

	 	*error		= 0;
		*num_phases	= Number_of_Phases;
	}











	void __cdecl General_Description (
		char	**string,		// Description string for problem
		int		*error,			// Error flag
		int		*str_len		// Length of description string
		)
	/* SYNOPSIS:
	 * returns a decription of the trajectory optimization problem
	 * (which may well be an empty string). */
	{
#ifdef DEBUG_MODE
printf("General Description\n");_flushall();
#endif

		*error	= 0;
		*string	= &Problem_Description[0];
		*str_len = strlen(*string);
	}












	void __cdecl Phase_Info (
		const int		*phase,					// Current phase number
		Type_Of_Phase	*phase_type,			// Type of phase (numeric, analytic)
		DLL_Item		*phase_desc,			// Phase description
		DLL_Unit_Item	*phase_time,			// Time description
		DLL_Item		*ini_cost,				// Initial cost description
		DLL_Item		*integ_cost,			// Integral cost description
		DLL_Item		*term_cost,				// Terminal cost description
		int				*error					// Error flag
		)
	/* SYNOPSIS:
	 * returns for the specified phase, the type of the phase, its
	 * name and description, name and description of the
	 * independant variable and for each cost term. */
	{
#ifdef DEBUG_MODE
printf("Phase Info\n");_flushall();
#endif

		*error = 0;

		/*	what kind of phase */
		*phase_type					= Numeric;

		Phase_Description[*phase-1].Name_len	= strlen(Phase_Description[*phase-1].Name);
		Phase_Description[*phase-1].Desc_len	= strlen(Phase_Description[*phase-1].Desc);

		phase_desc->Desc			= Phase_Description[*phase-1].Desc;
		phase_desc->Name			= Phase_Description[*phase-1].Name;

		phase_desc->Desc_len        = Phase_Description[*phase-1].Desc_len;
		phase_desc->Name_len        = Phase_Description[*phase-1].Name_len;

		/*	time info */
		Time_Description.Name_len	= strlen(Time_Description.Name);
		Time_Description.Desc_len	= strlen(Time_Description.Desc);
		Time_Description.Unit_len	= strlen(Time_Description.Unit);

		phase_time->Desc			= Time_Description.Desc;
		phase_time->Name			= Time_Description.Name;
		phase_time->Unit			= Time_Description.Unit;

		phase_time->Desc_len		= Time_Description.Desc_len;
		phase_time->Name_len		= Time_Description.Name_len;
		phase_time->Unit_len		= Time_Description.Unit_len;

		/* ini cost info */
		Initial_Cost_Item.Name_len	= strlen(Initial_Cost_Item.Name);
		Initial_Cost_Item.Desc_len	= strlen(Initial_Cost_Item.Desc);

		ini_cost->Desc					= Initial_Cost_Item.Desc;
		ini_cost->Name					= Initial_Cost_Item.Name;

		ini_cost->Desc_len				= Initial_Cost_Item.Desc_len;
		ini_cost->Name_len				= Initial_Cost_Item.Name_len;

		/*	lagrange cost info
	 	*	if there isn't any, set Integral_Cost := Unspecified.  */
		Integral_Cost_Item.Name_len	= strlen(Integral_Cost_Item.Name);
		Integral_Cost_Item.Desc_len	= strlen(Integral_Cost_Item.Desc);

		integ_cost->Desc			= Integral_Cost_Item.Desc;
		integ_cost->Name			= Integral_Cost_Item.Name;

		integ_cost->Desc_len		= Integral_Cost_Item.Desc_len;
		integ_cost->Name_len		= Integral_Cost_Item.Name_len;
		
		/* terminal cost info */
		Terminal_Cost_Item.Name_len	= strlen(Terminal_Cost_Item.Name);
		Terminal_Cost_Item.Desc_len	= strlen(Terminal_Cost_Item.Desc);
			
		term_cost->Desc				= Terminal_Cost_Item.Desc;
		term_cost->Name				= Terminal_Cost_Item.Name;

		term_cost->Desc_len			= Terminal_Cost_Item.Desc_len;
		term_cost->Name_len			= Terminal_Cost_Item.Name_len;
	}









	void __cdecl State_Info (
		const int				phase,				// Current phase number
		int						*con_item_size,		// Number of states
		DLL_Connected_Unit_Item	**con_item_vec,		// Vector of State descriptions
		int						*error				// Error flag
		)
		/* SYNOPSIS:
		 * returns (implicitly) the dimension and (explicitly) the name and
		 * description information for the states. The Is_Connected flags
		 * tell which states are connected to quantities of the previous
		 * phase; this information is only relevant if Phase >= 2. */
	{
		int i;

#ifdef DEBUG_MODE
printf("State Info\n");_flushall();
#endif

		*error = 0;

		*con_item_size = sizeof(States)/sizeof(DLL_Connected_Unit_Item);

		for (i=0; i<*con_item_size; i++) {
			States[i].Name_len = strlen(States[i].Name);
			States[i].Desc_len = strlen(States[i].Desc);
			States[i].Unit_len = strlen(States[i].Unit);
		}

		*con_item_vec = States;
	}









	void __cdecl Control_Info (
		const int			phase,				// Current phase number
		int					*con_item_size,		// Number of controls
		DLL_Control_Item	**con_item_vec,		// Vector of Control descriptions
		int					*error				// Error flag
		)
		/* SYNOPSIS:
		 * returns (implicitly) the dimension and (explicitly) the name and
		 * description information for the controls. The Is_Connected flags
		 * tell which controls are connected to quantities of the previous
		 * phase; this information is only relevant if Phase >= 2. */
	{
		int i;
		
#ifdef DEBUG_MODE
printf("Control Info\n");_flushall();
#endif

		*error = 0;

		*con_item_size = sizeof(Control_Vector)/sizeof(DLL_Control_Item);
			
		for (i=0; i<*con_item_size; i++) {
			Control_Vector[i].Name_len = strlen(Control_Vector[i].Name);
			Control_Vector[i].Desc_len = strlen(Control_Vector[i].Desc);
			Control_Vector[i].Unit_len = strlen(Control_Vector[i].Unit);
		}
			
		*con_item_vec = Control_Vector;
	}










	void __cdecl Integer_Parameter_Info (
		const int			phase,						// Current phase number
		int					*con_item_size,				// Number of Integer parameters
		DLL_Connected_Item	**con_item_vec,				// Vector of IPar. descriptions
		int					*error						// Error flag
		)
		/* SYNOPSIS:
		 * returns (implicitly) the dimension and (explicitly)
		 * the name and description information for the integer parameters.  */
	{
		int i;

#ifdef DEBUG_MODE
printf("Integer Parameter Info\n");_flushall();
#endif

		*error = 0;

		//*con_item_size = sizeof(IPars)/sizeof(DLL_Connected_Item);
		*con_item_size = 0;

		for (i=0; i<*con_item_size; i++) {
			IPars[i].Name_len = strlen(IPars[i].Name);
			IPars[i].Desc_len = strlen(IPars[i].Desc);
		}		

		*con_item_vec = IPars;
	}











	void __cdecl Real_Parameter_Info (
		const int				phase,				// Current phase number
		int						*con_item_size,		// Number of real parameters
		DLL_Connected_Unit_Item	**con_item_vec,		// Vector of RPar. descriptions
		int						*error				// Error flag
		)
		/* SYNOPSIS:
		 * returns (implicitly) the dimension and (explicitly) the name and
		 * description information for the real parameters. The Is_Connected
		 * flags tell which real parameters are connected to quantities of
		 * the previous phase; this information is only relevant if Phase >= 2. */
	{
		int i;
		
#ifdef DEBUG_MODE
printf("Real Parameter Info\n");_flushall();
#endif

		*error = 0;

		*con_item_size = sizeof(RPars)/sizeof(DLL_Connected_Unit_Item);
		
		for (i=0; i<*con_item_size; i++) {
			RPars[i].Name_len = strlen(RPars[i].Name);
			RPars[i].Desc_len = strlen(RPars[i].Desc);
			RPars[i].Unit_len = strlen(RPars[i].Unit);
		}
		*con_item_vec = RPars;
	}









	void __cdecl Initial_Boundary_Constraint_Info (
		const int			phase,				// Current phase number
		int					*con_item_size,		// Number of initial boundary constraints
		DLL_Constraint_Item	**con_item_vec,		// Vector of IBC descriptions
		int					*error				// Error flag
		)
		/* SYNOPSIS:
		 * returns (implicitly) the dimension and (explicitly) the name,
		 * description, constraint type and enforcement information for the
		 * initial boundary constraints of phase one.  */
	{
		int i;

#ifdef DEBUG_MODE
printf("Initial Boundary Constraint Info\n");_flushall();
#endif

		*error = 0;

		//*con_item_size = 0;
		*con_item_size = sizeof(IBCon)/sizeof(DLL_Constraint_Item);

		for (i=0; i<*con_item_size; i++) {
			IBCon[i].Name_len = strlen(IBCon[i].Name);
			IBCon[i].Desc_len = strlen(IBCon[i].Desc);
		}

		*con_item_vec = IBCon;
	}










	void __cdecl Final_Boundary_Constraint_Info (
		const int			phase,				// Current phase number
		int					*con_item_size,		// Number of final boundary constraints
		DLL_Constraint_Item	**con_item_vec,		// Vector of FBC descriptions
		int					*error				// Error flag
		)
		/* SYNOPSIS:
		 * returns (implicitly) the dimension and (explicitly) the name,
		 * description, constraint type and enforcement  information for the
		 * terminal (final) boundary constraints of the specified phase. */
	{
		int i;

#ifdef DEBUG_MODE
printf("Final Boundary Constraint Info\n");_flushall();
#endif

		*error = 0;

		//*con_item_size = 0;
		*con_item_size = sizeof(FBCon)/sizeof(DLL_Constraint_Item);

		for (i=0; i<*con_item_size; i++) {
			FBCon[i].Name_len = strlen(FBCon[i].Name);
			FBCon[i].Desc_len = strlen(FBCon[i].Desc);
		}
			
		*con_item_vec = FBCon;
	}







	void __cdecl Path_Constraint_Info (
		const int			phase,				// Current phase number
		int					*con_item_size,		// Number of path constraints
		DLL_Constraint_Item	**con_item_vec,		// Vector of PC descriptions
		int					*error				// Error flag
		)
		/* SYNOPSIS:
		 * returns (implicitly) the dimension and (explicitly) the name,
		 * description, constraint type and enforcement information for the
		 * path constraints of the specified phase. */
	{
		int i;
		
#ifdef DEBUG_MODE
printf("Path Constraint Info\n");_flushall();
#endif

		*error = 0;
	
		//*con_item_size = 0;
		*con_item_size = sizeof(PCon)/sizeof(DLL_Constraint_Item);

		for (i=0; i<*con_item_size; i++) {
			PCon[i].Name_len = strlen(PCon[i].Name);
			PCon[i].Desc_len = strlen(PCon[i].Desc);
		}
			
		*con_item_vec = PCon;
	}






	void __cdecl Parameter_Constraint_Info (
		const int			phase,				// Current phase number
		int					*con_item_size,		// Number of parameter constraints
		DLL_Constraint_Item	**con_item_vec,		// Vector of Par.C. descriptions
		int					*error				// Error flag
		)
		/* SYNOPSIS:
		 * returns (implicitly) the dimension and (explicitly) the name,
		 * description, constraint type and enforcement  information for the
		 * parameter constraints of the specified phase.  */
	{
#ifdef DEBUG_MODE
printf("Parameter Constraint Info\n");_flushall();
#endif

		*error = 0;
		
		//*con_item_size = 0;
		*con_item_size = sizeof(PrCon)/sizeof(DLL_Constraint_Item);


		PrCon[0].Name_len = strlen(PrCon[0].Name);
		PrCon[0].Desc_len = strlen(PrCon[0].Desc);

		*con_item_vec = PrCon;
	}