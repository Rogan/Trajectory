/*
 * $Revision: 3264 $ $Date: 2009-06-29 14:17:15 +0200 (Mo, 29 Jun 2009) $
 *
 * Copyright by Astos Solutions GmbH, Stuttgart, Germany
 * All rights reserved.
 * For details on copyright and terms of use see
 * http://www.astos.de/files_copyright_and_terms_of_use.html
 */

/* -------------------------------------------------------------------------------
 *  Gesop-Model.Auxiliary_Scalars
 *
 *  TEMPLATE: No model is defined.
 * ------------------------------------------------------------------------------- */ 
 
#include <stdio.h>
#include "gesop_model.h"


static DLL_Auto_Item AutoItem[1] =
    {
    // Name, Desc, Unit, Pattern, Eval_Type, Index
    {"", -1, "", -1, "None", -1, 0, 0, 0}
    };

static DLL_Value_Item ValueItem[1] =
    {
    // Name, Desc, Unit, Pattern, Value
    {"", -1, "", -1, "None", -1, 0, 0.0}
    };

static DLL_Event_Item EventItem[1] =
    {
    // Name, Desc, Eval_Value, Eval_Type, Index
    {"", -1, "", -1, 0.0, 0, 0}
    };

static DLL_Object_Item ObjectItem[1] =
    {
    // Name, Desc, Group, Kind, Num_AF, Num_Points, Index
    {"", -1, "", -1, "", -1, 0, 0, 0, 0}
    };

/*---------------
 *   Operations -
 *---------------*/
void __cdecl Algebraic_Scalars_Info(
    const int           *phase,         // Current phase number
    int                 *item_size,     // Number of algebraic scalars
    DLL_Pattern_Item    **item_vec,     // Vector of algebraic scalar descriptions
    int                 *error)         // Error flag
{
    *error = 0;
    *item_size = 0;     // @MODEL DEFINITION
    *item_vec = NULL;   // @MODEL DEFINITION
}

void __cdecl Algebraic_Scalars(
    const int               *phase,             // Current phase number
    const int               *dimx,              // Dimension of state vector X
    const int               *dimu,              // Dimension of control vector U
    const int               *dimip,             // Dimension of integer parameter vector
    const int               *dimrp,             // Dimension of real parameter vector
    const Phase_Info_Type   *fazinf,            // Record with additional phase info
    const double            *t,                 // Time t of evaluation
    const double            x[],                // State vector X at time t
    const double            u[],                // Control vector U at time t
    const double            udot[],             // Derivative of the control vector U
    const int               ipar[],             // Integer parameter vector of phase
    const double            rpar[],             // Real parameter vector of phase
    const int               *dimof,             // Dimension of scalar vector
    const int               *f_evaltime,        // Evaluation time flag
    double                  *f_itemvalues,      // Values of algebraic scalars
    int                     *error)             // Error flag
{
    *error = 0;
    // @MODEL DEFINITION
}

void __cdecl Algebraic_Auto_Scalars_Info(
    const int               *phase,             // Current phase number
    int                     *dimof,             // Dimension of scalar vector
    DLL_Auto_Item           **item_vec,         // Values of algebraic scalars
    int                     *error)             // Error flag
{
    int i;
    *error = 0;
    *dimof = 0;
    for (i=0; i<*dimof; i++) {
        AutoItem[i].Name_len = strlen(AutoItem[i].Name);
        AutoItem[i].Desc_len = strlen(AutoItem[i].Desc);
        AutoItem[i].Unit_len = strlen(AutoItem[i].Unit);
    }
    *item_vec = AutoItem;
}

void __cdecl Algebraic_Value_Scalars_Info(
    const int               *phase,             // Current phase number
    int                     *dimof,             // Dimension of scalar vector
    DLL_Value_Item          **item_vec,         // Values of algebraic scalars
    int                     *error)             // Error flag
{
    int i;
    *error = 0;
    *dimof = 0;
    for (i=0; i<*dimof; i++) {
        ValueItem[i].Name_len = strlen(ValueItem[i].Name);
        ValueItem[i].Desc_len = strlen(ValueItem[i].Desc);
        ValueItem[i].Unit_len = strlen(ValueItem[i].Unit);
    }
    *item_vec = ValueItem;
}

void __cdecl Algebraic_Event_Scalars_Info(
    const int               *phase,             // Current phase number
    int                     *dimof,             // Dimension of scalar vector
    DLL_Event_Item          **item_vec,         // Values of algebraic events
    int                     *error)             // Error flag
{
    int i;
    *error = 0;
    *dimof = 0;
    for (i=0; i<*dimof; i++) {
        EventItem[i].Name_len = strlen(EventItem[i].Name);
        EventItem[i].Desc_len = strlen(EventItem[i].Desc);
    }
    *item_vec = EventItem;
}

void __cdecl Algebraic_Items_Info(
    const int           *phase,         // Current phase number
    int                 *item_size,     // Number of algebraic items
    DLL_Object_Item     **item_vec,     // Vector of algebraic items descriptions
    int                 *error)         // Error flag
{
    int i;
    *error = 0;
    *item_size = 0;
    for (i=0; i<*item_size; i++) {
        ObjectItem[i].Name_len = strlen(ObjectItem[i].Name);
        ObjectItem[i].Desc_len = strlen(ObjectItem[i].Desc);
        ObjectItem[i].Group_len = strlen(ObjectItem[i].Group);
    }
    *item_vec = ObjectItem;
}

void __cdecl Algebraic_Items(
    const int               *phase,             // Current phase number
    const int               *dimip,             // Dimension of integer parameter vector
    const int               *dimrp,             // Dimension of real parameter vector
    const Phase_Info_Type   *fazinf,            // Record with additional phase info
    const int               ipar[],             // Integer parameter vector of phase
    const double            rpar[],             // Real parameter vector of phase
    const int               *index,
    const int               *dimaf,
    const int               *dimpt,
    double                  *itemvalues,       // Values of algebraic items (double_matrix)
    int                     *error)             // Error flag
{
    *error = 0;

    // @MODEL DEFINITION
}
