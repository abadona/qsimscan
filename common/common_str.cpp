//////////////////////////////////////////////////////////////////////////////
//// This software module is developed by SciDM (Scientific Data Management) in 1998-2015
//// 
//// This program is free software; you can redistribute, reuse,
//// or modify it with no restriction, under the terms of the MIT License.
//// 
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//// 
//// For any questions please contact Denis Kaznadzey at dkaznadzey@yahoo.com
//////////////////////////////////////////////////////////////////////////////

#define __common_str_cpp__
#include "common_str.h"

const char* TRUE_STR = "TRUE";
const char* FALSE_STR = "FALSE";
const char* YES_STR = "YES";
const char* NO_STR = "NO";
const char* Y_STR = "Y";
const char* N_STR = "N";
const char* T_STR = "T";
const char* F_STR = "F";
const char* EMPTY_STR = "";
const char* SPACE_STR = " ";
const char* SP2_STR = "  ";
const char* SP4_STR = "    ";
const char* SP6_STR = "      ";
const char* SP8_STR = "        ";
const char* SP10_STR = "          ";
const char* SP12_STR = "            ";
const char* SP14_STR = "              ";
const char* SP16_STR = "                ";
const char* ZERO_STR = "0";
const char* ONE_STR = "1";
const char* MINUS_ONE_STR = "-1";
const char* TWO_STR = "2";
const char* MINUS_TWO_STR = "-2";
const char* THREE_STR = "3";
const char* MINUS_THREE_STR = "-3";
const char* NAME_STR = "name";
const char* NUMBER_STR = "number";
const char* FILENAME_STR = "filename";
const char* FILENAMES_STR = "filenames";
const char* FILEMASKS_STR = "filemasks";
const char* INTEGER_STR = "integer";
const char* BOOLEAN_STR = "boolean";
const char* STRING_STR = "string";
const char* TEXT_STR = "text";
const char* FLOAT_STR = "float";
const char* DOUBLE_STR = "double";
const char* OBJNAME_STR = "object_name";
const char* UNKNOWN_STR = "unknown";

#include <cstring>
#include "portability.h"

const char* inverse_bs (const char* bs)
{
    if (bs == TRUE_STR || strcasecmp (bs, TRUE_STR) == 0)
        return FALSE_STR;
    else if (bs == FALSE_STR || strcasecmp (bs, FALSE_STR) == 0)
        return TRUE_STR;
    else
        return bs;
}

static const char* true_strs  [] = {TRUE_STR, YES_STR, T_STR, Y_STR};
static const char* false_strs [] = {FALSE_STR, NO_STR, F_STR, N_STR};

bool bool_eval (const char* bs, bool* value)
{
    int idx;
    for (idx = 0; idx < sizeof (true_strs) / sizeof (const char*); idx ++)
    {
        if (strcasecmp (true_strs [idx], bs) == 0)
        {
            if (value) *value = true;
            return true;
        }
    }
    for (idx = 0; idx < sizeof (false_strs) / sizeof (const char*); idx ++)
    {
        if (strcasecmp (false_strs [idx], bs) == 0)
        {
            if (value) *value = false;
            return true;
        }
    }
    return false;
}

const char* bool_str (bool val)
{
    if (val) return TRUE_STR;
    else return FALSE_STR;
}

