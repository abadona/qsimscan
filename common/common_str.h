
//////////////////////////////////////////////////////////////////////////////
// This software module is developed by SCIDM team in 1999-2015.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// For any questions please contact SciDM team by email at team@scidm.org
//////////////////////////////////////////////////////////////////////////////

#ifndef __common_str_h__
#define __common_str_h__

#ifndef __common_str_cpp__

extern const char* TRUE_STR;
extern const char* FALSE_STR;
extern const char* YES_STR;
extern const char* NO_STR;
extern const char* Y_STR;
extern const char* N_STR;
extern const char* T_STR;
extern const char* F_STR;
extern const char* EMPTY_STR;
extern const char* SPACE_STR;
extern const char* SP2_STR;
extern const char* SP4_STR;
extern const char* SP6_STR;
extern const char* SP8_STR;
extern const char* SP10_STR;
extern const char* SP12_STR;
extern const char* SP14_STR;
extern const char* SP16_STR;
extern const char* ZERO_STR;
extern const char* ONE_STR;
extern const char* MINUS_ONE_STR;
extern const char* TWO_STR;
extern const char* MINUS_TWO_STR;
extern const char* THREE_STR;
extern const char* MINUS_THREE_STR;
extern const char* NAME_STR;
extern const char* NUMBER_STR;
extern const char* FILENAME_STR;
extern const char* FILENAMES_STR;
extern const char* FILEMASKS_STR;
extern const char* INTEGER_STR;
extern const char* BOOLEAN_STR;
extern const char* STRING_STR;
extern const char* TEXT_STR;
extern const char* FLOAT_STR;
extern const char* DOUBLE_STR;
extern const char* OBJNAME_STR;
extern const char* UNKNOWN_STR;

#endif // __common_str_cpp__

#include <stddef.h>

const char* inverse_bs (const char* bs); // returns string expressing boolean value opposite to passed in
bool bool_eval (const char* bs, bool* value = NULL); // returns True if passed string contains valid stringified boolean value and sets value; otherwice returns false
const char* bool_str (bool val); // returns name for a boolean value

#endif
