
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
// For any questions please contact SciDM team by email at scidmteam@yahoo.com
//////////////////////////////////////////////////////////////////////////////

#include <sciminmax.h>
#include "align.h"

// #define DEBUG_TRACE

#ifdef DEBUG_TRACE
    #include <iostream>
    #include <iomanip>
#endif
/*
---------------------------------
!!! PERFORMANCE CRITICAL CODE !!!
---------------------------------

inner Y coordinate alignment loops:
align_y_loop_xxx ()

these functions can benefit from P2+ cmovxx instructions
(GCC, EGCS, Intel compilers or asm)

    cmp ...
    cmovxx ...
    addc dir, dir
        or
    cmovxx dir, ALIGN_XXX
*/


#if defined (__x86__)

#if defined (_MSC_VER)

#if !defined (_USE_INTEL_COMPILER)
#error This file code must be compiled with Intel compiler
#endif

int align_y_loop_so (ALIGN_VECT* ap, int *wp, int gip, int gep, int max_w, int len)
{
    //initialise bottom boundary
    int v = - (gip + gep);
    ALIGN_VECT* sent = ap + len - 3;

    __asm {
        mov		eax,	0     // eax <- zero
        mov		ecx,	v     // ecx <- initial gap penalty
        mov		esi,	ap    // esi, edi <- align_vect
        mov		edi,	ap

        align 16
loop___1:
        mov		ebx,	wp            // ebx <- weights for current X symbol
        mov		edx,	[esi].r       // edx <- current 'residue'
        // add		eax,	[ebx+edx*8]   // eax += weight for current y | x pair
        add      eax,    [ebx+edx*4]   // eax += weight for current y | x pair

        xor		ebx,	ebx
        mov		edx,	[esi].h
        cmp		ebx,	ecx			;U-V pairing
        cmovl	ebx,	ecx
        cmp		eax,	edx
        cmovl	eax,	edx
        cmp		ebx,	eax
        cmovl	ebx,	eax

        cmp		ebx,	max_w
        jle		no___upd
        mov		max_w, ebx
no___upd:
        mov		eax,	[edi]ap.w
        mov		[edi]ap.w, ebx

        sub		ebx,	gip

        cmp		edx,	ebx			;U-V pairing
        cmovl	edx,	ebx
        cmp		ecx,	ebx
        cmovl	ecx,	ebx
        mov		ebx,	gep
        sub		edx,	ebx
        sub		ecx,	ebx
        mov		[esi].h, edx

        mov		ebx,	wp
        add		esi,	12
        add		edi,	12
        mov		edx,	[esi].r
        add		eax,	[ebx+edx*4]

        xor		ebx,	ebx
        mov		edx,	[esi].h
        cmp		ebx,	ecx			;U-V pairing
        cmovl	ebx,	ecx
        cmp		eax,	edx
        cmovl	eax,	edx
        cmp		ebx,	eax
        cmovl	ebx,	eax

        cmp		ebx,	max_w
        jle		no___upd_1
        mov		max_w, ebx
no___upd_1:
        mov		eax,	[edi]ap.w
        mov		[edi]ap.w, ebx

        sub		ebx,	gip

        cmp		edx,	ebx			;U-V pairing
        cmovl	edx,	ebx
        cmp		ecx,	ebx
        cmovl	ecx,	ebx
        mov		ebx,	gep
        sub		edx,	ebx
        sub		ecx,	ebx
        mov		[esi].h, edx

        mov		ebx,	wp
        add		esi,	12
        add		edi,	12
        mov		edx,	[esi].r
        add		eax,	[ebx+edx*4]

        xor		ebx,	ebx
        mov		edx,	[esi].h
        cmp		ebx,	ecx			;U-V pairing
        cmovl	ebx,	ecx
        cmp		eax,	edx
        cmovl	eax,	edx
        cmp		ebx,	eax
        cmovl	ebx,	eax

        cmp		ebx,	max_w
        jle		no___upd_2
        mov		max_w, ebx
no___upd_2:
        mov		eax,	[edi]ap.w
        mov		[edi]ap.w, ebx

        sub		ebx,	gip

        cmp		edx,	ebx			;U-V pairing
        cmovl	edx,	ebx
        cmp		ecx,	ebx
        cmovl	ecx,	ebx
        mov		ebx,	gep
        sub		edx,	ebx
        sub		ecx,	ebx
        mov		[esi].h, edx

        mov		ebx,	wp
        add		esi,	12
        add		edi,	12
        mov		edx,	[esi].r
        add		eax,	[ebx+edx*4]

        xor		ebx,	ebx
        mov		edx,	[esi].h
        cmp		ebx,	ecx			;U-V pairing
        cmovl	ebx,	ecx
        cmp		eax,	edx
        cmovl	eax,	edx
        cmp		ebx,	eax
        cmovl	ebx,	eax

        cmp		ebx,	max_w
        jle		no___upd_3
        mov		max_w, ebx
no___upd_3:
        mov		eax,	[edi]ap.w
        mov		[edi]ap.w, ebx

        sub		ebx,	gip

        cmp		edx,	ebx			;U-V pairing
        cmovl	edx,	ebx
        cmp		ecx,	ebx
        cmovl	ecx,	ebx
        mov		ebx,	gep
        sub		edx,	ebx
        sub		ecx,	ebx
        mov		[esi].h, edx

        add		esi,	12
        add		edi,	12

        cmp		esi,	sent
        jl		loop___1
    }
    return max_w;
}

#define SCORING_LOOP_DEFINED
#elif defined (__GNUC__)

#if defined (__x86_64__)

// Struc ALIGN_VECT members access:
// w : +0
// h : +4
// r : +8
//int align_y_loop_so (ALIGN_VECT* ap, int *wp, int gip, int gep, int max_w, int len);

int __attribute__((noinline))  align_y_loop_so (ALIGN_VECT* ap, int *wp, int gip, int gep, int max_w, int len)
{
    // int v = - (gip + gep);
    int v = 0;
    ALIGN_VECT* sent = ap + len - 3;

// "m" (v), "m" (ap), "m" (wp), "m" (max_w), "m" (gip), "m" (gep), "m" (sent)
    asm volatile (
        //"nop\n\t"
        //"nop\n\t"
        "movl		$0, %%eax\n\t"
        "movl		%0, %%ecx\n\t"
        //"movl       %1, %%esi\n\t"
        "movq		%1, %%rsi\n\t"
        //"movl       %1, %%edi\n\t"
        "movq		%1, %%rdi\n\t"
        ".align 16\n\t"
    "loop___1:\n\t"
        //"movl		%2, %%ebx\n\t"
        "movq       %2, %%rbx\n\t"
        //"movl		0x8(%%esi), %%edx\n\t"
        "xorq       %%rdx, %%rdx\n\t"
        "movl       0x8(%%rsi), %%edx\n\t"
        //"addl		(%%ebx, %%edx, 0x8), %%eax\n\t"
        //"addl       (%%rbx, %%rdx, 0x8), %%eax\n\t"  // presumable address of weight correction
        "addl       (%%rbx, %%rdx, 0x4), %%eax\n\t"

        "xorl		%%ebx, %%ebx\n\t"
        //"movl		0x4(%%esi), %%edx\n\t"
        "movl       0x4(%%rsi), %%edx\n\t"
        "cmpl		%%ecx, %%ebx\n\t"			// ;U-V pairing
        "cmovll		%%ecx, %%ebx\n\t"
        "cmpl		%%edx, %%eax\n\t"
        "cmovll		%%edx, %%eax\n\t"
        "cmpl		%%eax, %%ebx\n\t"
        "cmovll		%%eax, %%ebx\n\t"
        "cmpl		%3,	%%ebx\n\t"
        "jle		no___upd\n\t"
        "movl		%%ebx, %3\n\t"
    "no___upd:\n\t"
        //"movl		(%%edi), %%eax\n\t"
        "movl       (%%rdi), %%eax\n\t"
        //"movl		%%ebx, (%%edi)\n\t"
        "movl       %%ebx, (%%rdi)\n\t"

        "subl		%4, %%ebx\n\t"

        "cmpl		%%ebx, %%edx\n\t"			//;U-V pairing
        "cmovll		%%ebx, %%edx\n\t"
        "cmpl		%%ebx, %%ecx\n\t"
        "cmovll		%%ebx, %%ecx\n\t"
        "movl		%5, %%ebx\n\t"
        "subl		%%ebx, %%edx\n\t"
        "subl		%%ebx, %%ecx\n\t"
        //"movl		%%edx, 0x4(%%esi)\n\t"
        "movl       %%edx, 0x4(%%rsi)\n\t"

        //"movl		%2, %%ebx\n\t"
        "movq       %2, %%rbx\n\t"
        //"addl		$12, %%esi\n\t"
        "addq       $12, %%rsi\n\t"
        //"addl		$12, %%edi\n\t"
        "addq       $12, %%rdi\n\t"
        //"movl		0x8(%%esi), %%edx\n\t"
        "xorq       %%rdx, %%rdx\n\t"
        "movl       0x8(%%rsi), %%edx\n\t"
        //"addl		(%%ebx, %%edx, 0x4), %%eax\n\t"
        "addl       (%%rbx, %%rdx, 0x4), %%eax\n\t"

        "xorl		%%ebx, %%ebx\n\t"
        //"movl		0x4(%%esi), %%edx\n\t"
        "movl       0x4(%%rsi), %%edx\n\t"
        "cmpl		%%ecx, %%ebx\n\t" //			;U-V pairing
        "cmovll		%%ecx, %%ebx\n\t"
        "cmpl		%%edx, %%eax\n\t"
        "cmovll		%%edx, %%eax\n\t"
        "cmpl		%%eax, %%ebx\n\t"
        "cmovll		%%eax, %%ebx\n\t"

        "cmpl		%3, %%ebx\n\t"
        "jle		no___upd_1\n\t"
        "movl		%%ebx, %3\n\t"
    "no___upd_1:\n\t"
        //"movl		(%%edi), %%eax\n\t"
        "movl       (%%rdi), %%eax\n\t"
        //"movl		%%ebx, (%%edi)\n\t"
        "movl       %%ebx, (%%rdi)\n\t"

        "subl		%4, %%ebx\n\t"

        "cmpl		%%ebx, %%edx\n\t" //			;U-V pairing
        "cmovll		%%ebx, %%edx\n\t"
        "cmpl		%%ebx, %%ecx\n\t"
        "cmovll		%%ebx, %%ecx\n\t"
        "movl		%5, %%ebx\n\t"
        "subl		%%ebx, %%edx\n\t"
        "subl		%%ebx, %%ecx\n\t"
        //"movl		%%edx, 0x4(%%esi)\n\t"
        "movl       %%edx, 0x4(%%rsi)\n\t"

        //"movl		%2, %%ebx\n\t"
        "movq       %2, %%rbx\n\t"
        //"addl		$12, %%esi\n\t"
        "addq       $12, %%rsi\n\t"
        //"addl		$12, %%edi\n\t"
        "addq       $12, %%rdi\n\t"
        //"movl		0x8(%%esi), %%edx\n\t"
        "xorq       %%rdx, %%rdx\n\t"
        "movl       0x8(%%rsi), %%edx\n\t"
        //"addl		(%%ebx, %%edx, 0x4), %%eax\n\t"
        "addl       (%%rbx, %%rdx, 0x4), %%eax\n\t"

        "xorl		%%ebx,	%%ebx\n\t"
        //"movl		0x4(%%esi), %%edx\n\t"
        "movl       0x4(%%rsi), %%edx\n\t"

        "cmpl		%%ecx, %%ebx\n\t" //			;U-V pairing
        "cmovll		%%ecx, %%ebx\n\t"
        "cmpl		%%edx, %%eax\n\t"
        "cmovll		%%edx, %%eax\n\t"
        "cmpl		%%eax, %%ebx\n\t"
        "cmovll		%%eax, %%ebx\n\t"

        "cmpl		%3, %%ebx\n\t"
        "jle		no___upd_2\n\t"
        "movl		%%ebx, %3\n\t"
    "no___upd_2:\n\t"
        //"movl		(%%edi), %%eax\n\t"
        "movl       (%%rdi), %%eax\n\t"
        //"movl		%%ebx, (%%edi)\n\t"
        "movl       %%ebx, (%%rdi)\n\t"

        "subl		%4, %%ebx\n\t"

        "cmpl		%%ebx, %%edx\n\t" //			;U-V pairing
        "cmovll		%%ebx, %%edx\n\t"
        "cmpl		%%ebx, %%ecx\n\t"
        "cmovll		%%ebx, %%ecx\n\t"
        "movl		%5, %%ebx\n\t"
        "subl		%%ebx, %%edx\n\t"
        "subl		%%ebx, %%ecx\n\t"
        //"movl		%%edx, 0x4(%%esi)\n\t"
        "movl       %%edx, 0x4(%%rsi)\n\t"

        //"movl		%2, %%ebx\n\t"
        "movq       %2, %%rbx\n\t"
        //"addl		$12, %%esi\n\t"
        "addq       $12, %%rsi\n\t"
        //"addl		$12, %%edi\n\t"
        "addq       $12, %%rdi\n\t"
        //"movl		0x8(%%esi), %%edx\n\t"
        "xorq       %%rdx, %%rdx\n\t"
        "movl       0x8(%%rsi), %%edx\n\t"
        //"addl		(%%ebx, %%edx, 0x4), %%eax\n\t"
        "addl       (%%rbx, %%rdx, 0x4), %%eax\n\t"

        "xorl		%%ebx, %%ebx\n\t"
        //"movl		0x4(%%esi), %%edx\n\t"
        "movl       0x4(%%rsi), %%edx\n\t"
        "cmpl		%%ecx, %%ebx\n\t" //			;U-V pairing
        "cmovll		%%ecx, %%ebx\n\t"
        "cmpl		%%edx, %%eax\n\t"
        "cmovll		%%edx, %%eax\n\t"
        "cmpl		%%eax, %%ebx\n\t"
        "cmovll		%%eax, %%ebx\n\t"

        "cmpl		%3, %%ebx\n\t"
        "jle		no___upd_3\n\t"
        "movl		%%ebx, %3\n\t"
    "no___upd_3:\n\t"
        //"movl		(%%edi), %%eax\n\t"
        "movl       (%%rdi), %%eax\n\t"
        //"movl		%%ebx, (%%edi)\n\t"
        "movl       %%ebx, (%%rdi)\n\t"

        "subl		%4, %%ebx\n\t"

        "cmpl		%%ebx, %%edx\n\t"	//		;U-V pairing
        "cmovll		%%ebx, %%edx\n\t"
        "cmpl		%%ebx, %%ecx\n\t"
        "cmovll		%%ebx, %%ecx\n\t"
        "movl		%5, %%ebx\n\t"
        "subl		%%ebx, %%edx\n\t"
        "subl		%%ebx, %%ecx\n\t"
        //"movl		%%edx, 0x4(%%esi)\n\t"
        "movl       %%edx, 0x4(%%rsi)\n\t"

        //"addl		$12, %%esi\n\t"
        "addq       $12, %%rsi\n\t"
        //"addl		$12, %%edi\n\t"
        "addq       $12, %%rdi\n\t"

        //"cmp		%6, %%esi\n\t"
        "cmp        %6, %%rsi\n\t"
        "jl		loop___1\n\t"

        :   //"=m" (max_w)
        :	"m" (v), "m" (ap), "m" (wp), "m" (max_w), "m" (gip), "m" (gep), "m" (sent)
        //:	"%eax", "%ebx", "%ecx", "%edx", "%esi", "%edi"
        :   "%eax", "%rbx", "%ecx", "%edx", "%rsi", "%rdi"
        );
    return max_w;
}

#else  // 32-bit mode - gcc assembler

// Struc ALIGN_VECT members access:
// w : +0
// h : +4
// r : +8
//int align_y_loop_so (ALIGN_VECT* ap, int *wp, int gip, int gep, int max_w, int len);

int __attribute__((noinline)) align_y_loop_so (ALIGN_VECT* ap, int *wp, int gip, int gep, int max_w, int len)
{
    ALIGN_VECT* sent = ap + len - 3;
    // int v = - (gip + gep);
    int v = 0;

// "m" (v), "m" (ap), "m" (wp), "m" (max_w), "m" (gip), "m" (gep), "m" (sent)
    asm volatile (
        "movl       $0, %%eax\n\t"
        "movl       %0, %%ecx\n\t"
        "movl       %1, %%esi\n\t"
        "movl       %1, %%edi\n\t"
        ".align 16\n\t"
    "loop___1:\n\t"
        "movl     %2, %%ebx\n\t"
        "movl     0x8(%%esi), %%edx\n\t"
        //"addl     (%%ebx,%%edx,0x8), %%eax\n\t"
        "addl     (%%ebx, %%edx, 0x4), %%eax\n\t"  // CORRECTION?

        "xorl       %%ebx, %%ebx\n\t"
        "movl     0x4(%%esi), %%edx\n\t"
        "cmpl       %%ecx, %%ebx\n\t"           // ;U-V pairing
        "cmovll     %%ecx, %%ebx\n\t"
        "cmpl       %%edx, %%eax\n\t"
        "cmovll     %%edx, %%eax\n\t"
        "cmpl       %%eax, %%ebx\n\t"
        "cmovll     %%eax, %%ebx\n\t"
        "cmpl       %3, %%ebx\n\t"
        "jle        no___upd\n\t"
        "movl       %%ebx, %3\n\t"
    "no___upd:\n\t"
        "movl     (%%edi), %%eax\n\t"
        "movl     %%ebx, (%%edi)\n\t"

        "subl       %4, %%ebx\n\t"

        "cmpl       %%ebx, %%edx\n\t"           //;U-V pairing
        "cmovll     %%ebx, %%edx\n\t"
        "cmpl       %%ebx, %%ecx\n\t"
        "cmovll     %%ebx, %%ecx\n\t"
        "movl       %5, %%ebx\n\t"
        "subl       %%ebx, %%edx\n\t"
        "subl       %%ebx, %%ecx\n\t"
        "movl     %%edx, 0x4(%%esi)\n\t"

        "movl     %2, %%ebx\n\t"
        "addl     $12, %%esi\n\t"
        "addl     $12, %%edi\n\t"
        "movl     0x8(%%esi), %%edx\n\t"
        "addl     (%%ebx, %%edx,0x4), %%eax\n\t"

        "xorl       %%ebx, %%ebx\n\t"
        "movl     0x4(%%esi), %%edx\n\t"
        "cmpl       %%ecx, %%ebx\n\t" //            ;U-V pairing
        "cmovll     %%ecx, %%ebx\n\t"
        "cmpl       %%edx, %%eax\n\t"
        "cmovll     %%edx, %%eax\n\t"
        "cmpl       %%eax, %%ebx\n\t"
        "cmovll     %%eax, %%ebx\n\t"

        "cmpl       %3, %%ebx\n\t"
        "jle        no___upd_1\n\t"
        "movl       %%ebx, %3\n\t"
    "no___upd_1:\n\t"
        "movl     (%%edi), %%eax\n\t"
        "movl     %%ebx, (%%edi)\n\t"

        "subl       %4, %%ebx\n\t"

        "cmpl       %%ebx, %%edx\n\t" //            ;U-V pairing
        "cmovll     %%ebx, %%edx\n\t"
        "cmpl       %%ebx, %%ecx\n\t"
        "cmovll     %%ebx, %%ecx\n\t"
        "movl       %5, %%ebx\n\t"
        "subl       %%ebx, %%edx\n\t"
        "subl       %%ebx, %%ecx\n\t"
        "movl     %%edx, 0x4(%%esi)\n\t"

        "movl     %2, %%ebx\n\t"
        "addl     $12, %%esi\n\t"
        "addl     $12, %%edi\n\t"
        "movl     0x8(%%esi), %%edx\n\t"
        "addl     (%%ebx, %%edx, 0x4), %%eax\n\t"

        "xorl       %%ebx,  %%ebx\n\t"
        "movl     0x4(%%esi), %%edx\n\t"

        "cmpl       %%ecx, %%ebx\n\t" //            ;U-V pairing
        "cmovll     %%ecx, %%ebx\n\t"
        "cmpl       %%edx, %%eax\n\t"
        "cmovll     %%edx, %%eax\n\t"
        "cmpl       %%eax, %%ebx\n\t"
        "cmovll     %%eax, %%ebx\n\t"

        "cmpl       %3, %%ebx\n\t"
        "jle        no___upd_2\n\t"
        "movl       %%ebx, %3\n\t"
    "no___upd_2:\n\t"
        "movl     (%%edi), %%eax\n\t"
        "movl     %%ebx, (%%edi)\n\t"

        "subl       %4, %%ebx\n\t"

        "cmpl       %%ebx, %%edx\n\t" //            ;U-V pairing
        "cmovll     %%ebx, %%edx\n\t"
        "cmpl       %%ebx, %%ecx\n\t"
        "cmovll     %%ebx, %%ecx\n\t"
        "movl       %5, %%ebx\n\t"
        "subl       %%ebx, %%edx\n\t"
        "subl       %%ebx, %%ecx\n\t"
        "movl     %%edx, 0x4(%%esi)\n\t"

        "movl     %2, %%ebx\n\t"
        "addl     $12, %%esi\n\t"
        "addl     $12, %%edi\n\t"
        "movl     0x8(%%esi), %%edx\n\t"
        "addl     (%%ebx, %%edx, 0x4), %%eax\n\t"

        "xorl       %%ebx, %%ebx\n\t"
        "movl     0x4(%%esi), %%edx\n\t"
        "cmpl       %%ecx, %%ebx\n\t" //            ;U-V pairing
        "cmovll     %%ecx, %%ebx\n\t"
        "cmpl       %%edx, %%eax\n\t"
        "cmovll     %%edx, %%eax\n\t"
        "cmpl       %%eax, %%ebx\n\t"
        "cmovll     %%eax, %%ebx\n\t"

        "cmpl       %3, %%ebx\n\t"
        "jle        no___upd_3\n\t"
        "movl       %%ebx, %3\n\t"
    "no___upd_3:\n\t"
        "movl     (%%edi), %%eax\n\t"
        "movl     %%ebx, (%%edi)\n\t"

        "subl       %4, %%ebx\n\t"

        "cmpl       %%ebx, %%edx\n\t"   //      ;U-V pairing
        "cmovll     %%ebx, %%edx\n\t"
        "cmpl       %%ebx, %%ecx\n\t"
        "cmovll     %%ebx, %%ecx\n\t"
        "movl       %5, %%ebx\n\t"
        "subl       %%ebx, %%edx\n\t"
        "subl       %%ebx, %%ecx\n\t"
        "movl     %%edx, 0x4(%%esi)\n\t"

        "addl     $12, %%esi\n\t"
        "addl     $12, %%edi\n\t"

        "cmp      %6, %%esi\n\t"
        "jl     loop___1\n\t"

        :   //"=m" (max_w)
        :   "m" (v), "m" (ap), "m" (wp), "m" (max_w), "m" (gip), "m" (gep), "m" (sent)
        :   "%eax", "%ebx", "%ecx", "%edx", "%esi", "%edi"
        );
    return max_w;
}


#endif

#define SCORING_LOOP_DEFINED
#endif

#endif

#if !defined (SCORING_LOOP_DEFINED)
#error "Please add scoring loop for this platform!"
#define SCORING_LOOP_DEFINED
// the c-code for max local alignment score calculation goes here
#endif

int align_y_loop_n2a_so (register ALIGN_VECT* ap1, register ALIGN_VECT* ap3, int *wp, int gip, int gep, int max_w, int len)
{
    //initialise bottom boundary
    int prev_w = 0;

    // int v = - (gip + gep);
    int v = 0;

    while (len-- > 0)
    {
        register int w = prev_w + wp[ap1->r];

        //w = max (0, h, v, prev_w + s(x, y))
        w = max_ (w, ap1->h);
        w = max_ (w, v);
        w = max_ (w, 0);

        //save max score
        max_w = max_ (max_w, w);

        //save w[x][y] for next pass
        prev_w = ap3->w;
        ap3->w = w;

        w -= gip;
        ap1->h = max_ (w, ap1->h) - gep;
        v = max_ (w, v) - gep;
        ap1++, ap3++;
    }
    return max_w;
}


/*
member inner loops used for finding alignment
fill backtrace matrix and store max score position
*/
/*
n2n, a2a alignment inner loop
*/
void ALIGN::align_y_loop (register ALIGN_VECT* ap_1, int* wp, int gip, int gep, char* bp, int x, int y, int len)
{
    //initialise bottom boundary
    int prev_w = 0;
    // int v = - (gip + gep);
    int v = 0;

#ifdef DEBUG_TRACE
    std::cout << std::endl;
    if (y > 0)
        std::cout << std::setw (4 * y) << std::left << " ";
#endif

    y += len - 1;
    while (len-- > 0)
    {
        register char dir = ALIGN_DIAG;
        register int w = prev_w + wp[ap_1->r];

        //w = max (0, h, v, prev_w + s(x, y))
        if (w < v)
            w = v, dir = ALIGN_DOWN;

        if (w < ap_1->h)
            w = ap_1->h, dir = ALIGN_LEFT;

        if (w <= 0)
            w = 0, dir = ALIGN_STOP;

#ifdef DEBUG_TRACE
        switch (dir)
        {
            case ALIGN_DOWN: std::cout << "-"; break;
            case ALIGN_LEFT: std::cout << "|"; break;
            case ALIGN_DIAG: std::cout << "\\"; break;
            case ALIGN_STOP: std::cout << " "; break;
        }
        std::cout << std::setw (3) << std::left << w << std::setw (0) << std::flush;
#endif
        //save max score position
        if (w > max_w)
            max_w = w, max_bp = bp, max_x = x, max_y = y - len;

        //save w[x][y] for next pass
        prev_w = ap_1->w;
        ap_1->w = w;

        //h = max (w - gip, h) - gep;
        w -= gip;
        if (w > ap_1->h)
            ap_1->h = w, dir |= ALIGN_HSKIP;
        ap_1->h -= gep;

        //v = max (w - gip, v) - gep;
        if (w > v)
            v = w, dir |= ALIGN_VSKIP;
        v -= gep;

        //save bactrace pointer (4bits / byte)
        *bp++ = dir;
        ap_1++;
    }
}


/*
heterogeneous sequences (n2a) alignment inner loop
*/
void ALIGN::align_y_loop_n2a (register ALIGN_VECT* ap_1, register ALIGN_VECT* ap_3, int* wp, int gip, int gep, char* bp, int x, int y, int len)
{
    //initialise bottom boundary
    int prev_w = 0;
    // int v = - (gip + gep);
    int v = 0;

    y += len - 1;
    while (len-- > 0)
    {
        register char dir = ALIGN_DIAG;
        register int w = prev_w + wp[ap_1->r];

        //w = max (0, h, v, prev_w + s(x, y))
        if (w < v)
            w = v, dir = ALIGN_DOWN;

        if (w < ap_1->h)
            w = ap_1->h, dir = ALIGN_LEFT;

        if (w <= 0)
            w = 0, dir = ALIGN_STOP;

        //save max score position
        if (w > max_w)
            max_w = w, max_bp = bp, max_x = x, max_y = y - len;

        //save w[x][y] for next pass
        prev_w = ap_3->w;
        ap_3->w = w;

        //h = max (w - gip, h) - gep;
        w -= gip;
        if (w > ap_1->h)
            ap_1->h = w, dir |= ALIGN_HSKIP;
        ap_1->h -= gep;

        //v = max (w - gip, v) - gep;
        if (w > v)
            v = w, dir |= ALIGN_VSKIP;
        v -= gep;

        //save bactrace pointer (4bits / byte)
        *bp++ = dir;
        ap_1++, ap_3++;
    }
}


/*
heterogeneous sequences (n2a) alignment inner loop
allows different gap penalties for x3 indels vs. x1 indels
*/
/*
void ALIGN::align_y_loop_n2a_3 (register ALIGN_VECT* ap_1, register ALIGN_VECT* ap_3, int *wp, char *bp, int len)
{
    //initialise bottom boundary
    int prev_w = 0;
    int v = - (gip + gep);

    while (len-- > 0)
    {
        register char dir = ALIGN_DIAG;
        register int w = prev_w + wp[ap_1->r];

        //w = max (0, h, h_3, v, prev_w + s(x, y))
        if (w < v)
            w = v, dir = ALIGN_DOWN;

        if (w < ap_1->h)
            w = ap_1->h, dir = ALIGN_LEFT;

        if (w < ap_3->h)
            w = ap_3->h, dir = ALIGN_LEFT_3;

        if (w <= 0)
            w = 0, dir = ALIGN_STOP;

        //save max score position
        if (w > max_w)
            max_w = w, max_bp = bp;

        //save w[x][y] for next pass
        prev_w = ap_3->w;
        ap_3->w = w;

        //h = max (w - gip, h) - gep;
        if (w - gip > ap_1->h)
            ap_1->h = w - gip, dir |= ALIGN_HSKIP;
        ap_1->h -= gep;

        //h = max (w - gip_3, h) - gep_3;
        w -= gip_3;
        if (w > ap_3->h)
            ap_3->h = w, dir |= ALIGN_HSKIP_3;
        ap_3->h -= gep_3;

        //v = max (w - gip_3, v) - gep;
        if (w > v)
            v = w, dir |= ALIGN_VSKIP;
        v -= gep_3;

        //save bactrace pointer (4bits / byte)
        *bp++ = dir;
        ap_1++, ap_3++;
    }
}
*/

/*
--------------------------------
END OF PERFORMANCE CRITICAL CODE
--------------------------------
*/

ALIGN::ALIGN (WEIGHTS<int, 24>* w, int max_ylen, int max_size)
{
    ALIGN::w = w;
    ALIGN::max_ylen = max_ylen;
    ALIGN::max_size = max_size;

    //allocate running Y-vector and backtrace matrix
    btrmx = (char*) malloc (max_size);
    ap_1 = (ALIGN_VECT*) malloc (max_ylen * sizeof (ALIGN_VECT));

    for (int i = 0; i < 3; i++)
        ap_3[i] = (ALIGN_VECT*) malloc (max_ylen * sizeof (ALIGN_VECT));

    if (!btrmx || !ap_1 || !ap_3[0] || !ap_3[1] || !ap_3[2])
        ERR("Unable to allocate memory for alignment");
}


ALIGN::ALIGN (WEIGHTS<int, 24>* w, int max_ylen)
{
    ALIGN::w = w;
    ALIGN::max_ylen = max_ylen;
    ALIGN::max_size = 0;
    btrmx = NULL;

    //allocate running Y-vector
    ap_1 = (ALIGN_VECT*) malloc (max_ylen * sizeof (ALIGN_VECT));
    for (int i = 0; i < 3; i++)
        ap_3[i] = (ALIGN_VECT*) malloc (max_ylen * sizeof (ALIGN_VECT));

    if (!ap_1 || !ap_3[0] || !ap_3[1] || !ap_3[2])
        ERR("unable to allocate memory for alignment vectors");
}


ALIGN::~ALIGN ()
{
    if (ap_1) free (ap_1);

    if (ap_3[0]) free (ap_3[0]);
    if (ap_3[1]) free (ap_3[1]);
    if (ap_3[2]) free (ap_3[2]);

    if (btrmx) free (btrmx);
}


/*
calculates best local alignment between sequence pair of the same type
returns maximum local alignment score
*/
int ALIGN::align (SEQ& xseq, SEQ& yseq)
{
    int x, y;
    char* bp = btrmx;

    //check if enough memory allocated for band alignment
    if (yseq.len > max_ylen || (max_size > 0 && xseq.len * yseq.len > max_size))
    {
        std::clog << std::endl << "align error: attempting to align sequences longer than declared, xlen = " << xseq.len << ", ylen = "  << yseq.len << std::endl << std::flush;
        return 0;
    }

    int gip = int (w->gip);
    int gep = int (w->gep);
    xstep = 1;
    bstep = yseq.len;
    xref = 0, yref = 0;

    //initialize left boundary
    //unpack Y sequence for faster processing
    for (y = 0; y < yseq.len; y++)
    {
        ap_1[y].w = 0;
        // ap_1[y].h = - (w->gip + w->gep);
        ap_1[y].h = 0;
        ap_1[y].r = yseq.get_code (y);
    }

    //find best local alignment
    max_w = 0, max_bp = btrmx;
//	max_w = 99, max_bp = btrmx;
    if (max_size)
        for (x = 0; x < xseq.len; x++, bp += bstep)
            align_y_loop (ap_1, w->mx[xseq.get_code(x)], gip, gep, bp, x, 0, yseq.len);
    else
        for (x = 0; x < xseq.len; x++, bp += bstep)
            max_w = align_y_loop_so (ap_1, w->mx[xseq.get_code(x)], gip, gep, max_w, yseq.len);

    return max_w;
}


/*
calculates best local alignment between sequence pair of different type using aa weighting
returns maximum local alignment score
*/
int ALIGN::align_na (AA_SEQ& yseq, NA_SEQ& xseq)
{
    int x, y;
    char* bp = btrmx;

    //check if enough memory allocated for band alignment
    if (yseq.len > max_ylen || (max_size > 0 && xseq.len * yseq.len > max_size))
    {
        std::clog << std::endl << "align error: attempting to align sequences longer than declared, xlen = " << xseq.len << ", ylen = "  << yseq.len << std::endl << std::flush;
        return 0;
    }

    int gip = int (w->gip);
    int gep = int (w->gep);
    xstep = 3;
    bstep = yseq.len;
    xref = 0, yref = 0;

    //initialize left boundary
    //unpack Y sequence for faster processing
    for (y = 0; y < yseq.len; y++)
    {
        ap_1[y].w = 0;
        // ap_1[y].h = - (gip + gep);
        ap_1[y].h = 0;
        ap_1[y].r = yseq.get_code (y);
        for (int i = 0; i < 3; i++)
        {
            ap_3[i][y].w = 0;
            //ap_3[i][y].h = - (gip + gep);
            ap_3[i][y].h = 0;
        }
    }

    //find best local alignment
    max_w = 0, max_bp = btrmx;
    if (max_size)
        for (x = 0; x < xseq.len; x++, bp += bstep)
            align_y_loop_n2a (ap_1, ap_3[x%3], w->mx[xseq.get_code(x)], gip, gep, bp, x, 0, yseq.len);
    else
        for (x = 0; x < xseq.len; x++, bp += bstep)
            max_w = align_y_loop_n2a_so (ap_1, ap_3[x%3], w->mx[xseq.get_code(x)], gip, gep, max_w, yseq.len);

    return max_w;
}


/*
calculates best local alignment between sequence pair of the same type
on a diagonal band (len, diag +- width)
returns maximum local alignment score
NOTE: batch xpos, ypos, len should be inside X and Y sequences, width > 0
*/
int ALIGN::align_band (SEQ& xseq, SEQ& yseq, int xpos, int ypos, int len, int width)
{
    char* bp = btrmx;

    //check if enough memory allocated for band alignment
    //and if batch variables are sane
    if (yseq.len > max_ylen || (max_size > 0 && len * width > max_size))
    {
        std::clog << std::endl << "align error: attempting to batch-align sequences longer than declared, ylen = " << yseq.len << ", len = "  << len << ", width = " << width << std::endl << std::flush;
        return 0;
    }

    int inilen = len;
    int gip = int (w->gip);
    int gep = int (w->gep);
    xstep = 1;
    bstep = width * 2;
    xref = xpos, yref = ypos;

    //initialize left boundary
    //unpack Y sequence for faster processing
    int ylast = min_ (yseq.len, yref + len + bstep);
    for (int i = max_ (0, yref - width); i < ylast; i++)
    {
        ap_1[i].w = 0;
        // ap_1[i].h = - (int) (w->gip + w->gep);
        ap_1[i].h = 0;
        ap_1[i].r = yseq.get_code (i);
    }

    //find best local alignment, save backtrace pointers
    max_w = 0, max_bp = btrmx, ypos -= width;

    int y_start, y_end, y_len; //, y;
    int *mx_row;
    char *bp_start;

    while (len-- > 0)
    {
        //clip Y vector versus batch boundaries
        y_start = max_ (0, ypos);
        y_end = min_ (yseq.len, ypos + bstep);
        // y_end = min_ (y_end, yseq.len);
        if (y_end > y_start)
        {
            y_len = y_end - y_start;

            bp_start = bp + y_start - ypos;
            mx_row = w->mx[xseq.get_code (xpos)];

            if (max_size)
                align_y_loop (ap_1 + y_start, mx_row, gip, gep, bp_start, xpos, y_start /* + yref*/, y_len);
            else
                max_w = align_y_loop_so (ap_1 + y_start, mx_row, gip, gep, max_w, y_len);
        }
        xpos++, ypos++, bp += bstep;
    }

    //force diagonal offset for backtrace matrix
    bstep--;

    return max_w;
}


/*
follows backtrace matrix, fills BATCH array, returns number of batches
*/
int ALIGN::backtrace (BATCH *b_ptr, int max_cnt, unsigned width)
{
    char *bp = max_bp;
    int state = *bp & 3;
    int x = max_x, y = max_y;
    int b_len = 1, b_cnt = 0;
    int lowest_y = max_ (0, (int) yref - (int) width);
    BATCH* b_start = b_ptr;

    //backtrace from highest score
    bool done = false;
    do
    {
        // if (semiglobal && !to_first && (*bp & ALIGN_ZERO))
        //    break;

        switch (state)
        {
            //follow v-trace down until ALIGN_VSKIP flag set
            case ALIGN_DOWN:
                y--;
                if (y >= yref)
                {
                    bp--;
                    if (*bp & ALIGN_VSKIP)
                        state = *bp & 3;
                }
                break;
            //follow h-trace left until ALIGN_HSKIP flag set
            case ALIGN_LEFT:
                x--;
                if (x >= xref)
                {
                    bp -= bstep;
                    if (*bp & ALIGN_HSKIP)
                        state = *bp & 3;
                }
                break;
            //follow diagonal until best score is achieved from v-trace or h-trace
            case ALIGN_DIAG:
                bp -= bstep + 1;
                if (x > xref && y > yref && state != (*bp & 3))
                {
                    state = *bp & 3;
                    if (b_cnt == max_cnt)
                        ers << "Not enough space to store all batches (>" << b_cnt-1 << ") of the alignment" << Throw;
                    b_cnt ++;
                    b_ptr->xpos = x;
                    b_ptr->ypos = y;
                    b_ptr->len = b_len;
                    b_ptr ++;
                    b_len = 0;
                }
                b_len ++;
                x --, y --;
                break;
            //end of alignment (w[x][y] was set to 0)
            case ALIGN_STOP:
                done = true;
                break;
        }
    }
    while (!done && x >= xref && y >= lowest_y);
    //if alignment ends at the edge of the matrix we get here
    if (state == ALIGN_DIAG)
    {
        if (b_cnt == max_cnt)
            ers << "Not enough space to store all batches (>" << b_cnt-1 << ") of the alignment" << Throw;
        b_cnt++;
        b_ptr->xpos = x+1;
        b_ptr->ypos = y+1;
        b_ptr->len = b_len-1;
        b_ptr ++;
    }

#if 0
    while (x >= xref + xstep && y > yref && y > 0)
    {
        switch (state)
        {
            //follow v-trace down until ALIGN_VSKIP flag set
            case ALIGN_DOWN:
                bp--;
                if (*bp & ALIGN_VSKIP)
                    state = *bp & 3;
                y--;
                break;

            //follow h-trace left until ALIGN_HSKIP flag set
            case ALIGN_LEFT:
                bp -= bstep;
                if (*bp & ALIGN_HSKIP)
                    state = *bp & 3;
                x--;
                break;

            //follow diagonal until best score is achieved from v-trace or h-trace
            case ALIGN_DIAG:
                bp -= xstep * bstep + 1;

                if (state != (*bp & 3))
                {
                    state = *bp & 3;

                    b_ptr->xpos = x;
                    b_ptr->ypos = y;
                    b_ptr->len = b_len;
                    b_ptr++;

                    if (++b_cnt >= max_cnt)
                    {
                        reverse<BATCH> (b_start, b_cnt);
                        return b_cnt;
                    }

                    b_len = 0;
                }
                x -= xstep, y --, b_len++;
                break;

            //end of alignment (w[x][y] was set to 0)
            case ALIGN_STOP:
                reverse<BATCH> (b_start, b_cnt);
                return b_cnt;
        }
    }

    //if alignment ends at the edge of the matrix we get here
    if (state == ALIGN_DIAG)
    {
        b_ptr->xpos = x;
        b_ptr->ypos = y;
        b_ptr->len = b_len;
        b_cnt++;
    }
#endif
    reverse<BATCH> (b_start, b_cnt);
    return b_cnt;
}
