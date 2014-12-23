#ifndef _mathutil_h
#define _mathutil_h
#ifdef __GNUC__
#pragma interface
#endif

/*
Copyright (C) 1988 Free Software Foundation
    written by Doug Lea (dl@rocky.oswego.edu)

This file is part of GNU CC.

GNU CC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY.  No author or distributor
accepts responsibility to anyone for the consequences of using it
or for whether it serves any particular purpose or works at all,
unless he says so in writing.  Refer to the GNU CC General Public
License for full details.

Everyone is granted permission to copy, modify and redistribute
GNU CC, but only under the conditions described in the
GNU CC General Public License.   A copy of this license is
supposed to have been given to you along with GNU CC so you
can know your rights and responsibilities.  It should be in a
file named COPYING.  Among other things, the copyright notice
and this notice must be preserved on all copies.
*/

#if ((defined __GNUC__ && __GNUC_MINOR__ >= 96) || (__GNUC__ >= 3) \
	|| (defined(__SUNPRO_CC) && __SUNPRO_CC >= 0x520))
#include <cmath>
using std::abs;
#else
inline double abs(double arg) { return (arg < 0.0)? -arg : arg; }
inline float abs(float arg) { return (arg < 0.0)? -arg : arg; }
inline long abs(long arg) { return (arg < 0)? -arg : arg; }
#endif
inline short abs(short arg) { return (arg < 0)? -arg : arg; }

inline long sqr(long arg) { return arg * arg; }
inline double sqr(double arg) { return arg * arg; }

inline unsigned char min(unsigned char a, unsigned char b) {return (a < b)?a:b;}
inline unsigned short min(unsigned short a, unsigned short b) {return (a < b)?a:b;}
inline int min(int a, int b) {return (a < b)?a:b;}
inline unsigned int min(unsigned int a, unsigned int b) {return (a < b)?a:b;}
inline unsigned long min(unsigned long a, unsigned long b) {return (a < b)?a:b;}
inline float min(float a, float b) {return (a < b)?a:b;}
inline double min(double a, double b) {return (a < b)?a:b;}

inline unsigned char max(unsigned char a, unsigned char b) {return (a > b)?a:b;}
inline unsigned short max(unsigned short a, unsigned short b) {return (a > b)?a:b;}
inline int max(int a, int b) {return (a > b)?a:b;}
inline unsigned int max(unsigned int a, unsigned int b) {return (a > b)?a:b;}
inline unsigned long max(unsigned long a, unsigned long b) {return (a > b)?a:b;}
inline float max(float a, float b) {return (a > b)?a:b;}
inline double max(double a, double b) {return (a > b)?a:b;}

#endif /* _mathutil_h */
