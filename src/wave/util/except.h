#ifndef _except_h
#define _except_h

static class ExceptInit {
public:
	ExceptInit(void);
	~ExceptInit(void);

private:
	static int init_count;
	static void (*old_new_handler)();
	static void (*old_terminate_handler)();
} except_init;

// NOTE:  The __SUNPRO_CC == 0x510 compiler (e.g., at IMPA as of 10/2002)
// has exceptions, but it also defines an incompatible exception structure
// in math.h.  It does have the header file rw/math.h (from Rogue Wave) that
// uses "#define exception math_exception" to try to get around this problem,
// but to take advantage of this would require changing every math.h in the
// code to rw/math.h.

#if (defined(__GNUC__) \
	&& ((__GNUC__ == 2 && __GNUC_MINOR__ >= 91) || __GNUC__ >= 3)) \
	|| (defined(__SUNPRO_CC) && __SUNPRO_CC >= 0x520)
#define except_h_HAS_EXCEPTIONS
#else
#undef except_h_HAS_EXCEPTIONS
#endif

#ifdef except_h_HAS_EXCEPTIONS

#include <exception>

#if defined(__GNUC__) && __GNUC__ >= 3
using std::exception;
#endif

#define THROW(except)	throw except

#else /* except_h_HAS_EXCEPTIONS */

#include <stdlib.h>
#include <iostream>

#if defined(__SUNPRO_CC) && __SUNPRO_CC <= 0x420
// struct exception is defined in math.h
#include <math.h>
#else
class exception {};
#endif

#ifdef __STDC__

#define THROW(except) \
 (cerr << "Run-time exception error; current exception: " \
  #except "\n\tNo handler for exception." << endl, abort(), 0)

#else // SUN C++ 1.0 does not have an ANSI preprocessor!?

#define THROW(except) \
 (cerr << "Run-time exception error:\n\tNo handler for exception." << endl, \
  abort(), 0)

#endif /* __STDC__ */

#endif /* except_h_HAS_EXCEPTIONS */

#endif /* _except_h */
