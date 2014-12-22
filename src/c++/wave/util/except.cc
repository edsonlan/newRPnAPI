#include <stdlib.h>
#include <new>
#include <iostream>
#include "except.h"
#if (__GNUC__ >= 3)
using std::set_terminate;
#endif

using namespace std;

int ExceptInit::init_count = 0;
void (*ExceptInit::old_new_handler)() = 0;
void (*ExceptInit::old_terminate_handler)() = 0;

static void
except_init_new_handler(void)
{
	cerr << "Run-time exception error; current exception: bad_alloc\n";
	cerr << "\tNo handler for exception." << endl;
	abort();
}

#ifdef except_h_HAS_EXCEPTIONS
static void
except_init_terminate_handler(void)
{
	cerr << "Run-time exception error; unknown exception\n";
	cerr << "\tNo handler for exception." << endl;
	abort();
}
#endif /* except_h_HAS_EXCEPTIONS */

ExceptInit::ExceptInit(void)
{
	if (++init_count == 1) {
		old_new_handler = set_new_handler(&except_init_new_handler);
#ifdef except_h_HAS_EXCEPTIONS
		old_terminate_handler
			= set_terminate(&except_init_terminate_handler);
#endif /* except_h_HAS_EXCEPTIONS */
	}
}

ExceptInit::~ExceptInit(void)
{
	if (--init_count == 0) {
		(void) set_new_handler(old_new_handler);
#ifdef except_h_HAS_EXCEPTIONS
		(void) set_terminate(old_terminate_handler);
#endif /* except_h_HAS_EXCEPTIONS */
	}
}
