// MipsPRO compiler defines _BOOL macro if bool implementation is in use
#if !( defined(__GNUC__) || defined(_BOOL) )

#include <iostream>
#include <iomanip.h>
#include <string.h>
#include "bool.h"

const bool true(1);
const bool false(0);

istream &operator>>(istream &is, bool &b)
{
	char c;
#ifdef ELABORATE
	char buf[6];
#endif /* ELABORATE */
	bool b_read;

	is >> c; // skip whitespace
	if (!is) return is;

	switch (c) {
	case '1':
		b_read = true;
		break;
	case '0':
		b_read = false;
		break;
#ifdef ELABORATE
	case 't':
		b_read = true;
		is.putback(c);
		is >> setw(5) >> buf;
		if (strlen(buf) > 1 && strcmp(buf, "true"))
			is.clear(ios::badbit);
		break;
	case 'T':
		b_read = true;
		is.putback(c);
		is >> setw(5) >> buf;
		if (strlen(buf) > 1 && strcmp(buf, "TRUE")
				    && strcmp(buf, "True"))
			is.clear(ios::badbit);
		break;
	case 'f':
		b_read = false;
		is.putback(c);
		is >> setw(6) >> buf;
		if (strlen(buf) > 1 && strcmp(buf, "false"))
			is.clear(ios::badbit);
		break;
	case 'F':
		b_read = false;
		is.putback(c);
		is >> setw(6) >> buf;
		if (strlen(buf) > 1 && strcmp(buf, "FALSE")
				    && strcmp(buf, "False"))
			is.clear(ios::badbit);
		break;
	case 'y':
		b_read = true;
		is.putback(c);
		is >> setw(4) >> buf;
		if (strlen(buf) > 1 && strcmp(buf, "yes"))
			is.clear(ios::badbit);
		break;
	case 'Y':
		b_read = true;
		is.putback(c);
		is >> setw(4) >> buf;
		if (strlen(buf) > 1 && strcmp(buf, "YES")
				    && strcmp(buf, "Yes"))
			is.clear(ios::badbit);
		break;
	case 'n':
		b_read = false;
		is.putback(c);
		is >> setw(3) >> buf;
		if (strlen(buf) > 1 && strcmp(buf, "no"))
			is.clear(ios::badbit);
		break;
	case 'N':
		b_read = false;
		is.putback(c);
		is >> setw(3) >> buf;
		if (strlen(buf) > 1 && strcmp(buf, "NO")
				    && strcmp(buf, "No"))
			is.clear(ios::badbit);
		break;
#endif /* ELABORATE */
	default:
		is.clear(ios::badbit);
		break;
	}

	if (is) b = b_read;
	return is;
}

ostream &operator<<(ostream &os, const bool &b)
{
	return os << (b ? "1" : "0");
}

#endif /* !( defined(__GNUC__) || defined(_BOOL) ) */
