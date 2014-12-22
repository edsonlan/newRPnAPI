#ifndef _bool_h
#define _bool_h

// "bool" is already defined by the GNU C++ compiler, version 2.7.2
// MipsPRO compiler defines _BOOL macro if bool implementation is in use
#if !(defined(__GNUC__) || defined(_BOOL))

class istream;
class ostream;


class bool {
public:
	bool(void) : value(0) {}
	bool(int value_) : value(value_ ? 1 : 0) {}
	operator int(void) const { return value; }

private:
	short value;
};

extern const bool true;
extern const bool false;

istream &operator>>(istream &is, bool &b);
ostream &operator<<(ostream &os, const bool &b);

#endif /* __GNUC__ || _BOOL */

#endif /* _bool_h */
