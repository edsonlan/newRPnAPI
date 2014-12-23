#ifndef _errorhandler_h
#define _errorhandler_h
#ifdef __GNUC__
#pragma interface
#endif

typedef void (*one_arg_error_handler_t)(const char*);
typedef void (*two_arg_error_handler_t)(const char*, const char*);

extern void default_one_arg_error_handler(const char*);
extern void default_two_arg_error_handler(const char*, const char*);

extern two_arg_error_handler_t lib_error_handler;

extern two_arg_error_handler_t
       set_lib_error_handler(two_arg_error_handler_t f);

#endif /* _errorhandler_h */
