/*
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) RpnPlugin.h
 */


#ifndef _RPNPLUGIN_H
#define	_RPNPLUGIN_H

class RpnPlugin {
public:

    virtual ~RpnPlugin();

};
//the types of the class factories
typedef RpnPlugin* create_t();
typedef void destroy_t(RpnPlugin*);

inline RpnPlugin::~RpnPlugin() {
}



#endif	/* _RPNPLUGIN_H */

