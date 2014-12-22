

#ifndef _RPnPluginManager_H
#define	_RPnPluginManager_H

#include "RpnPlugin.h"
#include "PluginService.h"
#include <sys/types.h>
#include <dirent.h>
#include <dlfcn.h>
#include <iostream>
#include <map>

using std::cout;
using std::cerr;
using std::string;
using std::map;


class RPnPluginManager {
private:
    static map<string, map <string,map<string, string> * >  > *configMap_; // Actual plugin configuration <pluginType,libname< <class,constructorMethod>>>
    static map<string, map <string, map<string, string> * > > *destroyMap_; // Actual plugin configuration <pluginType,libname< <class,destroyrMethod>>>
    static string * pluginDir_;

    
public:

    static void setPluginDir(const string);
    static void configPlugin(const string,const string,const string, const string,const string);
    static RpnPlugin * getPluginInstance(const string);
    static void unload(RpnPlugin *, const string);
    virtual ~RPnPluginManager();

};




#endif	/* _RPnPluginManager_H */

