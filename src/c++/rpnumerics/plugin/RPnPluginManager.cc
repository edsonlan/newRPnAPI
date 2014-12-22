#include "RPnPluginManager.h"

map<string, map <string, map<string, string> * > > * RPnPluginManager::configMap_ = NULL;
map<string, map <string, map<string, string> * > > * RPnPluginManager::destroyMap_ = NULL;

string * RPnPluginManager::pluginDir_ = new string();

RPnPluginManager::~RPnPluginManager() {

    delete pluginDir_;
    delete configMap_;
    delete destroyMap_;

}

void RPnPluginManager::unload(RpnPlugin * plugin, const string pluginType) {

    map<string, map <string, map <string, string> *> >::iterator it;

    it = configMap_->find(pluginType);

    if (it == configMap_->end()) {

        cerr << "Plugin " << pluginType << " not configured yet (UNLOAD) " << "\n";


    } else {

        map <string, map<string, string> * > pluginTypeMap = it->second;

        map <string, map<string, string> * >::iterator pluginMapIterator = pluginTypeMap.begin();

        string dirPath(*pluginDir_);

        dirPath += pluginMapIterator->first;
        PluginService service(dirPath);

        //        cout << "Unload Lib: " << dirPath << "\n";


        // Searching destroy method

        map<string, map <string, map <string, string> *> >::iterator it;

        it = destroyMap_->find(pluginType);

        if (it == destroyMap_->end()) {

            cerr << "Plugin " << pluginType << " has no destroy method configured yet" << "\n";

        } else {

            map <string, map<string, string> * > pluginTypeMap = it->second;

            map <string, map<string, string> * >::iterator pluginMapIterator = pluginTypeMap.begin();

            map <string, string> * classMap = pluginMapIterator->second;

            map <string, string> ::iterator classIterator = classMap->begin();

            service.unload(plugin,classIterator->second);
        }


    }
}

RpnPlugin * RPnPluginManager::getPluginInstance(const string pluginType) {

    map<string, map <string, map <string, string> *> >::iterator it;

    it = configMap_->find(pluginType);

    if (it == configMap_->end()) {

        cerr << "Plugin " << pluginType << " not configured yet (LOAD) " << "\n";
        return NULL;

    } else {

        map <string, map<string, string> * > pluginTypeMap = it->second;

        map <string, map<string, string> * >::iterator pluginMapIterator = pluginTypeMap.begin();

        map <string, string> * classMap = pluginMapIterator->second;

        map <string, string> ::iterator classIterator = classMap->begin();

        string dirPath(*pluginDir_);

        dirPath += pluginMapIterator->first;

        PluginService service(dirPath);

        return service.load(classIterator->second);

    }
}

void RPnPluginManager::setPluginDir(const string pluginDir) {
    delete pluginDir_;
    delete configMap_;
    pluginDir_ = new string(pluginDir);
    configMap_ = new map<string, map <string, map<string, string> * > > ();
    destroyMap_ = new map<string, map <string, map<string, string> * > > ();


}

void RPnPluginManager::configPlugin(const string pluginType, const string libName, const string className, const string constructorMethod,const string destructorMethod) {

    map<string, string> * classMap = new map<string, string > ();

    map<string, string> * classMapDestructor = new map<string, string > ();

    //Setting up constructor map

    classMap->operator[](className) = constructorMethod;

    configMap_->erase(pluginType);

    configMap_->operator[](pluginType)[libName] = classMap;

    //Setting up destructor map

    classMapDestructor->operator[](className) = destructorMethod;

    destroyMap_->erase(pluginType);

    destroyMap_->operator[](pluginType)[libName] = classMapDestructor;

}
