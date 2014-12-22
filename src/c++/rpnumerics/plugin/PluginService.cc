#include "PluginService.h"

PluginService::PluginService(const string libFileName) {


    libFileName_ = new string(libFileName);

    pluginlib_ = dlopen(libFileName_->c_str(), RTLD_LAZY);
     if (!pluginlib_) {
        cerr << "Cannot load library: " << dlerror() << '\n';
     
    }

}

PluginService::~PluginService() {
    delete libFileName_;
}


RpnPlugin * PluginService::load(const string createMethod) {

    create_t* create_plugin = (create_t*) dlsym(pluginlib_, createMethod.c_str());
    return create_plugin();
}

void PluginService::unload(RpnPlugin * plugin,const string destroyMethod) {

    destroy_t* destroy_plugin = (destroy_t*) dlsym(pluginlib_, destroyMethod.c_str());
    destroy_plugin(plugin);
    dlclose(pluginlib_);

}
