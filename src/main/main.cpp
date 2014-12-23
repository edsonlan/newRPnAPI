/* 
 * File:   main.cpp
 * Author: edsonlan
 *
 * Created on December 23, 2014, 12:21 PM
 */

#include <cstdlib>

#include "GraphGLPlotter.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    GUI::init(&argc,argv);
    RpNumerics::init();
    GraphGLPlotter t;
    
    GUI::start_loop();
    
    
    return 0;
}

