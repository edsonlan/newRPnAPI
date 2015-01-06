/* 
 * File:   main.cpp
 * Author: edsonlan
 *
 * Created on December 23, 2014, 12:21 PM
 */

#include <cstdlib>


#include "GraphGLPlotter.h"
#include "Collision.h"
#include "CollisionException.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    GUI::init(&argc,argv);
    RpNumerics::init();
    GraphGLPlotter t;
    
    
//    Collision c;
//    
//    try{
//        const RealVector& Point = c.getPoint();        
//    }
//    catch (CollisionException &e){
//        cout <<"Entrando em catch"<<endl;
//        cout<<e.what()<<endl;
// 
//    }
    
    
    GUI::start_loop();
    
    
    return 0;
}

