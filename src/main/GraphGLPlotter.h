/* 
 * File:   TesteGraphGL.h
 * Author: edsonlan
 *
 * Created on December 11, 2014, 4:28 PM
 */

#ifndef TESTEGRAPHGL_H
#define	TESTEGRAPHGL_H
#include "GWin2d.h"
#include "RpNumerics.h"


class GraphGLPlotter: public GWin2d {

    public:
        GraphGLPlotter();
        void repaint_proc(void);
        void button_proc(Button  button , double  xpos , double  ypos );
        
        
        GraphGLPlotter(const GraphGLPlotter& orig);
        virtual ~GraphGLPlotter();
        private:

};

#endif	/* TESTEGRAPHGL_H */

