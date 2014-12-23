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


class TesteGraphGL: public GWin2d {

    public:
        TesteGraphGL();
        void repaint_proc(void);
        void button_proc(Button  button , double  xpos , double  ypos );
        
        
        TesteGraphGL(const TesteGraphGL& orig);
        virtual ~TesteGraphGL();
        private:

};

#endif	/* TESTEGRAPHGL_H */

