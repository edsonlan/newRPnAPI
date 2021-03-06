/* 
 * File:   TesteGraphGL.cpp
 * Author: edsonlan
 * 
 * Created on December 11, 2014, 4:28 PM
 */

#include <FL/Fl.H>

#include "GraphGLPlotter.h"
#include "RealVector.h"
#include "CoordsArray.h"
#include "RPnResult.h"


#include "LevelConfiguration.h"
#include "EquationLevelMethod.h"

GraphGLPlotter::GraphGLPlotter() : GWin2d(Frame::root(), Panel::HEIGHT_ONE_ROW,
GUI::default_x + GUI::default_width + 100, GUI::default_y,
GUI::default_height, GUI::default_width) {


    Bounds2d mybounds;
    mybounds.xmin = 0.0;
    mybounds.xmax = 1.0;
    mybounds.ymin = 0.0;
    mybounds.ymax = 1.0;



    bounds(mybounds);



}

GraphGLPlotter::GraphGLPlotter(const GraphGLPlotter& orig) : GWin2d(Frame::root(), Panel::HEIGHT_ONE_ROW,
GUI::default_x + GUI::default_width + 100, GUI::default_y,
GUI::default_height, GUI::default_width) {
}

void GraphGLPlotter::repaint_proc(void) {

    std::vector<std::vector<RealVector> > b;
    RpNumerics::physics->boundary()->physical_boundary(b);
    
    for (int i = 0; i < b.size(); i++) {
        if (b[i].size() > 1) {

            double coords[4];

            coords[0] = b[i][0].component(0);
            coords[1] = b[i][0].component(1);
            coords[2] = b[i][1].component(0);
            coords[3] = b[i][1].component(1);

            CoordsArray coordsArray;

            coordsArray.init(2, 3, coords);

            prims().polyline(coordsArray);

        }
    }

    for (int i = 0; i < RpNumerics::viewList->getSize(); i++) {


        RPnResult * curve = RpNumerics::viewList->get(i);
        
        const vector<RealVector> data = curve->getData().getElements();

        if (curve != 0) {

            for (int i = 0; i < data.size()/2; i++) {
                
                const RealVector p1 = curve->getData().getElement(2*i);

                const RealVector p2 = curve->getData().getElement(2*i+1);

                double coords[4];
                coords[0] = p1.component(0);
                coords[1] = p1.component(1);
                coords[2] = p2.component(0);
                coords[3] = p2.component(1);

                CoordsArray coordsArray;
                coordsArray.init(2, 2, coords);
                prims().polyline(coordsArray);

            }


        }
        
        else {
            cout <<"Curve is null"<<endl;
        }

    }

}

void GraphGLPlotter::button_proc(Button button, double xpos, double ypos) {
    if (button & BUTTON_PRESS) {


        RealVector inDC(2);

        inDC.component(0) = xpos;
        inDC.component(1) = ypos;

        LevelConfiguration config(inDC, 1);

        EquationLevelMethod levelcurve(config);

        RPnResult * curve = levelcurve.calc();

        RpNumerics::viewList->add(curve);

        repaint_proc();
    }




}

GraphGLPlotter::~GraphGLPlotter() {
    
    RpNumerics::clean();


}

