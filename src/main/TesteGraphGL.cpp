/* 
 * File:   TesteGraphGL.cpp
 * Author: edsonlan
 * 
 * Created on December 11, 2014, 4:28 PM
 */

#include <FL/Fl.H>

#include "TesteGraphGL.h"
#include "RealVector.h"
#include "CoordsArray.h"

TesteGraphGL::TesteGraphGL() : GWin2d(Frame::root(), Panel::HEIGHT_ONE_ROW,
GUI::default_x + GUI::default_width + 100, GUI::default_y,
GUI::default_height, GUI::default_width) {


    Bounds2d mybounds;
    mybounds.xmin = 0.0;
    mybounds.xmax = 1.0;
    mybounds.ymin = 0.0;
    mybounds.ymax = 1.0;



    bounds(mybounds);
    
//    button_proc(BUTTON_PRESS, 0.5,0.5);

}

TesteGraphGL::TesteGraphGL(const TesteGraphGL& orig) : GWin2d(Frame::root(), Panel::HEIGHT_ONE_ROW,
GUI::default_x + GUI::default_width + 100, GUI::default_y,
GUI::default_height, GUI::default_width) {
}

void TesteGraphGL::repaint_proc(void) {
    
    

    
    std::vector<std::vector<RealVector> > b;
    RpNumerics::physics->boundary()->physical_boundary(b);
    
    cout<<"chamando repaint"<<endl;

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

//        RpnSolution* solution = RpNumerics::viewList->get(i);
//
//        for (int j = 0; j < solution->getCoords()->size() / 2; j++) {
//            double coords[4];
//            coords[0] = solution->getCoords()->at(2 * j)(0);
//            coords[1] = solution->getCoords()->at(2 * j)(1);
//            coords[2] = solution->getCoords()->at(2 * j + 1)(0);
//            coords[3] = solution->getCoords()->at(2 * j + 1)(1);
//
//            CoordsArray coordsArray;
//            coordsArray.init(2, 2, coords);
//            prims().polyline(coordsArray);
//


        }


}


    

      
  


void TesteGraphGL::button_proc(Button button, double xpos, double ypos) {
    if (button & BUTTON_PRESS) {


//        RealVector inDC(2);
//
//        inDC.component(0) = xpos;
//        inDC.component(1) = ypos;
//
//
//        ThreePhaseFlowEquationFunctionLevelCurve * functionCurve = new ThreePhaseFlowEquationFunctionLevelCurve(RpNumerics::physics->flux(), RpNumerics::physics->gridvalues());
//
//        EquationLevelSolution *levelcurve = new EquationLevelSolution(inDC, 1, functionCurve);
//
//        levelcurve->calc();
//        
//        RpNumerics::viewList->add(levelcurve);
//        
//        std::cout << "Chamando DC " << xpos << "  " << ypos << std::endl;

       
        repaint_proc();

    }


//    cout << "tamanho da lista  : " << RpNumerics::viewList->getSize() << endl;

}




TesteGraphGL::~TesteGraphGL() {
    std::cout << "Chamando TesteGraphGL" << std::endl;

}

