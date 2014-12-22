/**
 * IMPA - Fluid Dynamics Laboratory
 *
 * RPn Project
 *
 * @(#) JNIHugoniotCurveCalc.cc
 **/



/*!
	
TODO:
	
NOTE : 

@ingroup JNI
 */


#include "canvas_NativeCanvas.h"
#include "jawt_md.h"
#include <iostream>


#include "DoubleMatrix.h"
#include "RpNumerics.h"
#include "ImplicitHugoniotCurve.h"


#include "ShockCurve.h"
#include "Graphics2DOpenGL.h"

#include "CoreyQuadSubPhysics.h"
#include "Iso2DTransform.h"
#include "ThreePhaseFlowMobilityLevelCurve.h"
#include "Segment.h"

#include <vector>

using namespace std;

JNIEXPORT void JNICALL Java_canvas_NativeCanvas_paint
(JNIEnv * env, jobject canvas, jobject graphics) {

    jclass testCanvasClass = env->FindClass("canvas/NativeCanvas");

    jfieldID frameID = env->GetFieldID(testCanvasClass, "id_", "I"); //TODO Util para a pegar a projecao do painel

    jint frameNumber = env->GetIntField(canvas, frameID);

    Graphics2DOpenGL graphics2DOpenGL(env, canvas);

    Transform2D * transform = RpNumerics::getTransform();

    graphics2DOpenGL.setCoordSystem(0, 1, 0, 1, 0, 1);

    std::vector<std::vector<RealVector> > b;
    RpNumerics::physics->boundary()->physical_boundary(b);


    for (int i = 0; i < b.size(); i++) {
        if (b[i].size() > 1) {

            graphics2DOpenGL.drawLine(b[i]);

        }
    }

    cout << "Tamanho da lista: " << RpNumerics::viewList->getSize() << endl;

    for (int i = 0; i < RpNumerics::viewList->getSize(); i++) {

        RPnResult * curve = RpNumerics::viewList->get(i);
        
//        cout<<"quantidade de dados: "<<curve->getData()->size()<<endl;

        if (curve!=0) {
            
            const Data * data = curve->getData() ;
            
            
            for (int i = 0; i < curve->getData()->getElements()->size(); i++) {

                
                
                
                

            
            
//            for (const list<Data*>::iterator it=data->getData()->begin(); it!=data->getData()->end(); ++it) {
                
                Segment * segment = (Segment *) curve->getData()->getElement(i);
                
                
                
                
//                cout<<*segment<<endl;
                
                                Point * p1 = (Point *)segment->getElement(0);
                                Point * p2 = (Point *)segment->getElement(1);
                                
                //
                //            RealVector p1 =curve->getSegments()->at(2*j);
                //            RealVector p2 =curve->getSegments()->at(2*j+1);
                
                                graphics2DOpenGL.drawSegment(p1->getCoord()(0),p1->getCoord()(1),p2->getCoord()(0),p2->getCoord()(1),0,0);
                //
                //
            }
        }



    }


}







