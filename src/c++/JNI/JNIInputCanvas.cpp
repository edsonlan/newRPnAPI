#include "include/canvas_InputCanvas.h"
#include "Graphics2DOpenGL.h"
#include "RpNumerics.h"
#include "../graphicsmodel/RPnMethod.h"
#include "ThreePhaseFlowEquationFunctionLevelCurve.h"
#include "RarefactionCurve.h"
#include "LSODE.h"
#include "RpSolution.h"

/*
 * Class:     canvas_InputCanvas
 * Method:    paint
 * Signature: (Ljava/awt/Graphics;)V
 */
JNIEXPORT void JNICALL Java_canvas_InputCanvas_paint
(JNIEnv *env, jobject canvas, jobject graphics) {


    Graphics2DOpenGL graphics2DOpenGL(env, canvas);

    graphics2DOpenGL.setCoordSystem(0, 1, 0, 1, 0, 1);


    std::vector<std::vector<RealVector> > b;
    RpNumerics::physics->boundary()->physical_boundary(b);


    for (int i = 0; i < b.size(); i++) {
        if (b[i].size() > 1) {

            for (int j = 0; j < b[i].size(); j++) {

//                transform->viewTransform(b[i].at(j));

            }

            graphics2DOpenGL.drawLine(b[i]);


        }
    }


}


