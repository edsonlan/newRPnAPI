
#include "control_Plotter.h"
#include "RpNumerics.h"
#include "ThreePhaseFlowEquationFunctionLevelCurve.h"
#include "EquationLevelMethod.h"
#include "ThreePhaseFlowPermeabilityLevelCurve.h"


JNIEXPORT void JNICALL Java_control_Plotter_cleanList
(JNIEnv * env, jobject cls) {


    RpNumerics::viewList->clean();

}

/*
 * Class:     control_Plotter
 * Method:    plot
 * Signature: (Ljava/lang/String;II)V
 */
JNIEXPORT void JNICALL Java_control_Plotter_plot
(JNIEnv * env, jobject obj, jstring curveName, jint xDC, jint yDC) {


    const char *command;

    command = env->GetStringUTFChars(curveName, NULL);

    string commandString(command);



    Transform2D * transform = RpNumerics::getTransform();


    RealVector inDC(3);


    inDC(0) = xDC;
    inDC(1) = yDC;
    inDC(2) = 1;

    cout << "Coordenadas do device : " << inDC << endl;


    transform->inverseTransform(inDC);
    
    
    cout<<"Coordenadas do mundo: "<<inDC<<endl;




    if (commandString.compare("param") == 0) {

        ThreePhaseFlowSubPhysics * phy = (ThreePhaseFlowSubPhysics *) RpNumerics::physics;

        phy->mug()->value(1 - inDC(0) - inDC(1));
        phy->muo()->value(inDC(1));
        phy->muw()->value(inDC(0));

        RpNumerics::physics->gridvalues()->clear_computations();
        
        for (int i = 0; i < RpNumerics::viewList->getSize(); i++) {

            RpNumerics::viewList->get(i)->realc();


        }


    } else

        if (commandString.compare("level") == 0) {
            
            LevelConfiguration  config(inDC, 1);
            
            EquationLevelMethod  levelcurve(config);
            
            
            RPnResult * curve =levelcurve.calc();
            
            
            RpNumerics::viewList->add(curve);
            
        



    } else if (commandString.compare("hugoniot") == 0) {
        
        
        
        cout<<"Em hugoniot"<<endl;
        
        
        


    }








}