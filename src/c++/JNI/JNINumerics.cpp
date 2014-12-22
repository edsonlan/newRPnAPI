#include "main_Numerics.h"
#include "RpNumerics.h"

/*
 * Class:     main_Numerics
 * Method:    initNumerics
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_main_Numerics_initNumerics
  (JNIEnv * env , jclass cls) {
    
    RpNumerics::init();
    
    
}

/*
 * Class:     main_Numerics
 * Method:    clean
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_main_Numerics_clean
  (JNIEnv * env , jclass cls ){
    
    RpNumerics::clean();
    
}