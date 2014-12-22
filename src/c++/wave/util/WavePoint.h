

#ifndef _WavePoint_H
#define	_WavePoint_H

//#include "Vector2.h"

#include "RealVector.h"

class WavePoint:public RealVector {
    
    public:
        
        WavePoint(int , double);
        ~WavePoint();
        
        double getSigma();
        
        private:
            
            
};


#endif	/* _WavePoint_H */

