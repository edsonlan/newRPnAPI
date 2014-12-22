#ifdef __GNUC__
#pragma implementation
#endif
#include "curvegen.h"
#include <math.h>

#define EPSILON (1e-5)

CurveGenerator::CurveGenerator(float diff,float sample,
                               float X,float Y,
                               float W,float H,
			       float* x,float* y,
			       int x_stride,int y_stride) : px_(x),py_(y),
     x_stride_(x_stride),y_stride_(y_stride),
     diff_inc_(diff<sample?diff:sample),
     sample_inc_(sample),
     xmin_(X), xmax_(X+W),
     ymin_(Y), ymax_(Y+H),
     w_(W), h_(H),
     min_iter_(long(sample/diff))
{
}

CurveGenerator::~CurveGenerator(void)
{
}

void CurveGenerator::setup_tangent(float s)
{
    if (fabs(s)<1)
    {
       cdy_=diff_inc_;
       cdx_=s*cdy_;
    }
    else
    {
       cdx_=diff_inc_;
       cdy_=cdx_/s;
    }
    if ( (ldy_>0 && cdy_<0) || (ldy_<0 && cdy_>0) )
    {
        cdx_=-cdx_;
        cdy_=-cdy_;
    }
}

int CurveGenerator::sample(void)
{
    float _slope;
    long num_iter=0;
    float last_sampled_slope=ldx_/ldy_;

    lx_=cx_;ly_=cy_;
    ldy_=cdy_;
  
    for (;;)
    {
        _slope=slope(cx_,cy_);
        //slope == dx/dy em (cx_,cy_)

        if (cy_&&(cdy_>0))
            if (fabs(_slope-(cx_/cy_))<EPSILON) return 1;
            //retorna 1 se houver um choque.
            //A busca por choques so ocorre em caracteristicas rastreadas 
            //de baixo para cima.
            //Returns 1 if there is a shock.
            //Search for shocks is done only on characteriscts traces
            //from the bottom up.

        setup_tangent(_slope);

        cx_+=cdx_; cy_+=cdy_;

        if (num_iter++>=min_iter_)
        {
            if (fabs(last_sampled_slope-_slope)>EPSILON)
               break;
        }
/*
Apos as primeiras min_iter_ iterations, comeca a verificar se houve 
variacao significativa na declividade da caracteristica. 
Se houve, retorna o ponto de amostragem. Senao, continua.

After the first min_iter_ iterations, starts checking whether there
has been significant variation of the characteristic's slope.
If so, returns the sample point. If not, continues tracing.
*/
        if (!in_boundings()) break;
        //retorna se o rastreamento sair da janela
        //Breaks if the trace goes outside the viewport
    }
    return 0;
}

int CurveGenerator::in_boundings(void)
{
   return (cx_>xmin_ && cy_>ymin_ && cx_<xmax_ && cy_<ymax_);
}
   float *px,*py;
   int n_points_,x_stride_,y_stride_;

void CurveGenerator::write_point(float x,float y)
{
    px_[n_points_*x_stride_]=x;
    py_[n_points_*y_stride_]=y;
    n_points_++;
}

int CurveGenerator::trace(void)
{

    float initial_slope=slope(ix_,iy_);
    cx_=ix_;cy_=iy_;
    ldx_=idx_;ldy_=idy_;

    setup_tangent(initial_slope);
    //gera o vetor de incremento inicial
    //initial increment vector

    ldx_=cdx_;ldy_=cdy_;

    n_points(0);
    write_point(ix_,iy_);
    do
    {
        if(sample()) return 1; 
	//retorna 1 se ocorreu um choque
	//returns 1 if a shock occurs
        write_point(cx_,cy_);
    }while(in_boundings());
    //continua ate a caracteristica sair da janela
    //continues until the characteristic goes outside the window
    return 0;
}

