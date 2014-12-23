#ifndef _curvegen_h
#define _curvegen_h
//#ifdef __GNUC__
//#pragma interface
//#endif

#include <math.h>

class CurveGenerator{
//rastreia curvas caracteristicas
//traces characterics

   float ix_,iy_;
   //ponto inicial para a proxima curva caracteristica
   //Starting point for the next characteristic

   float idx_,idy_;
   //Vetor diretor em (ix_,iy_)
   //Direction hint vector at (ix_,iy_)


   float *px_,*py_;
   int n_points_,x_stride_,y_stride_;
   //arrays para armazenar os pontos
   //array where curves are stored

   float diff_inc_,sample_inc_;
   //incremento diferencial (para rastreio)
   //e de amostragem (distancia minima entre os pontos a serem obtidos)
   //differential increment (for the rastering) e sampling increment (minimum
   //distance between the points to be sampled)

   float xmin_,xmax_,ymin_,ymax_,w_,h_;
   //dados do viewport (janela)
   //viewport(window) data

   float lx_,ly_,
            //Ultimo ponto de amostragem
            //Last sampled point

            cx_,cy_,
            //Ponto atual de rastreio
            //Current raster position

            cdx_,cdy_,
            //Vetor tangente diferencial em (CurrentX,CurrentY)
            //Tangent differential vector at (CurrentX,CurrentY)

            ldx_,ldy_,
            //Vetor tangente diferencial no ultimo ponto de rastreamento
            //Tangent differential vector at the previous raster position

            min_iter_;
            //Numero minimo de iteracoes do rastreio para que o novo ponto
            //possa estar a distancia SampleIncrement do ultimo ponto de
            //amostragem
            //Min number of iterations for the raster loop to obtain a point
            //that can be at a SampleIncrement distance from the last sampled
            //point.

   void setup_tangent(float);
   //Atualiza cdx_ e cdy_ de acordo com slope() e ldy_
   //Assigns new values to cdx_ and cdy_ based on slope() and ldy_

   int sample(void);
   //Rastreia o proximo ponto de amostragem a partir de cx_,cy_ cdx_ e cdy_
   //Traces the next sample point basend on cx_,cy_ cdx_ and cdy_

   int in_boundings(void);
   //Retorna 0 se (cx,cy) estiver fora do viewport
   //Returns 0 if (cx,cy) is outside the viewport

protected:

   CurveGenerator(float diff,float sample,
	          float xmin,float ymin,
	          float width,float height,
		  float* x,float* y,
		  int x_stride,int y_stride);
   virtual ~CurveGenerator();

   float xmin() const {return xmin_;}
   float xmax() const {return xmax_;}
   float ymin() const {return ymin_;}
   float ymax() const {return ymax_;}
   float w() const {return w_;}
   float h() const {return h_;}

   float sample_increment() const {return sample_inc_;}
   float diff_increment() const {return diff_inc_;}

   int n_points() const {return n_points_;}
   void n_points(int n) {n_points_=n;}
   void write_point(float,float);
   float px(int i) {return px_[i*x_stride_];}
   float py(int i) {return py_[i*y_stride_];}

   void setup_trace(float x,float y,float dx,float dy)
     {ix_=x;iy_=y;idx_=dx;idy_=dy;}

   virtual float slope(float,float)=0;
   //funcao que da a declividade dx/dy em (x,y)
   //depende do interpolador, por isso e virtual
   //function that returns dx/dy at (x,y)
   //depends on the interpolator, and is therefore virtual

public:

   int trace(void);
};
#endif
