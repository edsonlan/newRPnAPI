#ifndef _interpol_h
#define _interpol_h
//#ifdef __GNUC__
//#pragma interface
//#endif

class Interpolator{

public:

   enum Method {Linear};

   int n_points_;
   float *domain_,*image_;

   float min_,max_;
   //Minimo e maximo de domain_ 

private:

   Method method_;

   float *slopes_;
   //Guarda os coeficientes angulares para interpolacao linear
   //Keeps the slope values for linear interpolation

public:

   Interpolator();
   Interpolator(int,float *,float *);
   //IMPORTANTE: este objeto nao mantem copias internas dos arrays, 
   //apenas referencias
   //IMPORTANT: this object does not keep copies of the arrays,
   //only references.

   ~Interpolator();

   void load(int,float *,float *);

   void init(Method m=Linear);
   //Pre-calcula os coeficientes angulares
   //Pre-calculates all linear coeficients

   int is_invalid(float);
   //Retorna 0 se o parameter estiver no intervalo de interpolacao,
   //-1 se estiver abaixo e 1 se estiver acima.
   //Returns 0 if the parameter is inside the interpolation interval,
   //-1 if it's below, and 1 if it's above.

   float eval(float);
   //Retorna o valor interpolado
   //Returns de interpolated value
};

#endif

