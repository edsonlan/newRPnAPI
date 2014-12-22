//#ifdef __GNUC__
//#pragma implementation
//#endif
#include <stddef.h>
#include <math.h>
#include "interpol.h"

Interpolator::Interpolator()
    :n_points_(0),
     domain_(NULL),
     image_(NULL),
     slopes_(NULL)
{
}

Interpolator::Interpolator(int s,float *d,float *i)
    :n_points_(s),
     domain_(d),
     image_(i),
     slopes_(NULL)
{
   init();
}

Interpolator::~Interpolator()
{
   if (slopes_) delete slopes_;
}

void Interpolator::load(int s,float *d,float *i)
{
   n_points_=s;		
   domain_=d;//new float[Size];
   image_=i;//new float[Size];
   init();
}

void Interpolator::init(Method m)
{
   int i;
   if (slopes_) delete slopes_;
   switch (method_=m)
   {
       case Linear:
       {
           min_=domain_[0];max_=domain_[n_points_-1];
           slopes_=new float[n_points_-1];
           for (i=1;i<n_points_;i++)
           {
               slopes_[i-1]=(image_[i]-image_[i-1])/(domain_[i]-domain_[i-1]);
           }
       }
       break;
   }
}

int Interpolator::is_invalid(float t)
{
   return t<min_?-1:t>max_?1:0;
}

float Interpolator::eval(float t)
{
   int i;
   if (t<=domain_[0]) return image_[0];
   if (t>=domain_[n_points_-1]) return image_[n_points_-1];
   switch (method_)
   {
       case Linear:
       {
           for (i=0;i<(n_points_-1);i++)
               if (domain_[i+1]>t)
                   return (image_[i]+slopes_[i]*(t-domain_[i]));
       }
   }
   return 0.;
}

