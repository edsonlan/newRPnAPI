//#ifdef __GNUC__
//#pragma implementation
//#endif
#include "plotter.h"
#include <stdio.h>
#include <iostream>
#include <string.h>

#define EPSILON (1e-5)

RiemannPlotter::RiemannPlotter(float xmin,float xmax,float ymax,
                               float distance,
                               int s, float *d,float *i,
                               void (*plot)(int),
			       float* x,float* y,int x_stride,int y_stride,
                               float diff,float sample)
               :CurveGenerator(diff,sample,
                               xmin,0,xmax-xmin,ymax,
			       x,y,x_stride,y_stride),
                interpol_(s,d,i),
		distance_(distance),
		plot_(plot)
{

    printf("Aqui\n");

//        RiemannPlotter Temp((float)grid[0], (float)grid[n - 1], (float)time,
//                            0.045,
//                            n, (float*)grid, (float*)speed,
//                            0, //PlotInterface,
//			    __RPLOTTER_XVals, __RPLOTTER_YVals,
//                            1, 1,
//                            0.001, 0.05);
}

RiemannPlotter::~RiemannPlotter(void)
{
}

void RiemannPlotter::reset(void)
{
    right_=gap_=0;
    current_x_start_=last_top_=xmin();
    last_right_=ymax();
}

int RiemannPlotter::setup_next(void)
{
    printf("    setup_next:\n");

    if (right_)
    {
	rgap_start_-=rgap_step_;
	if (rgap_start_>rgap_end_)
	{
	    setup_trace(0,rgap_start_,1,1);
	    return 1;
	}
	else {right_=0;last_right_=rgap_end_-rgap_step_;}
    }
    if (gap_)
    {
        gap_start_+=gap_step_;
        if (gap_start_<gap_end_)
        {
	    setup_trace(gap_start_,ymax(),0,-1);
            return 1;
        }
        else {gap_=0;last_top_=gap_end_+gap_step_;}
    }

    if ( !right_ && !gap_ && current_x_start_<xmax() )
    {
        do current_x_start_+=distance_;
        while (fabs(current_x_start_)<EPSILON);
        //Nao e tracada uma caracteristica partindo da origem
        //No characteristic is plotted from the origin
        
        if (current_x_start_<xmax())
        {
	    setup_trace(current_x_start_,0,0,1);
            return 1;
        }

        if (last_top_<xmax())
        {
            gap_=1;
            gap_start_=last_top_;
            gap_end_=xmax();
            gap_step_=distance_;
            if (setup_next()) return 1;
        }

	if (fabs(xmin())<EPSILON && last_right_ > EPSILON )
	{
	    right_=1;
	    rgap_start_=last_right_;
	    rgap_end_=0;
	    rgap_step_=distance_/slope(0,last_right_);
	    return setup_next();
	}
    }
    return 0;
}

void RiemannPlotter::get_interception(void)
{
    if (!n_points()) return;
    if (!gap_)
    {
	if (ymax() <= py(n_points()-1))
	{
	    current_top_=px(n_points()-1);
	    if (current_top_-last_top_ > distance_*2)
	    {
		gap_start_=last_top_;
		gap_end_=current_top_;
		gap_=(int)((gap_end_-gap_start_)/distance_);
		gap_step_=(gap_end_-gap_start_)/(gap_+1);
		gap_end_-=gap_step_/2;
		//evita eventuais curvas coladas
		//prevents eventual overlapping characteristics
	    }
	    last_top_=current_top_;
	}
    }
    else if (fabs(xmin())<EPSILON)
    {
	if (px(n_points()-1) <= xmin() &&
	    py(n_points()-1) > sample_increment())
	{
	    current_right_=py(n_points()-1);
	    if (last_right_-current_right_ > distance_*2)
	    {
		rgap_start_=last_right_;
		rgap_end_=current_right_;
		right_=(int)((rgap_start_-rgap_end_)/(distance_/slope(0,current_right_)));
		rgap_end_=(rgap_end_-rgap_start_)/(right_+1);
		rgap_end_+=rgap_step_/2;
            }
            last_right_=current_right_;
        }
    }
}

void RiemannPlotter::gen_family(std::vector<std::vector<RealVector> > &c){
//    c.clear();

    printf("Will generate characterstics.\n");

    int shock;
    reset();
    while(setup_next())
    {
        printf("    ****\n");
        shock=trace();
	get_interception();
        if (n_points()>1) {
            printf("n = %d\n", n_points());
            std::vector<RealVector> temp;

            for (int i = 0; i < n_points(); i++) {
                RealVector r(2);
                r.component(0) = (double)px(i);
                r.component(1) = (double)py(i);

                temp.push_back(r);
            }

            if (temp.size() > 0) c.push_back(temp);
        }
    }

    return;
}

//void RiemannPlotter::gen_family(void)
//{
//    int shock;
//    reset();
//    while(setup_next())
//    {
//        shock=trace();
//	get_interception();
//        if (n_points()>1) (*plot_)(n_points());
//    }
//}

float RiemannPlotter::slope(float x,float y)
{

   float top=y?h()*(x/y):(x>0?xmax():xmin());
   return interpol_.eval(top);
}
