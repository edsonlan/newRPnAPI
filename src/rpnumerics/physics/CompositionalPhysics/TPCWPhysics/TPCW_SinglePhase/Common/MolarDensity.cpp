#include "MolarDensity.h"


MolarDensity::MolarDensity(int t, double Pressure) : type(t), 
                                                     P(Pressure), 
                                                     R(8.314472){

    Tcc = 304.10;
    Tcw = 647.3 ;

    Pcc =  73.7*1e5;  // This pressure must not be changed.  The standard pressure unit of this cod is Pascals.
    Pcw = 221.2*1e5;  // This pressure must not be changed.  The standard pressure unit of this cod is Pascals.

    bc     = 0.077796*R*Tcc/Pcc ;
    log_bc = log(bc) ;
    bw     = 0.077796*R*Tcw/Pcw ;
    log_bw = log(bw) ;

    q1 = -0.4347  ;
    q2 = -0.003654;
    omega_w  =  0.344  ;
    kappa_1w = -0.06635;
    omega_c  =  0.23894;
    kappa_1c =  0.04285;

    Delta_G0_12 = 3909.50;
    Delta_G1_12 =   18.90;
    Delta_G0_21 = 1473.60;
    Delta_G1_21 =   16.98;

    // Auxiliar constants for the function Molar Density (in tau12 and tau21)
    Aux_G12 = Delta_G0_12/R;
    Aux_G21 = Delta_G0_21/R;

    // Precision for conditional equal to zero. TODO: Ver o valor que deve ser empregado.
    zero_epsilon = pow(0.1, 14);

    alpha12     = -0.3;  //TODO: Notice that alpha12 is defined as -alpha12

    C0_c = -1.49632e-6;
    C1_c =  4.70560e-9;
    C0_w = -1.07169e-7;
    C1_w = -1.01200e-8;
//    cout<<" CTOR Molar density: "<<this<<endl;

}

MolarDensity::~MolarDensity(){
    
//      cout<<"     DTOR Molar density: "<<this<<endl;
}

// TODO: ESTE CODIGO DEVERIA ESTAR NUMA CLASSE AUXILIAR DE UTILITARIOS

// Cubic root. Original at:
//
//     ftp://ftp.netlib.org/fdlibm/s_cbrt.c
//
// and some #define's from here:
//
//     http://www.netlib.org/fdlibm/fdlibm.h
// 
double MolarDensity::cbrt(double x){
    if      (x > 0.0)  return exp(log(x)/3.0);
    else if (x == 0.0) return 0.0;
    else               return -exp(log(-x)/3.0);

//    #ifdef __LITTLE_ENDIAN
//    #define __HI(x) *(1+(int*)&x)
//    #define __LO(x) *(int*)&x
//    #define __HIp(x) *(1+(int*)x)
//    #define __LOp(x) *(int*)x
//    #else
//    #define __HI(x) *(int*)&x
//    #define __LO(x) *(1+(int*)&x)
//    #define __HIp(x) *(int*)x
//    #define __LOp(x) *(1+(int*)x)
//    #endif

//    #ifdef __STDC__
//    static const unsigned 
//    #else
//    static unsigned 
//    #endif
//    B1 = 715094163, /* B1 = (682-0.03306235651)*2**20 */
//    B2 = 696219795; /* B2 = (664-0.03306235651)*2**20 */

//    #ifdef __STDC__
//    static const double
//    #else
//    static double
//    #endif
//    C =  5.42857142857142815906e-01, /* 19/35     = 0x3FE15F15, 0xF15F15F1 */
//    D = -7.05306122448979611050e-01, /* -864/1225 = 0xBFE691DE, 0x2532C834 */
//    E =  1.41428571428571436819e+00, /* 99/70     = 0x3FF6A0EA, 0x0EA0EA0F */
//    F =  1.60714285714285720630e+00, /* 45/28     = 0x3FF9B6DB, 0x6DB6DB6E */
//    G =  3.57142857142857150787e-01; /* 5/14      = 0x3FD6DB6D, 0xB6DB6DB7 */

//    int hx;
//    double r,s,t=0.0,w;
//    unsigned sign;

//    hx = __HI(x);		/* high word of x */
//////printf("HI(%f) = %d\n", x, hx);
//    sign=hx&0x80000000; 		/* sign= sign(x) */
//    hx  ^=sign;
//    if(hx>=0x7ff00000) return(x+x); /* cbrt(NaN,INF) is itself */
//    if((hx|__LO(x))==0) 
//        return(x);		/* cbrt(0) is itself */

//    __HI(x) = hx;	/* x <- |x| */
//    /* rough cbrt to 5 bits */
//    if(hx<0x00100000) 		/* subnormal number */
//        {__HI(t)=0x43500000; 		/* set t= 2**54 */
//            t*=x; __HI(t)=__HI(t)/3+B2;
//        }
//        else
//            __HI(t)=hx/3+B1;	

//    /* new cbrt to 23 bits, may be implemented in single precision */
//    r=t*t/x;
//    s=C+r*t;
//    t*=G+F/(s+E+D/s);	

//    /* chopped to 20 bits and make it larger than cbrt(x) */ 
//    __LO(t)=0; __HI(t)+=0x00000001;


//    /* one step newton iteration to 53 bits with error less than 0.667 ulps */
//    s=t*t;		/* t*t is exact */
//    r=x/s;
//    w=t+t;
//    r=(r-t)/(w+r);	/* r-s is exact */
//    t=t+t*r;

//    /* retore the sign bit */
//    __HI(t) |= sign;
//    return(t);
}

// Slightly modified from the original at:
//
//    http://www.dreamincode.net/code/snippet2915.htm
//
// The function returns:
//
//     0: There are three real roots, and they are all the same.
//     1: There are three real roots, and they are all different.
//     2: There are one real root and two complex roots (one being the other's conjugate).
//
int MolarDensity::CubicRoots(const double a, const double b, const double c, const double d, 
                             std::complex<double> &r0, 
                             std::complex<double> &r1, 
                             std::complex<double> &r2){

    double ROOTTHREE = 1.73205080756888;

    // Find the discriminant
    double f, g, h;
    f = (3.0 * c / a - pow(b, 2) / pow(a, 2)) / 3.0;
    g = (2.0 * pow(b, 3) / pow(a, 3) - 9.0 * b * c / pow(a, 2) + 27.0 * d / a) / 27.0;
    h = pow(g, 2) / 4.0 + pow(f, 3) / 27.0;

    // Evaluate discriminant
    // TODO: Aqui entra o zero_epsilon, conferir que o valor seja adecuado.
    if ( ( fabs(f) < zero_epsilon) && (fabs(g) < zero_epsilon) && (fabs(h) < zero_epsilon) ){
        // 3 equal roots
        double x;

        // when f, g, and h all equal 0 the roots can be found by the following line
        x = -cbrt(d / a);

        r0 = std::complex<double>(x, 0.0);
        r1 = std::complex<double>(x, 0.0);
        r2 = std::complex<double>(x, 0.0);

        return 0;
    }
    else if (h <= 0.0){
        // 3 real roots
        double q, i, j, k, l, m, n, p;
        // complicated maths making use of the method
        i = pow(pow(g, 2) / 4.0 - h, 0.5);
        j = cbrt(i);
        k = acos(-(g / (2.0 * i)));
        m = cos(k / 3.0);
        n = ROOTTHREE * sin(k / 3.0);
        p = -(b / (3.0 * a));

        r0 = std::complex<double>(2.0 * j * m + p, 0.0);
        r1 = std::complex<double>(-j * (m + n) + p, 0.0);
        r2 = std::complex<double>(-j * (m - n) + p, 0.0);

        return 1;
    }
    else if (h > 0.0){
        // 1 real root and 2 complex roots
        double r, s, t, u, p;
        // complicated maths making use of the method
        r = -(g / 2.0) + pow(h, 0.5);
        s = cbrt(r);
        t = -(g / 2.0) - pow(h, 0.5);
        u = cbrt(t);
        p = -(b / (3.0 * a));

        r0 = std::complex<double>((s + u) + p, 0.0);
        r1 = std::complex<double>(-(s + u) / 2.0 + p, (s - u) * ROOTTHREE / 2.0);
        r2 = std::complex<double>(-(s + u) / 2.0 + p, -(s - u) * ROOTTHREE / 2.0);

        return 2;
    }
}


//
// epsilon_c
//

int  MolarDensity::epsilon_c_jet(const double T, int degree, JetMatrix &epsilon_cj){

    if (degree >= 0){

        double T_red          = T/Tcc;
        double sqrt_Tr        = sqrt(T_red);
        double epsilon_factor = 0.457235 * Tcc / 0.077796; // It represents ac/(bc*R)
        double kappa_0c       = 0.378893 + 1.4897153*omega_c - 0.17131848*omega_c*omega_c + 0.0196554*omega_c*omega_c*omega_c;
        double inv_T          = 1.0/T;
         
        double alpha_c = 1.0 + ( kappa_0c + kappa_1c * (1.0 + sqrt_Tr) * (0.7 - T_red) ) * (1.0 - sqrt_Tr);

        double epsilon_c = epsilon_factor * alpha_c * alpha_c * inv_T;

        //TODO: Fazer as derivadas explicitas daqui em diante...

//        double epsilon_c = ( Tcc * pow((sqrt(T/Tcc)-1.0)*(omega_c*1.4897153 - (omega_c*omega_c)*0.17131848 + (omega_c*omega_c*omega_c)*0.0196554 - kappa_1c*(sqrt(T/Tcc)+1.0) * (T/Tcc-0.7) + 0.378893) - 1.0,2.0)*(0.457235/0.077796))/T;

        epsilon_cj.set(0,epsilon_c) ;
        
        if (degree >= 1){

            double depsilon_c_dT = 1.0/(T*T)*Tcc*pow((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0,2.0)*(-1.647362700496599E16/2.802896292887321E15)-(Tcc*(((kappa_1c*(sqrt(T/Tcc)+1.0))/Tcc+(kappa_1c*1.0/sqrt(T/Tcc)*(T/Tcc-7.0/1.0E1)*(1.0/2.0))/Tcc)*(sqrt(T/Tcc)-1.0)-(1.0/sqrt(T/Tcc)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcc)*((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*(3.294725400993198E16/2.802896292887321E15))/T;
            
            epsilon_cj.set(0, 0, depsilon_c_dT);

            if (degree == 2){
                                
                double d2epsilon_c_dT2 = (Tcc*pow(((kappa_1c*(sqrt(T/Tcc)+1.0))/Tcc+(kappa_1c*1.0/sqrt(T/Tcc)*(T/Tcc-7.0/1.0E1)*(1.0/2.0))/Tcc)*(sqrt(T/Tcc)-1.0)-(1.0/sqrt(T/Tcc)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcc,2.0)*(3.294725400993198E16/2.802896292887321E15))/T+1.0/(T*T*T)*Tcc*pow((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0,2.0)*(3.294725400993198E16/2.802896292887321E15)-(Tcc*((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*((sqrt(T/Tcc)-1.0)*(1.0/(Tcc*Tcc)*kappa_1c*1.0/sqrt(T/Tcc)-1.0/(Tcc*Tcc)*kappa_1c*1.0/pow(T/Tcc,3.0/2.0)*(T/Tcc-7.0/1.0E1)*(1.0/4.0))+1.0/(Tcc*Tcc)*1.0/pow(T/Tcc,3.0/2.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/4.0)+(((kappa_1c*(sqrt(T/Tcc)+1.0))/Tcc+(kappa_1c*1.0/sqrt(T/Tcc)*(T/Tcc-7.0/1.0E1)*(1.0/2.0))/Tcc)*1.0/sqrt(T/Tcc))/Tcc)*(3.294725400993198E16/2.802896292887321E15))/T+1.0/(T*T)*Tcc*(((kappa_1c*(sqrt(T/Tcc)+1.0))/Tcc+(kappa_1c*1.0/sqrt(T/Tcc)*(T/Tcc-7.0/1.0E1)*(1.0/2.0))/Tcc)*(sqrt(T/Tcc)-1.0)-(1.0/sqrt(T/Tcc)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcc)*((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*(6.589450801986396E16/2.802896292887321E15);

                epsilon_cj.set(0, 0, 0, d2epsilon_c_dT2 );

                }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}


//
// epsilon_w
//

int  MolarDensity::epsilon_w_jet(const double T, int degree, JetMatrix &epsilon_wj){
    //std::cout << "Inside epsilon_w_jet: T = " << T << ", degree = " << degree << std::endl;
    if (degree >= 0){

        double T_red          = T/Tcw;
        double sqrt_Tr        = sqrt(T_red);
        double epsilon_factor = 0.457235 * Tcw / 0.077796; // It represents ac/(bc*R)
        double kappa_0w       = 0.378893 + 1.4897153*omega_w - 0.17131848*omega_w*omega_w + 0.0196554*omega_w*omega_w*omega_w;
        double inv_T          = 1.0/T;
         
        double alpha_w = 1.0 + ( kappa_0w + kappa_1w * (1.0 + sqrt_Tr) * (0.7 - T_red) ) * (1.0 - sqrt_Tr);

        double epsilon_w = epsilon_factor * alpha_w * alpha_w * inv_T;

        //TODO: Fazer as derivadas explicitas daqui em diante...

        epsilon_wj.set(0,epsilon_w) ;
 
        if (degree >= 1){

            double depsilon_w_dT = 1.0/(T*T)*Tcw*pow((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0,2.0)*(-1.647362700496599E16/2.802896292887321E15)-(Tcw*(((kappa_1w*(sqrt(T/Tcw)+1.0))/Tcw+(kappa_1w*1.0/sqrt(T/Tcw)*(T/Tcw-7.0/1.0E1)*(1.0/2.0))/Tcw)*(sqrt(T/Tcw)-1.0)-(1.0/sqrt(T/Tcw)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcw)*((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*(3.294725400993198E16/2.802896292887321E15))/T;
 
            epsilon_wj.set(0, 0, depsilon_w_dT);

            if (degree == 2){
                                
                double d2epsilon_w_dT2 = (Tcw*pow(((kappa_1w*(sqrt(T/Tcw)+1.0))/Tcw+(kappa_1w*1.0/sqrt(T/Tcw)*(T/Tcw-7.0/1.0E1)*(1.0/2.0))/Tcw)*(sqrt(T/Tcw)-1.0)-(1.0/sqrt(T/Tcw)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcw,2.0)*(3.294725400993198E16/2.802896292887321E15))/T+1.0/(T*T*T)*Tcw*pow((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0,2.0)*(3.294725400993198E16/2.802896292887321E15)-(Tcw*((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*((sqrt(T/Tcw)-1.0)*(1.0/(Tcw*Tcw)*kappa_1w*1.0/sqrt(T/Tcw)-1.0/(Tcw*Tcw)*kappa_1w*1.0/pow(T/Tcw,3.0/2.0)*(T/Tcw-7.0/1.0E1)*(1.0/4.0))+1.0/(Tcw*Tcw)*1.0/pow(T/Tcw,3.0/2.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/4.0)+(((kappa_1w*(sqrt(T/Tcw)+1.0))/Tcw+(kappa_1w*1.0/sqrt(T/Tcw)*(T/Tcw-7.0/1.0E1)*(1.0/2.0))/Tcw)*1.0/sqrt(T/Tcw))/Tcw)*(3.294725400993198E16/2.802896292887321E15))/T+1.0/(T*T)*Tcw*(((kappa_1w*(sqrt(T/Tcw)+1.0))/Tcw+(kappa_1w*1.0/sqrt(T/Tcw)*(T/Tcw-7.0/1.0E1)*(1.0/2.0))/Tcw)*(sqrt(T/Tcw)-1.0)-(1.0/sqrt(T/Tcw)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcw)*((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*(6.589450801986396E16/2.802896292887321E15);

                epsilon_wj.set(0, 0, 0, d2epsilon_w_dT2 );

                }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}






/*
//
// epsilon_c
//

int  MolarDensity::epsilon_c_jet(const double T, int degree, JetMatrix &epsilon_cj){

    if (degree >= 0){
         
// We can divide this calculation into little pieces in order to verify first this jet.

        double epsilon_c = (Tcc*pow((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0,2.0)*(1.647362700496599E16/2.802896292887321E15))/T;

        epsilon_cj(0,epsilon_c) ;
        
        if (degree >= 1){

            double depsilon_c_dT = 1.0/(T*T)*Tcc*pow((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0,2.0)*(-1.647362700496599E16/2.802896292887321E15)-(Tcc*(((kappa_1c*(sqrt(T/Tcc)+1.0))/Tcc+(kappa_1c*1.0/sqrt(T/Tcc)*(T/Tcc-7.0/1.0E1)*(1.0/2.0))/Tcc)*(sqrt(T/Tcc)-1.0)-(1.0/sqrt(T/Tcc)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcc)*((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*(3.294725400993198E16/2.802896292887321E15))/T;
            
            epsilon_cj(0, 0, depsilon_c_dT);

            if (degree == 2){
                                
                double d2epsilon_c_dT2 = (Tcc*pow(((kappa_1c*(sqrt(T/Tcc)+1.0))/Tcc+(kappa_1c*1.0/sqrt(T/Tcc)*(T/Tcc-7.0/1.0E1)*(1.0/2.0))/Tcc)*(sqrt(T/Tcc)-1.0)-(1.0/sqrt(T/Tcc)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcc,2.0)*(3.294725400993198E16/2.802896292887321E15))/T+1.0/(T*T*T)*Tcc*pow((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0,2.0)*(3.294725400993198E16/2.802896292887321E15)-(Tcc*((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*((sqrt(T/Tcc)-1.0)*(1.0/(Tcc*Tcc)*kappa_1c*1.0/sqrt(T/Tcc)-1.0/(Tcc*Tcc)*kappa_1c*1.0/pow(T/Tcc,3.0/2.0)*(T/Tcc-7.0/1.0E1)*(1.0/4.0))+1.0/(Tcc*Tcc)*1.0/pow(T/Tcc,3.0/2.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/4.0)+(((kappa_1c*(sqrt(T/Tcc)+1.0))/Tcc+(kappa_1c*1.0/sqrt(T/Tcc)*(T/Tcc-7.0/1.0E1)*(1.0/2.0))/Tcc)*1.0/sqrt(T/Tcc))/Tcc)*(3.294725400993198E16/2.802896292887321E15))/T+1.0/(T*T)*Tcc*(((kappa_1c*(sqrt(T/Tcc)+1.0))/Tcc+(kappa_1c*1.0/sqrt(T/Tcc)*(T/Tcc-7.0/1.0E1)*(1.0/2.0))/Tcc)*(sqrt(T/Tcc)-1.0)-(1.0/sqrt(T/Tcc)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcc)*((sqrt(T/Tcc)-1.0)*(omega_c*(6.709081269968127E15/4.503599627370496E15)-(omega_c*omega_c)*(3.086199370758719E15/1.801439850948198E16)+(omega_c*omega_c*omega_c)*(5.665283335412355E15/2.882303761517117E17)-kappa_1c*(sqrt(T/Tcc)+1.0)*(T/Tcc-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*(6.589450801986396E16/2.802896292887321E15);

                epsilon_cj(0, 0, 0, d2epsilon_c_dT2 );

                }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}


//
// epsilon_w
//

int  MolarDensity::epsilon_w_jet(const double T, int degree, JetMatrix &epsilon_wj){

    if (degree >= 0){
         
        double epsilon_w = (Tcw*pow((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0,2.0)*(1.647362700496599E16/2.802896292887321E15))/T;
        epsilon_wj(0,epsilon_w) ;
 
        if (degree >= 1){

            double depsilon_w_dT = 1.0/(T*T)*Tcw*pow((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0,2.0)*(-1.647362700496599E16/2.802896292887321E15)-(Tcw*(((kappa_1w*(sqrt(T/Tcw)+1.0))/Tcw+(kappa_1w*1.0/sqrt(T/Tcw)*(T/Tcw-7.0/1.0E1)*(1.0/2.0))/Tcw)*(sqrt(T/Tcw)-1.0)-(1.0/sqrt(T/Tcw)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcw)*((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*(3.294725400993198E16/2.802896292887321E15))/T;
 
            epsilon_wj(0, 0, depsilon_w_dT);

            if (degree == 2){
                                
                double d2epsilon_w_dT2 = (Tcw*pow(((kappa_1w*(sqrt(T/Tcw)+1.0))/Tcw+(kappa_1w*1.0/sqrt(T/Tcw)*(T/Tcw-7.0/1.0E1)*(1.0/2.0))/Tcw)*(sqrt(T/Tcw)-1.0)-(1.0/sqrt(T/Tcw)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcw,2.0)*(3.294725400993198E16/2.802896292887321E15))/T+1.0/(T*T*T)*Tcw*pow((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0,2.0)*(3.294725400993198E16/2.802896292887321E15)-(Tcw*((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*((sqrt(T/Tcw)-1.0)*(1.0/(Tcw*Tcw)*kappa_1w*1.0/sqrt(T/Tcw)-1.0/(Tcw*Tcw)*kappa_1w*1.0/pow(T/Tcw,3.0/2.0)*(T/Tcw-7.0/1.0E1)*(1.0/4.0))+1.0/(Tcw*Tcw)*1.0/pow(T/Tcw,3.0/2.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/4.0)+(((kappa_1w*(sqrt(T/Tcw)+1.0))/Tcw+(kappa_1w*1.0/sqrt(T/Tcw)*(T/Tcw-7.0/1.0E1)*(1.0/2.0))/Tcw)*1.0/sqrt(T/Tcw))/Tcw)*(3.294725400993198E16/2.802896292887321E15))/T+1.0/(T*T)*Tcw*(((kappa_1w*(sqrt(T/Tcw)+1.0))/Tcw+(kappa_1w*1.0/sqrt(T/Tcw)*(T/Tcw-7.0/1.0E1)*(1.0/2.0))/Tcw)*(sqrt(T/Tcw)-1.0)-(1.0/sqrt(T/Tcw)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)*(1.0/2.0))/Tcw)*((sqrt(T/Tcw)-1.0)*(omega_w*(6.709081269968127E15/4.503599627370496E15)-(omega_w*omega_w)*(3.086199370758719E15/1.801439850948198E16)+(omega_w*omega_w*omega_w)*(5.665283335412355E15/2.882303761517117E17)-kappa_1w*(sqrt(T/Tcw)+1.0)*(T/Tcw-7.0/1.0E1)+6.825529494453157E15/1.801439850948198E16)-1.0)*(6.589450801986396E16/2.802896292887321E15);

                epsilon_wj(0, 0, 0, d2epsilon_w_dT2 );

                }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}*/



//
// tau12
//

int  MolarDensity::tau12_jet(const double T, int degree, JetMatrix &tau12j){

    if (degree >= 0){
         
        double tau12 = Aux_G12/T + Delta_G1_12/R;
        tau12j.set(0,tau12); 

        if (degree >= 1){

            double dtau12_dT = -Aux_G12/(T*T);
            tau12j.set(0, 0, dtau12_dT);

            if (degree == 2){
                                
                double d2tau12_dT2 = 2.0*Aux_G12/(T*T*T) ;
                tau12j.set(0, 0, 0, d2tau12_dT2 );
                }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}

//
// tau21
//

int  MolarDensity::tau21_jet(const double T, int degree, JetMatrix &tau21j){

    if (degree >= 0){
         
        double tau21 = Aux_G21/T + Delta_G1_21/R;
        tau21j.set(0,tau21); 

        if (degree >= 1){

            double dtau21_dT = -Aux_G21/(T*T);
            tau21j.set(0, 0, dtau21_dT);

            if (degree == 2){
                                
                double d2tau21_dT2 = 2.0*Aux_G21/(T*T*T);
                tau21j.set(0, 0, 0, d2tau21_dT2 );

                }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}


//
// G12
//

int  MolarDensity::G12_jet(const double T, int degree, JetMatrix &G12j){

    if (degree >= 0){
         
        double G12 = exp(alpha12*(Aux_G12/T+Delta_G1_12/R));
        G12j.set(0,G12); 

        if (degree >= 1){

            double dG12_dT = -alpha12*Aux_G12*exp(alpha12*(Aux_G12/T+Delta_G1_12/R))/(T*T);
            G12j.set(0, 0, dG12_dT);

            if (degree == 2){
                                
                double d2G12_dT2 = alpha12*Aux_G12*(alpha12*Aux_G12/T + 2.0)*exp(alpha12*(Aux_G12/T+Delta_G1_12/R))/(T*T*T);

                G12j.set(0, 0, 0, d2G12_dT2 );

                }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}


//
// G21
//

int  MolarDensity::G21_jet(const double T, int degree, JetMatrix &G21j){

    if (degree >= 0){
         
        double G21 = exp(alpha12*(Aux_G21/T+Delta_G1_21/R));
        G21j.set(0,G21); 

        if (degree >= 1){

            double dG21_dT = -alpha12*Aux_G21*exp(alpha12*(Aux_G21/T+Delta_G1_21/R))/(T*T);
            G21j.set(0, 0, dG21_dT);

            if (degree == 2){
                                
                double d2G21_dT2 = alpha12*Aux_G21*(alpha12*Aux_G21/T + 2.0)*exp(alpha12*(Aux_G21/T+Delta_G1_21/R))/(T*T*T);

                G21j.set(0, 0, 0, d2G21_dT2 );

                }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}

//
// Below, the JETS of Gamma1 and Gamma2.  They change depending on the phase type.
// To use these jets we must load the jets of tau12, tau21, G12 and G21.
//
									
/** TODO: Notation conflict. BE CAREFUL: ************************
 *								*
 * For the NRTL activity coefficient it is necesary to evaluate	*
 *   G^E                 tau21 G21           tau12 G12		*
 * ------- = x1 x2 ( ----------------- + ----------------- )	*
 *   R T                x1 + x2 G21         x2 + x1 G12		*
 * The two fractions on the RHS are need and we call them here,	*
 * respectively, [Gamma1] and [Gamma2].				*
 *    The actual [gamma1] and [gamma2], which are exponential	*
 * expressions based on these [Gamma#] are computed elsewhere	*
 * in this file to be used in the file VLE_Flash_TPCW.cpp	*
 *								*
*****************************************************************/

//
// Gamma1
//

int MolarDensity::Gamma1_jet(const double x, const double T, int degree, JetMatrix &Gamma1j){

    JetMatrix tau21j(1);
    tau21_jet(T, degree, tau21j);

    JetMatrix  G21j(1);
    G21_jet(T, degree, G21j);

    if (degree >= 0){

        double frac, SgnDiff, SgnPar; 
      
        double tau21 = tau21j.get(0);
        double G21   =   G21j.get(0);

        //std::cout << "    G21 = " << G21 << ", tau21 = " << tau21 << std::endl;

        if (type == MOLAR_DENSITY_VAPOR) {
            frac    = 1.0/( (1.0 - x) + x*G21 );
            SgnDiff = -x;                        // This is the "sign" of the differential dfrac_dT
            SgnPar  = -1.0;                      // This is the "sign" inside d2Gamma1_dxdT
        }
        else { /* type == MOLAR_DENSITY_LIQUID */
            frac    = 1.0/( x + (1.0 - x)*G21 );
            SgnDiff = x - 1.0;                   // This is the "sign" of the differential dfrac_dT
            SgnPar  = 1.0;                       // This is the "sign" inside d2Gamma1_dxdT and dfrac_dx
        }

        double Gamma1 = tau21*G21*frac;
        Gamma1j.set(0,Gamma1);

        //std::cout << "    Gamma1 = " << Gamma1 << std::endl;
        
        if (degree >= 1){
            double dG21_dT   =   G21j.get(0, 0); //std::cout << "    dG21_dT = " << dG21_dT << std::endl;
            double dtau21_dT = tau21j.get(0, 0); //std::cout << "    dtau21_dT = " << dtau21_dT << std::endl;


            double frac2    = frac*frac;
            double dfrac_dT = SgnDiff*dG21_dT*frac2;
            double dfrac_dx = - SgnPar*(1.0 - G21)*frac2;

            double dGamma1_dx = tau21*G21*dfrac_dx;
            double dGamma1_dT = ( dtau21_dT*G21 + tau21*dG21_dT )*frac + tau21*G21*dfrac_dT;

            Gamma1j.set(0, 0, dGamma1_dx);
            Gamma1j.set(0, 1, dGamma1_dT);

            if (degree == 2){
                                
                double d2tau21_dT2 = tau21j.get(0,0,0);
                double d2G21_dT2   =   G21j.get(0,0,0);

                double frac3 = frac*frac2;

                double d2Gamma1_dx2  = 2.0*tau21*G21*( 1.0 - G21 )*( 1.0 - G21 )*frac3;
                double d2Gamma1_dxdT = ( dtau21_dT*G21 + tau21*dG21_dT )*dfrac_dx + SgnPar*tau21*G21*dG21_dT*( frac2 - 2.0*( 1.0 - G21 )*SgnDiff*frac3 );
                double d2Gamma1_dTdx =  d2Gamma1_dxdT;
                double d2Gamma1_dT2  = ( d2tau21_dT2*G21 + 2.0*dtau21_dT*dG21_dT + tau21*d2G21_dT2 )*frac + 2.0*( dtau21_dT*G21 + tau21*dG21_dT )*dfrac_dT + tau21*G21*SgnDiff*( d2G21_dT2 + 2.0*SgnDiff*dG21_dT*dG21_dT*frac )*frac2;

                Gamma1j.set(0, 0, 0, d2Gamma1_dx2 );
                Gamma1j.set(0, 0, 1, d2Gamma1_dxdT);
                Gamma1j.set(0, 1, 0, d2Gamma1_dTdx);
                Gamma1j.set(0, 1, 1, d2Gamma1_dT2 );
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}


//
// Gamma2
//

int MolarDensity::Gamma2_jet(const double x, const double T, int degree, JetMatrix &Gamma2j){

    JetMatrix tau12j(1);
    tau12_jet(T, degree, tau12j);

    JetMatrix  G12j(1);
    G12_jet(T, degree, G12j);

    if (degree >= 0){

        double frac, SgnDiff, SgnPar; 
      
        double tau12 = tau12j.get(0);
        double G12   =   G12j.get(0);

        if (type == MOLAR_DENSITY_VAPOR) {
            frac    = 1.0/( x + (1.0 - x)*G12 );
            SgnDiff = x - 1.0;                   // This is the "sign" of the differential dfrac_dT
            SgnPar  = 1.0;                       // This is the "sign" inside d2Gamma2_dxdT and dfrac_dx
        }
        else { /* type == MOLAR_DENSITY_LIQUID */
            frac    = 1.0/( (1.0 - x) + x*G12 );
            SgnDiff = -x;                        // This is the "sign" of the differential dfrac_dT
            SgnPar  = -1.0;                      // This is the "sign" inside d2Gamma2_dxdT
        }

        double Gamma2 = tau12*G12*frac;
        Gamma2j.set(0,Gamma2);
        
        if (degree >= 1){
            double dtau12_dT = tau12j.get(0,0);
            double dG12_dT   =   G12j.get(0,0);

            double frac2    = frac*frac;
            double dfrac_dT = SgnDiff*dG12_dT*frac2;
            double dfrac_dx = - SgnPar*(1.0 - G12)*frac2;

            double dGamma2_dx = tau12*G12*dfrac_dx;
            double dGamma2_dT = ( dtau12_dT*G12 + tau12*dG12_dT )*frac + tau12*G12*dfrac_dT;

            Gamma2j.set(0, 0, dGamma2_dx);
            Gamma2j.set(0, 1, dGamma2_dT);

            if (degree == 2){
                                
                double d2tau12_dT2 = tau12j.get(0,0,0);
                double d2G12_dT2   =   G12j.get(0,0,0);

                double frac3 = frac*frac2;

                double d2Gamma2_dx2  = 2.0*tau12*G12*( 1.0 - G12 )*( 1.0 - G12 )*frac3;
                double d2Gamma2_dxdT = ( dtau12_dT*G12 + tau12*dG12_dT )*dfrac_dx + SgnPar*tau12*G12*dG12_dT*( frac2 - 2.0*( 1.0 - G12 )*SgnDiff*frac3 ); 
                double d2Gamma2_dTdx =  d2Gamma2_dxdT;
                double d2Gamma2_dT2  = ( d2tau12_dT2*G12 + 2.0*dtau12_dT*dG12_dT + tau12*d2G12_dT2 )*frac + 2.0*( dtau12_dT*G12 + tau12*dG12_dT )*dfrac_dT + tau12*G12*SgnDiff*( d2G12_dT2 + 2.0*SgnDiff*dG12_dT*dG12_dT*frac )*frac2;

                Gamma2j.set(0, 0, 0, d2Gamma2_dx2 );
                Gamma2j.set(0, 0, 1, d2Gamma2_dxdT);
                Gamma2j.set(0, 1, 0, d2Gamma2_dTdx);
                Gamma2j.set(0, 1, 1, d2Gamma2_dT2 );
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}


//
// gamma1
//


int MolarDensity::gamma1_jet(const double x, const double T, int degree, JetMatrix &gamma1j){

    JetMatrix G12j(1);
    G12_jet(T, degree, G12j);

    JetMatrix G21j(1);
    G21_jet(T, degree, G21j);

    JetMatrix Gamma1j(2);
    Gamma1_jet(x, T, degree, Gamma1j);

    JetMatrix Gamma2j(2);
    Gamma2_jet(x, T, degree, Gamma2j);

    if (degree >= 0){
        double x1, x2, Signal;
        if (type == MOLAR_DENSITY_VAPOR) { // All x derivatives are x1 derivatives by definition
            x1 = (1.0 - x);
            x2 = x;
            Signal = -1.0;
        }
        else { /* type == MOLAR_DENSITY_LIQUID */
            x1 = x;
            x2 = (1.0 - x);
            Signal = 1.0;
        }

        double G12 = G12j.get(0);  
        double G21 = G21j.get(0);  

        double Gamma1 = Gamma1j.get(0);  
        double Gamma2 = Gamma2j.get(0);

        double x22  = x2*x2;
        double div1 = 1.0/(x1 + x2*G21);
        double div2 = 1.0/(x2 + x1*G12);
        double interior = G21*Gamma1*div1 + Gamma2*div2;

        double gamma1 = exp( x22*interior );
        gamma1j.set(0,gamma1) ;
        
        if (degree >= 1){ //TODO: Ainda nao foram testados os graus maiores do que 0 para este JET
            double dG12_dT = G12j.get(0,0) ;        
            double dG21_dT = G21j.get(0,0) ;        

            double dGamma1_dx = Gamma1j.get(0,0);
            double dGamma1_dT = Gamma1j.get(0,1);
            double dGamma2_dx = Gamma2j.get(0,0);
            double dGamma2_dT = Gamma2j.get(0,1);

            double div1_2 = div1 * div1;
            double div2_2 = div2 * div2;

            double ddiv1_dx = -Signal*div1_2*(1.0 - G21);
            double ddiv1_dT = -div1_2*x2*dG21_dT;
            double ddiv2_dx =  Signal*div2_2*(1.0 - G12);
            double ddiv2_dT = -div2_2*x1*dG12_dT;

            double dinterior_dx = G21*(dGamma1_dx*div1 + Gamma1*ddiv1_dx) + dGamma2_dx*div2 + Gamma2*ddiv2_dx;
            double dinterior_dT = dG21_dT*Gamma1*div1 + G21*(dGamma1_dT*div1 + Gamma1*ddiv1_dT) + dGamma2_dT*div2 + Gamma2*ddiv2_dT;
            double aux_interior = x22*dinterior_dx - Signal*2.0*x2*interior;

            double dgamma1_dx = aux_interior*gamma1;
            double dgamma1_dT = x22*dinterior_dT*gamma1;

            gamma1j.set(0, 0, dgamma1_dx);
            gamma1j.set(0, 1, dgamma1_dT);

            if (degree == 2){
                double d2G12_dT2 = G12j.get(0,0,0);
                double d2G21_dT2 = G21j.get(0,0,0);

                double d2Gamma1_dx2  = Gamma1j.get(0,0,0);
                double d2Gamma1_dxdT = Gamma1j.get(0,0,1);
                double d2Gamma1_dTdx = Gamma1j.get(0,1,0);
                double d2Gamma1_dT2  = Gamma1j.get(0,1,1); 
                double d2Gamma2_dx2  = Gamma2j.get(0,0,0);
                double d2Gamma2_dxdT = Gamma2j.get(0,0,1);
                double d2Gamma2_dTdx = Gamma2j.get(0,1,0); 
                double d2Gamma2_dT2  = Gamma2j.get(0,1,1);

                double div1_3 = div1_2 * div1;
                double div2_3 = div2_2 * div2;

                double d2div1_dx2  = 2.0*div1_3*(1.0 - G21)*(1.0 - G21);
                double d2div1_dxdT = Signal*(2.0*div1*(1.0 - G21)*x2 + 1.0)*div1_2*dG21_dT;
                double d2div1_dT2  = (2.0*div1*x2*dG21_dT*dG21_dT - d2G21_dT2)*x2*div1_2;
                double d2div2_dx2  = 2.0*div2_3*(1.0 - G12)*(1.0 - G12);
                double d2div2_dxdT = -Signal*(2.0*div2*(1.0 - G12)*x1 + 1.0)*div2_2*dG12_dT;
                double d2div2_dT2  = (2.0*div2*x1*dG12_dT*dG12_dT - d2G12_dT2)*x1*div2_2;

                double d2interior_dx2  = G21*(d2Gamma1_dx2*div1 + 2.0*dGamma1_dx*ddiv1_dx + Gamma1*d2div1_dx2)
                                       + d2Gamma2_dx2*div2 + 2.0*dGamma2_dx*ddiv2_dx + Gamma2*d2div2_dx2;
                double d2interior_dxdT = dG21_dT*(dGamma1_dx*div1 + Gamma1*ddiv1_dx)
                                       + G21*(d2Gamma1_dxdT*div1 + dGamma1_dT*ddiv1_dx + dGamma1_dx*ddiv1_dT + Gamma1*d2div1_dxdT)
                                       + d2Gamma2_dxdT*div2 + dGamma2_dT*ddiv2_dx + dGamma2_dx*ddiv2_dT + Gamma2*d2div2_dxdT;
                double d2interior_dT2  = d2G21_dT2*Gamma1*div1 + 2.0*dG21_dT*(dGamma1_dT*div1 + Gamma1*ddiv1_dT)
                                       + G21*(d2Gamma1_dT2*div1 + 2.0*dGamma1_dT*ddiv1_dT + Gamma1*d2div1_dT2)
                                       + d2Gamma2_dT2*div2 + 2.0*dGamma2_dT*ddiv2_dT + Gamma2*d2div2_dT2;

                double d2gamma1_dx2  = (aux_interior*aux_interior + 2.0*interior - 4.0*Signal*x2*dinterior_dx + x22*d2interior_dx2)*gamma1;
                double d2gamma1_dxdT = (x22*d2interior_dxdT - Signal*2.0*x2*dinterior_dT)*gamma1 + x22*dinterior_dT*dgamma1_dx;
                double d2gamma1_dTdx = d2gamma1_dxdT;
                double d2gamma1_dT2  = x22*(d2interior_dT2 + x22*dinterior_dT*dinterior_dT)*gamma1;

                gamma1j.set(0, 0, 0, d2gamma1_dx2 );
                gamma1j.set(0, 0, 1, d2gamma1_dxdT);
                gamma1j.set(0, 1, 0, d2gamma1_dTdx);
                gamma1j.set(0, 1, 1, d2gamma1_dT2 );
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}


//
// gamma2
//

int MolarDensity::gamma2_jet(const double x, const double T, int degree, JetMatrix &gamma2j){

    JetMatrix G12j(1);
    G12_jet(T, degree, G12j);

    JetMatrix G21j(1);
    G21_jet(T, degree, G21j);

    JetMatrix Gamma1j(2);
    Gamma1_jet(x, T, degree, Gamma1j);

    JetMatrix Gamma2j(2);
    Gamma2_jet(x, T, degree, Gamma2j);

    if (degree >= 0){
        double x1, x2, Signal;
        if (type == MOLAR_DENSITY_VAPOR) { // All x derivatives are x1 derivatives by definition
            x1 = (1.0 - x);
            x2 = x;
            Signal = -1.0;
        }
        else { /* type == MOLAR_DENSITY_LIQUID */
            x1 = x;
            x2 = (1.0 - x);
            Signal = 1.0;
        }

        double G12 = G12j.get(0);  
        double G21 = G21j.get(0);  

        double Gamma1 = Gamma1j.get(0);  
        double Gamma2 = Gamma2j.get(0);

        double x12  = x1*x1;
        double div1 = 1.0/(x2 + x1*G12);
        double div2 = 1.0/(x1 + x2*G21);
        double interior = G12*Gamma2*div1 + Gamma1*div2;

        double gamma2 = exp( x12*interior );
        gamma2j.set(0,gamma2) ;
        
        if (degree >= 1){ //TODO: Ainda nao foram testados os graus maiores do que 0 para este JET
            double dG12_dT = G12j.get(0,0) ;        
            double dG21_dT = G21j.get(0,0) ;        

            double dGamma1_dx = Gamma1j.get(0,0);
            double dGamma1_dT = Gamma1j.get(0,1);
            double dGamma2_dx = Gamma2j.get(0,0);
            double dGamma2_dT = Gamma2j.get(0,1);

            double div1_2 = div1 * div1;
            double div2_2 = div2 * div2;

            double ddiv1_dx =  Signal*div1_2*(1.0 - G12);
            double ddiv1_dT = -div1_2*x1*dG12_dT;
            double ddiv2_dx = -Signal*div2_2*(1.0 - G21);
            double ddiv2_dT = -div2_2*x2*dG21_dT;

            double dinterior_dx = G12*(dGamma2_dx*div1 + Gamma2*ddiv1_dx) + dGamma1_dx*div2 + Gamma1*ddiv2_dx;
            double dinterior_dT = dG12_dT*Gamma2*div1 + G12*(dGamma2_dT*div1 + Gamma2*ddiv1_dT) + dGamma1_dT*div2 + Gamma1*ddiv2_dT;
            double aux_interior = x12*dinterior_dx + Signal*2.0*x1*interior;

            double dgamma2_dx = aux_interior*gamma2;
            double dgamma2_dT = x12*dinterior_dT*gamma2;

            gamma2j.set(0, 0, dgamma2_dx);
            gamma2j.set(0, 1, dgamma2_dT);

            if (degree == 2){
                double d2G12_dT2 = G12j.get(0,0,0);
                double d2G21_dT2 = G21j.get(0,0,0);

                double d2Gamma1_dx2  = Gamma1j.get(0,0,0);
                double d2Gamma1_dxdT = Gamma1j.get(0,0,1);
                double d2Gamma1_dTdx = Gamma1j.get(0,1,0);
                double d2Gamma1_dT2  = Gamma1j.get(0,1,1); 
                double d2Gamma2_dx2  = Gamma2j.get(0,0,0);
                double d2Gamma2_dxdT = Gamma2j.get(0,0,1);
                double d2Gamma2_dTdx = Gamma2j.get(0,1,0); 
                double d2Gamma2_dT2  = Gamma2j.get(0,1,1);

                double div1_3 = div1_2 * div1;
                double div2_3 = div2_2 * div2;

                double d2div1_dx2  = 2.0*div1_3*(1.0 - G12)*(1.0 - G12);
                double d2div1_dxdT = -Signal*(2.0*div1*(1.0 - G12)*x1 + 1.0)*div1_2*dG12_dT;
                double d2div1_dT2  = (2.0*div1*x1*dG12_dT*dG12_dT - d2G12_dT2)*x1*div1_2;
                double d2div2_dx2  = 2.0*div2_3*(1.0 - G21)*(1.0 - G21);
                double d2div2_dxdT =  Signal*(2.0*div2*(1.0 - G21)*x2 + 1.0)*div2_2*dG21_dT;
                double d2div2_dT2  = (2.0*div2*x2*dG21_dT*dG21_dT - d2G21_dT2)*x2*div2_2;
      
                double d2interior_dx2  = G12*(d2Gamma2_dx2*div1 + 2.0*dGamma2_dx*ddiv1_dx + Gamma2*d2div1_dx2)
                                       + d2Gamma1_dx2*div2 + 2.0*dGamma1_dx*ddiv2_dx + Gamma1*d2div2_dx2;
                double d2interior_dxdT = dG12_dT*(dGamma2_dx*div1 + Gamma2*ddiv1_dx)
                                       + G12*(d2Gamma2_dxdT*div1 + dGamma2_dT*ddiv1_dx + dGamma2_dx*ddiv1_dT + Gamma2*d2div1_dxdT)
                                       + d2Gamma1_dxdT*div2 + dGamma1_dT*ddiv2_dx + dGamma1_dx*ddiv2_dT + Gamma1*d2div2_dxdT;
                double d2interior_dT2  = d2G12_dT2*Gamma2*div1 + 2.0*dG12_dT*(dGamma2_dT*div1 + Gamma2*ddiv1_dT)
                                       + G12*(d2Gamma2_dT2*div1 + 2.0*dGamma2_dT*ddiv1_dT + Gamma2*d2div1_dT2) + d2Gamma1_dT2*div2
                                       + 2.0*dGamma1_dT*ddiv2_dT + Gamma1*d2div2_dT2;

                double d2gamma2_dx2  = (aux_interior*aux_interior + 2.0*interior + Signal*4.0*x1*dinterior_dx + x12*d2interior_dx2)*gamma2;
                double d2gamma2_dxdT = (Signal*2.0*x1*dinterior_dT + x12*d2interior_dxdT)*gamma2 + x12*dinterior_dT*dgamma2_dx;
                double d2gamma2_dTdx = d2gamma2_dxdT;
                double d2gamma2_dT2  = x12*(d2interior_dT2 + x12*dinterior_dT*dinterior_dT)*gamma2;

                gamma2j.set(0, 0, 0, d2gamma2_dx2 );
                gamma2j.set(0, 0, 1, d2gamma2_dxdT);
                gamma2j.set(0, 1, 0, d2gamma2_dTdx);
                gamma2j.set(0, 1, 1, d2gamma2_dT2 );
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}



// JET GE.  It corresponds to Gibbs Free Energy, with linearization of the delta_G's.  Optimization
// procedure performed by Ali Akbar Eftekhari at TUDelft, The Netherlands. a.a.eftekhari@tudelft.nl
//
// Expressions may differ for vapor and liquid.
//

int  MolarDensity::GE_jet(const double x, const double T, int degree, JetMatrix &gej){

    // Here we load the sub-jets of GE
   
    JetMatrix Gamma1j(2);
    Gamma1_jet(x, T, degree, Gamma1j);  
 
    JetMatrix Gamma2j(2);
    Gamma2_jet(x, T, degree, Gamma2j);   

 
    if (degree >= 0){
        double Gamma1 = Gamma1j.get(0);
        double Gamma2 = Gamma2j.get(0);  
        double x2 = x*x;
        double GE = (x-x2)*( Gamma1 + Gamma2 );

        gej.set(0,GE);
        
        if (degree >= 1){           
            double dGamma1_dx = Gamma1j.get(0,0);
            double dGamma1_dT = Gamma1j.get(0,1);

            double dGamma2_dx = Gamma2j.get(0,0);
            double dGamma2_dT = Gamma2j.get(0,1);

            double dGE_dx = (1.0-2.0*x)*( Gamma1 + Gamma2 ) + (x-x2)*( dGamma1_dx + dGamma2_dx ); 
            double dGE_dT = (x- x2)*( dGamma1_dT + dGamma2_dT );
                   
            gej.set(0, 0, dGE_dx);
            gej.set(0, 1, dGE_dT);

            if (degree == 2){
                 
                double d2Gamma1_dx2  = Gamma1j.get(0,0,0);
                double d2Gamma1_dxdT = Gamma1j.get(0,0,1);
                double d2Gamma1_dTdx = Gamma1j.get(0,1,0);
                double d2Gamma1_dT2  = Gamma1j.get(0,1,1);
 
                double d2Gamma2_dx2  = Gamma2j.get(0,0,0);
                double d2Gamma2_dxdT = Gamma2j.get(0,0,1);
                double d2Gamma2_dTdx = Gamma2j.get(0,1,0); 
                double d2Gamma2_dT2  = Gamma2j.get(0,1,1);

                double d2GE_dx2 =  -2.0*(Gamma1 + Gamma2) + 2.0*(1.0-2.0*x)*( dGamma1_dx + dGamma2_dx ) + (x-x2)*(d2Gamma1_dx2 + d2Gamma2_dx2 );
                double d2GE_dxdT = (1.0-2.0*x)*( dGamma1_dT + dGamma2_dT ) + (x-x2)*( d2Gamma1_dxdT + d2Gamma2_dxdT );
                double d2GE_dTdx = d2GE_dxdT ;
                double d2GE_dT2  = (x-x2)*( d2Gamma1_dT2 + d2Gamma2_dT2 ); 

                gej.set(0, 0, 0, d2GE_dx2);
                gej.set(0, 0, 1, d2GE_dxdT);
                gej.set(0, 1, 0, d2GE_dTdx);
                gej.set(0, 1, 1, d2GE_dT2);
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}


// b_jet  This calculates the mixing coefficient b depending on the physical situation.

int  MolarDensity::vapor_b_jet(const double x, int degree, JetMatrix &bj){
    ////printf("    MolarDensity::vapor_b_jet\n");
    if (degree >= 0){
        double b = (1.0 - x)*bc + x*bw; 
        bj.set(0, b);

        if (degree >= 1){
            double db_dx = bw - bc;
            bj.set(0, 0, db_dx);

            if (degree == 2){
                double d2b_dx2 = 0.0;
                bj.set(0, 0, 0, d2b_dx2);
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}

int  MolarDensity::liquid_b_jet(const double x, int degree, JetMatrix &bj){
    ////printf("    MolarDensity::liquid_b_jet 1\n");
    if (degree >= 0){
        double b = x*bc + (1.0 - x)*bw;     ////printf("    MolarDensity::liquid_b_jet 2. Size = %d\n", bj.size());
        bj.set(0, b);    ////printf("    MolarDensity::liquid_b_jet 3\n");

        if (degree >= 1){
            double db_dx = bc - bw;
            bj.set(0, 0, db_dx);

            if (degree == 2){
                double d2b_dx2 = 0.0;
                bj.set(0, 0, 0, d2b_dx2);
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}

int  MolarDensity::b_jet(const double x, int degree, JetMatrix &bj){
    int info; ////printf("    MolarDensity::b_jet\n");
    if (type == MOLAR_DENSITY_VAPOR) info = vapor_b_jet(x, degree, bj);
    else                             info = liquid_b_jet(x, degree, bj);

    return info;
}




//
//
// The L_jet calculates the JETS for the L term.  Notice that its expression may change
// depending if we are in a LIQUID or VAPOR situation.
// 
// L = -q1[ (1-yw)*ec(T) + yw*ew(T) ] -q2[ (1-yw)*ec(T)^2 + yw*ew(T)^2 ] - GE(yw,T) + [ (1-yw)*ln bc + yw*ln bw ]     - ln b(yw)   . VAPOR
//
// L = -q1[ (xc)*ec(T) + (1-xc)*ew(T) ] -q2[ (xc)*ec(T)^2 + (1-xc)*ew(T)^2 ] - GE(yw,T) + [ xc*ln bc + (1-xc)*ln bw ] - ln b(xc)   . LIQUID
// 
//

int  MolarDensity::vapor_L_jet(const double x, const double T, int degree, JetMatrix &Lj){

    JetMatrix gej(2);
    GE_jet(x,T, degree, gej);    

    JetMatrix  bj(1);
    b_jet(x,degree,bj);

    JetMatrix  ecj(1);
    epsilon_c_jet(T,degree,ecj);

    JetMatrix  ewj(1);
    epsilon_w_jet(T,degree,ewj);
     

    if (degree >= 0){

        double GE    = gej.get(0);
        double b     = bj.get(0);
        double log_b = log(b);
        double ec    = ecj.get(0);
        double ew    = ewj.get(0);        
        double L     = -q1*( (1.0-x)*ec + x*ew ) - q2*( (1.0-x)*ec*ec + x*ew*ew ) - GE - ( (1.0-x)*log(b/bc) + x*log(b/bw) ) ;

        //cout << "L (vapor) " << L << endl;

        Lj.set(0,L);
        
        if (degree >= 1){
            
            double dGE_dx = gej.get(0,0); 
            double dGE_dT = gej.get(0,1);
            double  db_dx = bj.get(0,0);
            double dec_dT = ecj.get(0,0); 
            double dew_dT = ewj.get(0,0);
            
            double dL_dx  = - q1*(ew - ec) - q2*( ew*ew-ec*ec )-dGE_dx+(log_bw-log_bc)-(db_dx)/b ; 
            double dL_dT  = -q1*( (1.0-x)*dec_dT+x*dew_dT ) - q2*( 2.0*(1.0-x)*ec*dec_dT + 2.0*x*ew*dew_dT ) - dGE_dT ;
                          
            Lj.set(0, 0, dL_dx);
            Lj.set(0, 1, dL_dT);

            if (degree == 2){
                double d2GE_dx2  = gej.get(0,0,0) ;
                double d2GE_dxdT = gej.get(0,0,1) ;
                double d2GE_dTdx = gej.get(0,1,0) ;
                double d2GE_dT2  = gej.get(0,1,1) ;

                double d2b_dx2   = bj.get(0, 0, 0) ;

                double d2ec_dT2  = ecj.get(0, 0, 0) ;
                double d2ew_dT2  = ewj.get(0, 0, 0) ;
     
                double d2L_dx2  = -d2GE_dx2 + ( 1.0/(b*b) )*db_dx*db_dx - ( 1.0/b )*d2b_dx2 ;
                double d2L_dxdT = -q1*( dew_dT - dec_dT ) - q2*( 2.0*ew*dew_dT - 2.0*ec*dec_dT ) - d2GE_dxdT ;
                double d2L_dTdx =  d2L_dxdT;
                double d2L_dT2  = -q1*( (1.0-x)*d2ec_dT2 + x*d2ew_dT2 ) - 2.0*q2*( (1.0-x)*(dec_dT)*(dec_dT) + (1.0-x)*ec*d2ec_dT2 + x*(dew_dT)*(dew_dT) + x*ew*d2ew_dT2 ) - d2GE_dT2 ; 

                Lj.set(0, 0, 0, d2L_dx2);
                Lj.set(0, 0, 1, d2L_dxdT);
                Lj.set(0, 1, 0, d2L_dTdx);
                Lj.set(0, 1, 1, d2L_dT2);
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}


int  MolarDensity::liquid_L_jet(const double x, const double T, int degree, JetMatrix &Lj){

    JetMatrix gej(2);
    GE_jet(x,T, degree, gej);   //std::cout << "liquid_L_jet after GE_jet" << std::endl;

    JetMatrix  bj(1);
    b_jet(x,degree,bj);//std::cout << "liquid_L_jet after b_jet" << std::endl;

    JetMatrix  ewj(1);
    epsilon_w_jet(T,degree,ewj);//std::cout << "liquid_L_jet after epsilon_w_jet" << std::endl;

    JetMatrix  ecj(1);
    epsilon_c_jet(T,degree,ecj);//std::cout << "liquid_L_jet after epsilon_c_jet" << std::endl;

    if (degree >= 0){

        double GE    = gej.get(0) ;
        double b     = bj.get(0) ;
        double log_b = log(b) ;
        double ec    = ecj.get(0) ;
        double ew    = ewj.get(0) ;        
        double L     = -q1*( x*ec + (1.0-x)*ew ) - q2*( x*ec*ec + (1.0-x)*ew*ew ) - GE + ( x*log_bc + (1.0-x)*log_bw ) - log_b  ;

        Lj.set(0,L);
        
        if (degree >= 1){
            
            double dGE_dx = gej.get(0,0); 
            double dGE_dT = gej.get(0,1);
            double  db_dx = bj.get(0,0);
            double dec_dT = ecj.get(0,0); 
            double dew_dT = ewj.get(0,0);
            
            double dL_dx  = - q1*( ec - ew ) - q2*( ec*ec-ew*ew ) - dGE_dx + ( log_bc-log_bw ) - ( db_dx )/b; 
            double dL_dT  = - q1*( x*dec_dT  + (1.0-x)*dew_dT ) - 2.0*q2*( x*ec*dec_dT + ( 1.0-x )*ew*dew_dT ) - dGE_dT;
                          
            Lj.set(0, 0, dL_dx);
            Lj.set(0, 1, dL_dT);

            if (degree == 2){

                double d2GE_dx2  = gej.get(0,0,0);
                double d2GE_dxdT = gej.get(0,0,1);
                double d2GE_dTdx = gej.get(0,1,0);
                double d2GE_dT2  = gej.get(0,1,1);

                double d2b_dx2   = bj.get(0,0,0);

                double d2ec_dT2  = ecj.get(0,0,0);
                double d2ew_dT2  = ewj.get(0,0,0);
     
                double d2L_dx2  = -d2GE_dx2 + ( 1.0/(b*b) )*db_dx*db_dx - ( 1.0/b )*d2b_dx2;
                double d2L_dxdT = -q1*( dec_dT - dew_dT ) - 2.0*q2*( ec*dec_dT - ew*dew_dT ) - d2GE_dxdT;
                double d2L_dTdx = d2L_dxdT;
                double d2L_dT2  = -q1*( x*d2ec_dT2+(1.0-x)*d2ew_dT2) - 2.0*q2*( x*(dec_dT)*(dec_dT) + x*ec*d2ec_dT2 + (1.0-x)*(dew_dT)*(dew_dT)+(1.0-x)*ew*d2ew_dT2 ) - d2GE_dT2; 
 
                Lj.set(0, 0, 0, d2L_dx2 );
                Lj.set(0, 0, 1, d2L_dxdT);
                Lj.set(0, 1, 0, d2L_dTdx);
                Lj.set(0, 1, 1, d2L_dT2 );               
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}


int MolarDensity::L_jet(const double x, const double T, int degree, JetMatrix &Lj){
    int info;
    if (type == MOLAR_DENSITY_VAPOR) info =  vapor_L_jet(x,T,degree, Lj);
    else                             info = liquid_L_jet(x,T,degree, Lj);

    return info;
}


// The Quadratic Equation Q

int MolarDensity::Q_jet(const double x, const double T, const double epsilon, 
                        JetMatrix &Lj, int degree, JetMatrix &qj){

    if (degree >= 0){
        double L = Lj.get(0);
        double epsilon2 = epsilon*epsilon;

        double Q = q2*epsilon2 + q1*epsilon + L;
        qj.set(0, Q);

        if (degree >= 1){
            double dL_dx = Lj.get(0, 0);
            double dL_dT = Lj.get(0, 1);

            double dQ_dx, dQ_dT, dQ_depsilon;

            dQ_dx = dL_dx ; 
            dQ_dT = dL_dT ; 

            dQ_depsilon = 2.0*q2*epsilon + q1;

            qj.set(0, 0, dQ_dx);
            qj.set(0, 1, dQ_dT);
            qj.set(0, 2, dQ_depsilon);

            if (degree == 2){
                double d2L_dx2  = Lj.get(0, 0, 0);
                double d2L_dxdT = Lj.get(0, 0, 1);
                double d2L_dTdx = Lj.get(0, 1, 0);
                double d2L_dT2  = Lj.get(0, 1, 1);

                double d2Q_dx2,  d2Q_dxdT, d2Q_dxdepsilon;
                double d2Q_dTdx, d2Q_dT2,  d2Q_dTdepsilon;
                double d2Q_depsilondx, d2Q_depsilondT, d2Q_depsilon2;

                d2Q_dx2        = d2L_dx2  ;
                d2Q_dxdT       = d2L_dxdT ;                
                d2Q_dxdepsilon = 0.0 ;

                d2Q_dTdx       = d2Q_dxdT ;

                d2Q_dT2        = d2L_dT2  ; 
                d2Q_dTdepsilon = 0.0 ;

                d2Q_depsilondx = d2Q_dxdepsilon;

                d2Q_depsilondT = d2Q_dTdepsilon;

                d2Q_depsilon2  = 2.0*q2;

                qj.set(0, 0, 0, d2Q_dx2);
                qj.set(0, 0, 1, d2Q_dxdT);
                qj.set(0, 0, 2, d2Q_dxdepsilon);

                qj.set(0, 1, 0, d2Q_dTdx);
                qj.set(0, 1, 1, d2Q_dT2);
                qj.set(0, 1, 2, d2Q_dTdepsilon);

                qj.set(0, 2, 0, d2Q_depsilondx);
                qj.set(0, 2, 1, d2Q_depsilondT);
                qj.set(0, 2, 2, d2Q_depsilon2);

            }
            else return -1; // ABORTED_PROCEDURE
        }
    }

    return 2; // SUCCESSFUL_PROCEDURE
}


// epsilon is found implicitly as the largest root of the quadratic equation
//
// q2*epsilon^{2} + q1*epsilon + L(x,T) = 0
//
// For the vapor phase the value of L has to be chosen accordingly,  e.g. x=yw
// For the liquid phase the value of L has to be chosen accordingly, e.g. x=xc
// This JET does not distinguish the difference between liquid and vapor, because
// this choice is established when choosing L.

int MolarDensity::Epsilon_jet(const double x, const double T, int degree, JetMatrix &epsj){
    //printf("Inside MolarDensity::Epsilon_jet()\n");
    // Find the correct root of the quadratic.  For this we first fill the jet of L.
    JetMatrix Lj(2);
    L_jet(x, T, degree, Lj);
    //std::cout << "Epsilon_jet after Lj" << std::endl;

    if (degree >= 0){
        double L = Lj.get(0);

        double A = q2;
        double B = q1;
        double C = L;
        double D = B*B - 4.0*A*C ;

        if (D < 0.0) {
            //printf("D = %f\n", D);
            return -1; // ABORTED_PROCEDURE
        }

        double sqrt_D = sqrt(D);

        double epsilonminus = (  ( -B - sqrt_D  )/(2.0*A)  ) ;
        double epsilonplus  = (  ( -B + sqrt_D  )/(2.0*A)  ) ;

        double epsilon;
        if( epsilonminus >= epsilonplus ) epsilon = epsilonminus;
        else                              epsilon = epsilonplus;

        // Fill Q, where we know that Q(x, T, epsilon(x,T)) = 0
        JetMatrix qj(3);
        Q_jet(x, T, epsilon, Lj, degree, qj);

        epsj.set(0, epsilon);

        if (degree >= 1){
            double dQ_dx = qj.get(0, 0);
            double dQ_dT = qj.get(0, 1);
            double dQ_depsilon = qj.get(0, 2);           

            double inv_dQ_depsilon = 1.0/dQ_depsilon;
            double depsilon_dx = -dQ_dx*inv_dQ_depsilon;
            double depsilon_dT = -dQ_dT*inv_dQ_depsilon;

            epsj.set(0, 0, depsilon_dx);
            epsj.set(0, 1, depsilon_dT);

            if (degree == 2){
                double d2Q_dx2  = qj.get(0, 0, 0);
                double d2Q_dxdT = qj.get(0, 0, 1);
                double d2Q_dxdepsilon = qj.get(0, 0, 2);

                double d2Q_dTdx = qj.get(0, 1, 0);
                double d2Q_dT2  = qj.get(0, 1, 1);
                double d2Q_dTdepsilon = qj.get(0, 1, 2);

                double d2Q_depsilondx = qj.get(0, 2, 0);
                double d2Q_depsilondT = qj.get(0, 2, 1);
                double d2Q_depsilon2  = qj.get(0, 2, 2);
                

                double d2epsilon_dx2, d2epsilon_dxdT, d2epsilon_dTdx, d2epsilon_dT2;

                d2epsilon_dx2  = -( d2Q_dx2 + (2.0*d2Q_dxdepsilon + d2Q_depsilon2*depsilon_dx)*depsilon_dx )*inv_dQ_depsilon;
                d2epsilon_dxdT = -( d2Q_dTdx + d2Q_dxdepsilon*depsilon_dT + (d2Q_depsilondT + d2Q_depsilon2*depsilon_dT)*depsilon_dx )*inv_dQ_depsilon;
                d2epsilon_dTdx =    d2epsilon_dxdT;
                d2epsilon_dT2  = -( d2Q_dT2  + (2.0*d2Q_dTdepsilon + d2Q_depsilon2*depsilon_dT)*depsilon_dT )*inv_dQ_depsilon;

                epsj.set(0, 0, 0, d2epsilon_dx2);
                epsj.set(0, 0, 1, d2epsilon_dxdT);
                epsj.set(0, 1, 0, d2epsilon_dTdx);
                epsj.set(0, 1, 1, d2epsilon_dT2);
            }
        }
        return 2; // SUCCESSFUL_PROCEDURE;
    }
    else return -1; // ABORTED_PROCEDURE
}

int MolarDensity::A_jet(const double x, const double T, int degree, JetMatrix &Aj){
    //std::cout << "A_jet" << std::endl;

    JetMatrix bj(1); // b automatically distinguishes between vapor and liquid.
    b_jet(x, degree, bj);
    //std::cout << "A_jet after b_jet" << std::endl;

    JetMatrix epsj(2);
    Epsilon_jet(x, T, degree, epsj); // Epsilon distinguishes vapor and liquid inside of L
    //std::cout << "A_jet after Epsilon_jet" << std::endl;

    ////std::cout << "A_jet before degree" << std::endl;

    if (degree >= 0){
        double inv_R  = 1.0/R;
        double inv_T  = 1.0/T;
        double inv_T2 = inv_T*inv_T;

        double eps = epsj.get(0); // Epsilon
        double b = bj.get(0);

        double A = eps*b*P*inv_R*inv_T;

        Aj.set(0, A);

        if (degree >= 1){
            double deps_dx = epsj.get(0, 0);
            double deps_dT = epsj.get(0, 1);

            double db_dx = bj.get(0, 0);

            double dA_dx, dA_dT;

            dA_dx = P*inv_R*inv_T*(deps_dx*b + eps*db_dx);
            dA_dT = P*b*inv_R*(deps_dT*inv_T - eps*inv_T2);

            Aj.set(0, 0, dA_dx);
            Aj.set(0, 1, dA_dT);

            if (degree == 2){
                double d2eps_dx2  = epsj.get(0, 0, 0);
                double d2eps_dxdT = epsj.get(0, 0, 1);
                double d2eps_dTdx = epsj.get(0, 1, 0);
                double d2eps_dT2  = epsj.get(0, 1, 1);

                double d2b_dx2 = bj.get(0, 0, 0);

                double d2A_dx2, d2A_dxdT, d2A_dTdx, d2A_dT2;

                d2A_dx2  = P*inv_R*inv_T*(d2eps_dx2*b + 2.0*deps_dx*db_dx + eps*d2b_dx2);

                d2A_dxdT = -P*inv_R*inv_T2*(deps_dx*b + eps*db_dx) + 
                           P*inv_R*inv_T*(d2eps_dxdT*b + deps_dT*db_dx);

                d2A_dTdx = d2A_dxdT;

                d2A_dT2  = P*b*inv_R*(d2eps_dT2*inv_T - 2.0*inv_T2*deps_dT + 2.0*eps*inv_T*inv_T2);

                Aj.set(0, 0, 0, d2A_dx2);
                Aj.set(0, 0, 1, d2A_dxdT);
                Aj.set(0, 1, 0, d2A_dTdx);
                Aj.set(0, 1, 1, d2A_dT2);
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}

int MolarDensity::B_jet(const double x, const double T, int degree, JetMatrix &Bj){
    JetMatrix bj(1);
    b_jet(x, degree, bj);

    if (degree >= 0){
        double b = bj.get(0);

        double inv_R  = 1.0/R;
        double inv_T  = 1.0/T;
        double inv_T2 = inv_T*inv_T;
        double B = b*P*inv_R*inv_T;

        Bj.set(0, B);
        
        if (degree >= 1){
            double db_dx = bj.get(0, 0);

            double dB_dx, dB_dT;

            dB_dx =  P*inv_R*inv_T*db_dx;
            dB_dT = -P*inv_R*inv_T2*b;

            Bj.set(0, 0, dB_dx);
            Bj.set(0, 1, dB_dT);

            if (degree == 2){
                double d2b_dx2 = bj.get(0, 0, 0);

                double d2B_dx2, d2B_dxdT, d2B_dTdx, d2B_dT2;

                d2B_dx2  = P*inv_R*inv_T*d2b_dx2;

                d2B_dxdT = -P*inv_R*inv_T2*db_dx;

                d2B_dTdx = d2B_dxdT;

                d2B_dT2  = 2.0*P*inv_R*b*inv_T*inv_T2;

                Bj.set(0, 0, 0, d2B_dx2);
                Bj.set(0, 0, 1, d2B_dxdT);
                Bj.set(0, 1, 0, d2B_dTdx);
                Bj.set(0, 1, 1, d2B_dT2);
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}





int MolarDensity::F_jet(const double x, const double T, const double Z, 
                        JetMatrix &Aj, JetMatrix &Bj, 
                        int degree, 
                        JetMatrix &fj){
    if (degree >= 0){
        double A = Aj.get(0);
        double B = Bj.get(0);
        double B2 = B*B;

        double Z2 = Z*Z;

        double F;
        
        F = Z2*Z - (1.0 - B)*Z2 + (A - 2.0*B - 3.0*B2)*Z - B*(A - B - B2);
        fj.set(0, F);

        if (degree >= 1){
            double dA_dx = Aj.get(0, 0);
            double dA_dT = Aj.get(0, 1);

            double dB_dx = Bj.get(0, 0);
            double dB_dT = Bj.get(0, 1);

            double dF_dx, dF_dT, dF_dZ;

            dF_dx = dB_dx*Z2 + (dA_dx - 2.0*dB_dx - 6.0*B*dB_dx)*Z - 
                    dB_dx*(A - B - B2) - B*(dA_dx - dB_dx - 2.0*B*dB_dx);

            dF_dT = dB_dT*Z2 + (dA_dT - 2.0*dB_dT - 6.0*B*dB_dT)*Z - 
                    dB_dT*(A - B - B2) - B*(dA_dT - dB_dT - 2.0*B*dB_dT);

            dF_dZ = 3.0*Z2 - 2.0*(1.0 - B)*Z + (A - 2.0*B - 3.0*B2);

            fj.set(0, 0, dF_dx);
            fj.set(0, 1, dF_dT);
            fj.set(0, 2, dF_dZ);

            if (degree == 2){
                double d2A_dx2  = Aj.get(0, 0, 0);
                double d2A_dxdT = Aj.get(0, 0, 1);
                double d2A_dTdx = Aj.get(0, 1, 0);
                double d2A_dT2  = Aj.get(0, 1, 1);

                double d2B_dx2  = Bj.get(0, 0, 0);
                double d2B_dxdT = Bj.get(0, 0, 1);
                double d2B_dTdx = Bj.get(0, 1, 0);
                double d2B_dT2  = Bj.get(0, 1, 1); 

                double d2F_dx2,  d2F_dxdT, d2F_dxdZ;
                double d2F_dTdx, d2F_dT2,  d2F_dTdZ;
                double d2F_dZdx, d2F_dZdT, d2F_dZ2;

                d2F_dx2 = d2B_dx2*Z2 + (d2A_dx2 - 2.0*d2B_dx2 - 6.0*dB_dx*dB_dx - 6.0*B*d2B_dx2)*Z - 
                          d2B_dx2*(A - B - B2) - dB_dx*(dA_dx - dB_dx - 2.0*B*dB_dx) - 
                          dB_dx*(dA_dx - dB_dx - 2.0*B*dB_dx) -
                          B*(d2A_dx2 - d2B_dx2 - 2.0*dB_dx*dB_dx - 2.0*B*d2B_dx2);

                d2F_dxdT = d2B_dxdT*Z2 + 
                           (d2A_dxdT - 2.0*d2B_dxdT - 6.0*dB_dx*dB_dT - 6.0*B*d2B_dxdT)*Z -
                           d2B_dxdT*(A - B - B2) -  dB_dx*(dA_dT - dB_dT - 2.0*B*dB_dT) -
                           dB_dT*(dA_dx - dB_dx - 2.0*B*dB_dx) -
                           B*(d2A_dxdT - d2B_dxdT - 2.0*dB_dT*dB_dx - 2.0*B*d2B_dxdT);

                d2F_dxdZ = 2.0*dB_dx*Z + dA_dx - 2.0*dB_dx - 6.0*B*dB_dx;

                d2F_dTdx = d2F_dxdT;

                d2F_dT2 = d2B_dT2*Z2 + (d2A_dT2 - 2.0*d2B_dT2 - 6.0*dB_dT*dB_dT - 6.0*B*d2B_dT2)*Z - 
                          d2B_dT2*(A - B - B2) - dB_dT*(dA_dT - dB_dT - 2.0*B*dB_dT) -
                          dB_dT*(dA_dT - dB_dT - 2.0*B*dB_dT) -
                          B*(d2A_dT2 - d2B_dT2 - 2.0*dB_dT*dB_dT - 2.0*B*d2B_dT2);

                d2F_dTdZ = 2.0*Z*dB_dT + dA_dT - 2.0*dB_dT - 6.0*B*dB_dT;

                d2F_dZdx = d2F_dxdZ;

                d2F_dZdT = d2F_dTdZ;

                d2F_dZ2 = 6.0*Z - 2.0*(1.0 - B);

                fj.set(0, 0, 0, d2F_dx2);
                fj.set(0, 0, 1, d2F_dxdT);
                fj.set(0, 0, 2, d2F_dxdZ);

                fj.set(0, 1, 0, d2F_dTdx);
                fj.set(0, 1, 1, d2F_dT2);
                fj.set(0, 1, 2, d2F_dTdZ);

                fj.set(0, 2, 0, d2F_dZdx);
                fj.set(0, 2, 1, d2F_dZdT);
                fj.set(0, 2, 2, d2F_dZ2);

            }
            else return -1; // ABORTED_PROCEDURE
        }
    }

    return 2; // SUCCESSFUL_PROCEDURE
}

//
// Z is found implicitly as the root of the following polynomial, given by the
// PRSV equation of state:
//
//    Z^3 - (1 - B)*Z^2 + (A - 2*B - 3*B^2)*Z - B*(A - B - B^2) = 0.
//
// For the vapor phase the largest root is selected.
// For the liquid phase the smallest root is selected instead.
//

int MolarDensity::Z_jet(const double x, const double T, int degree, JetMatrix &zj){

    //printf("x = %g, T = %g, degree = %d\n", x, T, degree);

    JetMatrix Aj(2), Bj(2);
    A_jet(x, T, degree, Aj); //printf("After A_jet\n");
    B_jet(x, T, degree, Bj); //printf("After B_jet\n");

    if (degree >= 0){
        ////printf("degree >= 0\n");
        double A = Aj.get(0);
        double B = Bj.get(0);

        //printf("Before sqrt3\n");

        // Find the roots of the polynomial.
        std::complex<double> r[3];
        int info_cubic_roots = CubicRoots(1.0, -(1.0 - B), (A - 2.0*B - 3.0*B*B), -B*(A - B - B*B), r[0], r[1], r[2]);

        std::vector<double> vroots_real(3);
        for (int i = 0; i < 3; i++) vroots_real[i] = r[i].real();

        std::vector<double> vroots_imag(3);
        for (int i = 0; i < 3; i++) vroots_imag[i] = r[i].imag();

        double max_no = vroots_real[0];
        double min_no = vroots_real[0];

        /**
         * TODO(1): In the original routine made by Ali Akbar Eftekhari in Matlab, it
         *          is pretended to take only real roots from the cubic polynomial.
         *          Even if the roots (for the supercritical CO2 case) always have two
         *          complex conjugated roots, the Matlab routine take the real part of
         *          these complex roots as [min_no]. Here we understand that the
         *          modification is small, so we let intentionally to happen.
         * TODO(2): For testing other cases, notice that the CubicRoots function returns
         *          the value 2 when complex conjugated roots are found. The two [if] in
         *          next lines must be uncommented. Recall to check that [min_no] must
         *          be small but larger than 0.0, and [max_no] almost 1.0.
         **/

       for (int i = 0; i < 3; i++) {
           if ( (vroots_real[i] > max_no) && (vroots_real[i] > 0.0) && (fabs(vroots_imag[i]) < zero_epsilon) ) {
                max_no = vroots_real[i];
           }
           if ( (vroots_real[i] < min_no) && (vroots_real[i] > 0.0) && (fabs(vroots_imag[i]) < zero_epsilon) ) {
                min_no = vroots_real[i];
           }
        }

        double vapor_z  = max_no;
        double liquid_z = min_no;

        double Z;
        if (type == MOLAR_DENSITY_VAPOR){
            Z = max_no;
        }
        else {
            Z = min_no;
        }

        //printf("    After sqrt3\n");

        // Fill F, etc.
        JetMatrix fj(3);
        F_jet(x, T, Z, Aj, Bj, degree, fj);

        zj.set(0, Z);

        if (degree >= 1){
            double dF_dx = fj.get(0, 0);
            double dF_dT = fj.get(0, 1);
            double dF_dZ = fj.get(0, 2);

            double dZ_dx, dZ_dT;
            double inv_dF_dZ = 1.0/dF_dZ;

            dZ_dx = -dF_dx*inv_dF_dZ ;
            dZ_dT = -dF_dT*inv_dF_dZ ;

            zj.set(0, 0, dZ_dx);
            zj.set(0, 1, dZ_dT);

            if (degree == 2){
                double d2F_dx2  = fj.get(0, 0, 0);
                double d2F_dxdT = fj.get(0, 0, 1);
                double d2F_dxdZ = fj.get(0, 0, 2);

                double d2F_dTdx = fj.get(0, 1, 0);
                double d2F_dT2  = fj.get(0, 1, 1);
                double d2F_dTdZ = fj.get(0, 1, 2);

                double d2F_dZdx = fj.get(0, 2, 0);
                double d2F_dZdT = fj.get(0, 2, 1);
                double d2F_dZ2  = fj.get(0, 2, 2);
              
                double d2Z_dx2, d2Z_dxdT, d2Z_dTdx, d2Z_dT2;

                d2Z_dx2  = - ( d2F_dx2 + ( 2.0*d2F_dxdZ + d2F_dZ2*dZ_dx )*dZ_dx )*inv_dF_dZ;
                d2Z_dxdT = -( d2F_dTdx + d2F_dZdx*dZ_dT + (d2F_dZdT + d2F_dZ2*dZ_dT)*dZ_dx )*inv_dF_dZ;
                d2Z_dTdx = d2Z_dxdT;
                d2Z_dT2  = - (d2F_dT2  + ( 2.0*d2F_dTdZ + d2F_dZ2*dZ_dT)*dZ_dT )*inv_dF_dZ;

                zj.set(0, 0, 0, d2Z_dx2);
                zj.set(0, 0, 1, d2Z_dxdT);
                zj.set(0, 1, 0, d2Z_dTdx);
                zj.set(0, 1, 1, d2Z_dT2);
            }
        }
        return 2; // SUCCESSFUL_PROCEDURE;
    }
    else return -1; // ABORTED_PROCEDURE
}



int MolarDensity::vapor_rho_jet(const double x, const double T, int degree, JetMatrix &r){
    JetMatrix zj(2);
    Z_jet(x, T, degree, zj);

    if (degree >= 0){
        double Z        = zj.get(0);
        double inv_ZRT  = 1.0/(Z*R*T);
        double inv_ZRT2 = inv_ZRT*inv_ZRT;
        double inv_ZRT3 = inv_ZRT2*inv_ZRT;

        double rho = P*inv_ZRT; 
        r.set(0, rho);

        if (degree >= 1){
            double dZ_dx = zj.get(0, 0);
            double dZ_dT = zj.get(0, 1);

            double drho_dx, drho_dT;

            drho_dx = -P*inv_ZRT2*(dZ_dx*R*T);
            drho_dT = -P*inv_ZRT2*(dZ_dT*R*T + Z*R);

            r.set(0, 0, drho_dx);
            r.set(0, 1, drho_dT);

            if (degree == 2){
                double d2Z_dx2  = zj.get(0, 0, 0);
                double d2Z_dxdT = zj.get(0, 0, 1);
                double d2Z_dTdx = zj.get(0, 1, 0);
                double d2Z_dT2  = zj.get(0, 1, 1);

                double d2rho_dx2, d2rho_dxdT, d2rho_dTdx, d2rho_dT2;

                d2rho_dx2  = 2.0*P*inv_ZRT3*(dZ_dx*R*T)*(dZ_dx*R*T) - 
                             P*inv_ZRT2*(d2Z_dx2*R*T);

                d2rho_dxdT = 2.0*P*inv_ZRT3*(dZ_dx*R*T)*(dZ_dT*R*T + Z*R) -
                             P*inv_ZRT2*(d2Z_dxdT*R*T + dZ_dx*R);

                d2rho_dTdx = d2rho_dxdT;

                d2rho_dT2  = 2.0*P*inv_ZRT3*(dZ_dT*R*T + Z*R)*(dZ_dT*R*T + Z*R) -
                             P*inv_ZRT2*(d2Z_dT2*R*T + 2.0*dZ_dT*R);

                r.set(0, 0, 0, d2rho_dx2);
                r.set(0, 0, 1, d2rho_dxdT);
                r.set(0, 1, 0, d2rho_dTdx);
                r.set(0, 1, 1, d2rho_dT2);
            }
        }
        return 2; // SUCCESSFUL_PROCEDURE ;
    }
    else return -1; // ABORTED_PROCEDURE
}


int  MolarDensity::liquid_C_jet(const double x, const double T, int degree, JetMatrix &Cj){

    if (degree >= 0){
        double C = x*(C0_c + C1_c*T) + (1.0 - x)*(C0_w + C1_w*T) ;
        Cj.set(0, C);

        if (degree >= 1){
            double dC_dx = C0_c - C0_w + (C1_c - C1_w)*T;
            double dC_dT = x*C1_c + (1.0 - x)*C1_w;

            Cj.set(0, 0, dC_dx);
            Cj.set(0, 1, dC_dT);

            if (degree == 2){
                double d2C_dx2, d2C_dxdT, d2C_dTdx, d2C_dT2;

                d2C_dx2  = 0.0;
                d2C_dxdT = C1_c - C1_w;
                d2C_dTdx = d2C_dxdT;
                d2C_dT2  = 0.0;

                Cj.set(0, 0, 0, d2C_dx2);
                Cj.set(0, 0, 1, d2C_dxdT);
                Cj.set(0, 1, 0, d2C_dTdx);
                Cj.set(0, 1, 1, d2C_dT2);
            }
            else return -1; // ABORTED_PROCEDURE;
        }
    }
    return 2; // SUCCESSFUL_PROCEDURE
}

int MolarDensity::liquid_rho_jet(const double x, const double T, int degree, JetMatrix &r){
    JetMatrix Cj(2); 
    liquid_C_jet(x, T, degree, Cj); 

    JetMatrix zj(2);//printf("Before Z_jet\n");
    Z_jet(x, T, degree, zj);//printf("After Z_jet\n");

    if (degree >= 0){
        double C           = Cj.get(0);
        double Z           = zj.get(0);
        double inv_ZRT_CP  = 1.0/(Z*R*T + C*P);
        double inv_ZRT_CP2 = inv_ZRT_CP*inv_ZRT_CP;
        double inv_ZRT_CP3 = inv_ZRT_CP2*inv_ZRT_CP;

        double rho = P*inv_ZRT_CP;
        r.set(0, rho);

        if (degree >= 1){
            double dZ_dx = zj.get(0, 0);
            double dZ_dT = zj.get(0, 1);

            double dC_dx = Cj.get(0, 0);
            double dC_dT = Cj.get(0, 1);

            double drho_dx, drho_dT;

            drho_dx = -P*inv_ZRT_CP2*(dZ_dx*R*T + dC_dx*P);
            drho_dT = -P*inv_ZRT_CP2*(dZ_dT*R*T + Z*R + dC_dT*P);

            r.set(0, 0, drho_dx);
            r.set(0, 1, drho_dT);

            if (degree == 2){
                double d2Z_dx2  = zj.get(0, 0, 0);
                double d2Z_dxdT = zj.get(0, 0, 1);
                double d2Z_dTdx = zj.get(0, 1, 0);
                double d2Z_dT2  = zj.get(0, 1, 1);

                double d2C_dx2  = Cj.get(0, 0, 0);
                double d2C_dxdT = Cj.get(0, 0, 1);
                double d2C_dTdx = Cj.get(0, 1, 0);
                double d2C_dT2  = Cj.get(0, 1, 1);

                double d2rho_dx2, d2rho_dxdT, d2rho_dTdx, d2rho_dT2;

                d2rho_dx2  = 2.0*P*inv_ZRT_CP3*(dZ_dx*R*T + dC_dx*P)*(dZ_dx*R*T + dC_dx*P) -
                             P*inv_ZRT_CP2*(d2Z_dx2*R*T + d2C_dx2*P);

                d2rho_dxdT = 2.0*P*inv_ZRT_CP3*(dZ_dT*R*T + Z*R + dC_dT*P)*(dZ_dx*R*T + dC_dx*P) -
                             P*inv_ZRT_CP2*(d2Z_dxdT*R*T + dZ_dx*R + d2C_dxdT*P);

                d2rho_dTdx = d2rho_dxdT;

                d2rho_dT2  = 2.0*P*inv_ZRT_CP3*(dZ_dT*R*T + Z*R + dC_dT*P)*(dZ_dT*R*T + Z*R + dC_dT*P) - 
                             P*inv_ZRT_CP2*(d2Z_dT2*R*T + 2.0*dZ_dT*R + d2C_dT2*P);

                r.set(0, 0, 0, d2rho_dx2);
                r.set(0, 0, 1, d2rho_dxdT);
                r.set(0, 1, 0, d2rho_dTdx);
                r.set(0, 1, 1, d2rho_dT2);
            }
        }
        return 2; // SUCCESSFUL_PROCEDURE;
    }
    else return -1; // ABORTED_PROCEDURE
}

int MolarDensity::rho_jet(const double x, const double T, int degree, JetMatrix &r){
    int info;
    if (type == MOLAR_DENSITY_VAPOR) info = vapor_rho_jet(x, T, degree, r);
    else                             info = liquid_rho_jet(x, T, degree, r);

    return info;    
}

