#ifndef _plotter_h
#define _plotter_h
//#ifdef __GNUC__
//#pragma interface
//#endif

#include "curvegen.h"
#include "interpol.h"

// Added by Morante, 05/09/2012
#include "RealVector.h"
#include <vector>

class RiemannPlotter:private CurveGenerator
{

    Interpolator interpol_;
    int gap_;
//1 se preenchendo lacunas, 0 do contrario.
//1 if filling gaps, 0 otherwise
    float gap_start_,gap_end_,gap_step_,current_x_start_,last_top_,current_top_;

    int right_;
//1 se preenchendo lacunas a direita, 0 do contrario.
//1 if filling right border gaps, 0 otherwise
    float rgap_start_,rgap_end_,rgap_step_,last_right_,current_right_;

    void reset(void);
    //Inicializa as variaveis de controle
    //Resets the control variables

    void get_interception(void);
/*
Processa a interseccao de uma caracteristica com o topo da janela, procurando 
por lacunas (gaps) a preencher. So e chamada apos gerar uma caracteristica
partindo da parte inferior do retangulo.

Processes the interssection of last plotted characteristic and the window's
top border, looking for gaps to fill. It is called only after generating
a characteristic starting at the window's lower border.
*/

    int setup_next(void);
/*
Prepara o gerador de curvas para rastrear a proxima caracteristica. 
Retorna 0 se nao ha mais curvas a rastrear.

Prepares to trace the next characteristics. Returns 0 if there are
no characteristics left to generate.
*/

    virtual float slope(float,float);

    float distance_;
    //Distancia maxima entre as curvas a serem plotadas
    //Max distance between curves to be plotted
    

    void (*plot_)(int);
/*
Ponteiro para a funcao de plot. Recebe o numero de pontos. Os valores de
(x,y) foram guardados nos arrays passados ao construtor, com o stride
especificado.

plot callback. Receives the number of points. The values of (x,y) have
been stored in the arrays passed to the constructor, using the requested
stride.
*/


public:
    virtual ~RiemannPlotter();
    RiemannPlotter(float xmin,float xmax,float ymax,
                   float distance,
                   int size, float *domain,float *image,
                   void (*plot)(int),float* x,float* y,
		   int x_stride,int y_stride,
                   float diff=1e-3,float sample=0.2);
/*
O construtor acima e a interface principal para as rotinas que geram
o grafico. size,domain e image sao passados para o const-
trutor do interpolador. x,y,x_stride e y_stride sao os arrays onde 
armazenar os dados para plotagem, efetuada chamando o callback plot. 
diff e sample sao passados para o construtor do gerador de curvas.

The above constructor is the main interface for the routines that
generate the graphic. size, domain and image are passed to the 
interpolators constructor, x,y,x_stride and y_stride are the arrays
where the plot data is stored, and ploting is done by the callback
plot. diff and sample are passed to the tracer's constructor.
*/

    // I modified this method to obtain the characteristics as a vector of points
    // Morante, 05/09/2012.
    //
    // void gen_family(void);
    void gen_family(std::vector<std::vector<RealVector> > &c);

    //Gera e desenha as caracteristicas
    //Creates and draws the characteristics
};
#endif
