#ifndef _DOUBLE_CONTACT_GPU_
#define _DOUBLE_CONTACT_GPU_

#include "ThreeImplicitFunctions.h"
#include "Contour2x2_Method_GPU.h"
#include "Double_Contact_Function.h"

class Double_Contact_GPU : public Double_Contact_Function {
    private:
    protected:
        RealVector lambda_left;
        Matrix<double> flux_left, accum_left;

        int left_family, right_family;

        const FluxFunction *lff;
    	const AccumulationFunction *laa;

        const FluxFunction *rff;
        const AccumulationFunction *raa;


        // For the GPU.
        //
        Matrix<float> gpu_lambda_left;
        Matrix<std::vector<float> > gpu_F_left, gpu_G_left;

        Matrix<float> gpu_lambda_right;
        Matrix<std::vector<float> > gpu_F_right, gpu_G_right;

        void gpu_allocate_arrays_in_rectangle(Double_Contact_GPU *dc, int length_left, Matrix<float> &val0, Matrix<float> &val1, Matrix<float> &val2);
        void gpu_allocate_arrays_in_triangle(Double_Contact_GPU *dc, int length_left, Matrix<float> &val0, Matrix<float> &val1, Matrix<float> &val2);
        void gpu_fill_arrays(int min_left, int max_left, Matrix<float> &val0, Matrix<float> &val1, Matrix<float> &val2);

        void gpu_check_sign(int min_left, int max_left, const Matrix<float> &val0, const Matrix<float> &val1, const Matrix<float> &val2, 
                            std::vector<short int> &index_l, std::vector<short int> &index_r);
    public:
        Double_Contact_GPU(){
            gv_left = gv_right = 0;

            lambda_left.resize(4);        
            flux_left.resize(2, 4);
            accum_left.resize(2, 4);

            left_family = right_family = 0;

            singular = false;
        }

        ~Double_Contact_GPU(){}

        bool prepare_cell(int i, int j);

        bool function_on_cell(double *val, int ir, int jr, int kl, int kr);

        void curve(const FluxFunction *lf, const AccumulationFunction *la, GridValues *lg, int lfam,
                   const FluxFunction *rf, const AccumulationFunction *ra, GridValues *rg, int rfam,
                   std::vector<RealVector> &left_curve, std::vector<RealVector> &right_curve);

        void gpu_curve(const FluxFunction *lf, const AccumulationFunction *la, GridValues *lg, int lfam,
                               const FluxFunction *rf, const AccumulationFunction *ra, GridValues *rg, int rfam,
                               std::vector<RealVector> &left_curve, std::vector<RealVector> &right_curve);
};

#endif // _DOUBLE_CONTACT_GPU_

