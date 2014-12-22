#ifndef _WAVECURVE_
#define _WAVECURVE_

#include "Curve.h"

#define WAVECURVE_POSITIVE_ARCLENGTH 1
#define WAVECURVE_NEGATIVE_ARCLENGTH (-1)

class WaveCurve {
    private:
    protected:
        void half_speed_map(int begin, int end, std::vector<double> &arc_length, std::vector<double> &speed, std::vector<RealVector> &eigenvalues) const;

        void half_speed_map(int type, int begin, int end, std::vector<Curve> &arclength_speed, std::vector<Curve> &arclength_eigenvalues) const;
    public:
        int family;

        int increase;

        ReferencePoint reference_point;

        std::vector<Curve> wavecurve;

        int beginnig_of_second_half;

        void init(const WaveCurve &orig);

        WaveCurve();
        WaveCurve(const WaveCurve &orig);
        WaveCurve(const WaveCurve *orig);
        ~WaveCurve();

        WaveCurve operator=(const WaveCurve &orig);

        void add(const Curve &c);

        // Insert a rarefaction-composite pair in the wavecurve BEFORE composite_insertion_point_index (mimicking std::vector::insert).
        //
        void insert_rarefaction_composite_pair(int composite_curve_index, int composite_insertion_point_index, const RealVector &rar_point, const RealVector &cmp_point);
        void insert_rarefaction_composite_pair(int composite_curve_index, int composite_insertion_point_index, const RealVector &rarcmp_point);

        // Validate the compatibility between the wavespeeds of two states in the wavecurves.
        //
        void validate_compatibility(WaveCurve &other_wavecurve, int curve_index, int point_index, int direction);

        // Clear a wavecurve.
        //
        void clear();

        // Speed map.
        //
        void speed_map(std::vector<double> &arc_length, std::vector<double> &speed, std::vector<RealVector> &eigenvalues, RealVector &reference_eigenvalues) const;

        void speed_map(std::vector<Curve> &arclength_speed, std::vector<Curve> &arclength_eigenvalues, std::vector<Curve> &arclength_reference_eigenvalues) const;
};

#endif // _WAVECURVE_

