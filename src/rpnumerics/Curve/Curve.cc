#include "Curve.h"


void Curve::init(const Curve &orig){
    type = orig.type;

    back_curve_index = orig.back_curve_index;

    family = orig.family;

    increase = orig.increase;

    reference_point = orig.reference_point;

    curve = orig.curve;
    xi = orig.xi;
    back_pointer = orig.back_pointer;
    back_curve_pointer = orig.back_curve_pointer;

    last_point = orig.last_point;
    final_direction = orig.final_direction;

    reason_to_stop = orig.reason_to_stop;

    speed = orig.speed;
    eigenvalues = orig.eigenvalues;
    extension_of_whom = orig.extension_of_whom;
    explicit_bifurcation_transition_index = orig.explicit_bifurcation_transition_index;

    compatible = orig.compatible;

    return;
}

Curve::Curve(){
    back_curve_index = -1;
}

Curve::Curve(const Curve &orig){
    init(orig);
}

Curve::Curve(const Curve *orig){
    init(*orig);
}

Curve::~Curve(){
}

Curve Curve::operator=(const Curve &orig){
    // Avoid self-assignment
    if (this != &orig) init(orig);

    return *this;
}

void Curve::clear(){
    curve.clear();
    speed.clear();
    eigenvalues.clear();
    back_pointer.clear();
    back_curve_pointer.clear();
    extension_of_whom.clear();
    explicit_bifurcation_transition_index.clear();
    compatible.clear();

    return;             
}

bool Curve::segment_intersects_box(const RealVector &seg0, const RealVector &seg1, const RealVector &b0, const RealVector &b1, int side_to_check, int &side){
    // These below are useless.
    //
    RealVector r;
    double alpha, beta;

    RealVector box0(2), box1(2);
    for (int i = 0; i < 2; i++){
        box0(i) = std::min(b0(i), b1(i));
        box1(i) = std::max(b0(i), b1(i));
    }

    // First test.
    //
    RealVector d0 = box0;
    RealVector d1(2); d1(0) = box1(0); d1(1) = box0(1);

    if (side_to_check != ALL_SIDES_EXCEPT_DOWN){
        if (segment_segment_intersection(seg0, seg1, d0, d1, r, alpha, beta)){
            side = DOWN;

            return true;
        }
    }

    // Second test.
    //
    RealVector d2 = box1;
    if (side_to_check != ALL_SIDES_EXCEPT_RIGHT){
        if (segment_segment_intersection(seg0, seg1, d1, d2, r, alpha, beta)){
            side = RIGHT;

            return true;
        }
    }

    // Third test.
    //
    RealVector d3(2); d3(0) = box0(0); d3(1) = box1(1);
    if (side_to_check != ALL_SIDES_EXCEPT_UP){
        if (segment_segment_intersection(seg0, seg1, d2, d3, r, alpha, beta)){
            side = UP;

            return true;
        }
    }

    // Fourth test.
    //
    if (side_to_check != ALL_SIDES_EXCEPT_LEFT){
        if (segment_segment_intersection(seg0, seg1, d3, d0, r, alpha, beta)){
            side = LEFT;

            return true;
        }
    }

    side = NONE;

    return false;
}

bool Curve::segment_intersects_box(const RealVector &seg0, const RealVector &seg1, const RealVector &box0, const RealVector &box1){
    int side;
    
    return segment_intersects_box(seg0, seg1, box0, box1, ALL_SIDES, side);
}

void Curve::disable_adjacent_cells(GridValues *g, int row, int col){
    int shift = 2; // Should be 1.

    int imin = std::max(row - shift, 0);
    int imax = std::min(row + shift, g->cell_type.rows() - 1);
        
    int jmin = std::max(col - shift, 0);
    int jmax = std::min(col + shift, g->cell_type.cols() - 1);

//    std::cout << "        g->cell_type.rows() = " << g->cell_type.rows() << std::endl;
//    std::cout << "        g->cell_type.cols() = " << g->cell_type.cols() << std::endl;
//    std::cout << "        row = " << row << ", imin = " << imin << ", imax = " << imax << std::endl;
//    std::cout << "        col = " << col << ", jmin = " << jmin << ", jmax = " << jmax << std::endl;

    for (int ii = imin; ii < imax; ii++){
        for (int jj = jmin; jj < jmax; jj++){
            g->cell_type(ii, jj) = CELL_IS_INVALID;
        }
    }

//    g->cell_type(row, col) = CELL_IS_INVALID;

    return;
}

void Curve::disable_intersecting_cells(GridValues *g){
    std::cout << "curve.size() = " << curve.size() << std::endl;

    for (int k = 0; k < curve.size() - 1; k++){
//        for (int i = 0; i < g->cell_type.rows(); i++){
//            for (int j = 0; j < g->cell_type.cols(); j++){
////                std::cout << "k = " << k << ", i = " << i << ", j = " << j << std::endl;

//                if (segment_intersects_box(curve[k], curve[k + 1], g->grid(i, j), g->grid(i + 1, j + 1))){
//                    g->cell_type(i, j) = CELL_IS_INVALID;
//                }
//            }
//        }

        int row, col;
        if (g->cell(curve[k], row, col)){
//            g->cell_type(row, col) = CELL_IS_INVALID;
//            std::cout << "    k = " << k << ", row = " << row << ", col = " << col << std::endl;
//            disable_adjacent_cells(g, row, col);

            disable_intersecting_cells(curve[k], curve[k + 1], g);
        }
    }

    return;
}



void Curve::disable_intersecting_cells(const RealVector &p, const RealVector &q, GridValues *g){
    int ray_row, ray_col;
    g->cell(p, ray_row, ray_col);


    std::cout << "p = " << p << std::endl;
    std::cout << "    row = " << ray_row << ", col = " << ray_col << std::endl;

    int q_row, q_col;
    g->cell(q, q_row, q_col);
    std::cout << "q = " << q << std::endl;
    std::cout << "    row = " << q_row << ", col = " << q_col << std::endl;

    int side_to_check = ALL_SIDES;

    while (ray_row != q_row || ray_col != q_col){
        g->cell_type(ray_row, ray_col) = CELL_IS_INVALID;

        // This assumes that the segment is contained in the grid due to the +1's below. Perhaps that could be changed.
        //
        int side;
        bool info = segment_intersects_box(p, q, g->grid(ray_row, ray_col), g->grid(ray_row + 1, ray_col + 1), side_to_check, side);
        if (!info) return;

        std::cout << "side = " << side << std::endl;

        if (side == LEFT){
            side_to_check = ALL_SIDES_EXCEPT_RIGHT;
            ray_col--;
        }
        else if (side == RIGHT){
            side_to_check = ALL_SIDES_EXCEPT_LEFT;
            ray_col++;
        }
        else if (side == UP){
            side_to_check = ALL_SIDES_EXCEPT_DOWN;
            ray_row--;
        }
        else if (side == DOWN){
            side_to_check = ALL_SIDES_EXCEPT_UP;
            ray_row++;
        }
        else { // NONE
            return;
        }

//        TestTools::pause();
    }

    g->cell_type(q_row, q_col) = CELL_IS_INVALID;

    return;
}

