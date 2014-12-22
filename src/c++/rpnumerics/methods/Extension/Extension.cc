#include "Extension.h"

void Extension::divide_curve(void *obj, bool (*f)(void*, const RealVector&), const Curve &original, std::vector<Curve> &out){
    out.clear();

    // Failsafe.
    //
    if (f == 0){
        out.push_back(original);

        return;
    }

    int n = original.curve.size();
    int pos = 0;
    while (pos < original.curve.size()){
        Curve temp_curve;

        // Find the first point of this part of the curve
        // withing the specified domain.
        //
        while (pos < n && !(*f)(obj, original.curve[pos])) pos++;

        if (pos >= n) break;

        while (pos < n && (*f)(obj, original.curve[pos])){
            temp_curve.curve.push_back(original.curve[pos]);
            temp_curve.back_pointer.push_back(original.back_pointer[pos]);
            temp_curve.speed.push_back(original.speed[pos]);

            pos++;
        }

        if (temp_curve.curve.size() > 0){
            temp_curve.type = original.type;
            out.push_back(temp_curve);
        }
    }

    return;
}

void Extension::extension_curve(const Curve &curve, 
                                void *fdo, bool (*fd)(void*, const RealVector&), 
                                void *fio, bool (*fi)(void*, const RealVector&), 
                                std::vector<Curve> &ext_curve){

    ext_curve.clear();

    // Find all the points that are within the domain. More than one curve may be
    // obtained, because the original curve may leave and enter the domain several
    // times.
    //
    std::vector<Curve> domain_curve;
    divide_curve(fdo, fd, curve, domain_curve);

    // Extend.
    //
    for (int i = 0; i < domain_curve.size(); i++){
        Curve temp;

        for (int j = 0; j < domain_curve[i].curve.size(); j++){
            RealVector ext_temp;
            int info = extension(domain_curve[i].curve[j], ext_temp);
            if (info == EXTENSION_OK){
                temp.curve.push_back(ext_temp);
                temp.speed.push_back(domain_curve[i].speed[j]);
                temp.back_pointer.push_back(domain_curve[i].back_pointer[j]);
            }
        }

        if (temp.curve.size() > 0){
            std::vector<Curve> temp_ext_curve;
            divide_curve(fio, fi, temp, temp_ext_curve);
            
            for (int k = 0; k < temp_ext_curve.size(); k++){
                // TODO: The speeds are missing. Create array of index at divide curve.
                temp_ext_curve[k].type = COMPOSITE_CURVE;
                ext_curve.push_back(temp_ext_curve[k]);
            }
        }
    }

    return;
}

void Extension::extension_curve(const std::vector<RealVector> &curve, const std::vector<RealVector> &domain_polygon, const std::vector<RealVector> &image_polygon, std::vector<RealVector> &ext_curve){
    // Find the convex hulls.
    //
    std::vector<RealVector> domain_convex_hull;
    if (domain_polygon.size() > 2) convex_hull(domain_polygon, domain_convex_hull);

    std::vector<RealVector> image_convex_hull;
    if (image_polygon.size() > 2) convex_hull(image_polygon, image_convex_hull);

    // Find the points of the curve that are within the domain.
    //
    std::vector<RealVector> final_curve;

    if (domain_convex_hull.size() > 0){
        for (int i = 0; i < curve.size(); i++){
            if (inside_convex_polygon(domain_convex_hull, curve[i])) final_curve.push_back(curve[i]);
        }
    }
    else final_curve = curve;

    RealVector ext_p;

    ext_curve.clear();
    for (int i = 0; i < final_curve.size(); i++){
        if (extension(final_curve[i], ext_p) == EXTENSION_OK){
            if (image_convex_hull.size() > 0){
                if (inside_convex_polygon(image_convex_hull, ext_p)) ext_curve.push_back(ext_p);
            }
            else ext_curve.push_back(ext_p);
        }
    }                

    return;
}

void Extension::extension_curve(const Curve &curve, const std::vector<RealVector> &domain_polygon, const std::vector<RealVector> &image_polygon, Curve &ext_curve){
    extension_curve(curve.curve, domain_polygon, image_polygon, ext_curve.curve);

    for (int i = 0; i < curve.speed.size(); i++) ext_curve.speed.push_back(curve.speed[i]);

    ext_curve.type = COMPOSITE_CURVE;

    return;
}

