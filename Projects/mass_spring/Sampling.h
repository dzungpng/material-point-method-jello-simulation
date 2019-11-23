#pragma once
#include "pcg32.h"

template<class T, int dim>
class Sampling {
public:
    using TV = Eigen::Matrix<T,dim,1>;

    /*
    ** Select only particles within the max and min corners from the input ogSamples. 
    */
    static void selectInBox(const std::vector<TV> ogSamples, std::vector<TV> &newSamples,
                            const TV minCorner, const TV maxCorner) {
        for(int i = 0; i < ogSamples.size(); i++) {
            if(!(ogSamples[i](0) < minCorner(0) || ogSamples[i](0) > maxCorner(0) ||
                ogSamples[i](1) < minCorner(1) || ogSamples[i](1) > maxCorner(1) ||
                ogSamples[i](2) < minCorner(2) || ogSamples[i](2) > maxCorner(2))) {
                    newSamples.push_back(ogSamples[i]);
            }
        }
    }

    /* 
    * compute 1D quadratic B spline weights 
    * x is assumed to be scaled in the index space (i.e., it is in a dx=1 grid)
    * w is 1x3 (row vector)
    */
    static void computeWeights1D(const T x, T &base_node, TV &wi) {
        base_node = floor(x-0.5) + 1;
        
        T d0 = x - base_node + (T)1;
        T z = (T)1.5 - d0;
        T z2 = z * z;
        wi(0) = (T)0.5 * z2;

        T d1 = d0 - (T)1;
        wi(1) = (T)0.75 - d1 * d1;

        T d2 = (T)1 - d1;
        T zz = (T)1.5 - d2;
        T zz2 = zz * zz;
        wi(2) = (T)0.5 * zz2;
    }

    // compute 1D quadratic B spline weights 
    // x is assumed to be scaled in the index space (i.e., it is in a dx=1 grid)
    // w is 1x3 (row vector)
    static void computeWeightsWithGradients1D(const T x, TV &w, TV &dw, T &base_node) {
        base_node = floor(x-0.5) + 1;
        w = TV::Zero();
        dw = TV::Zero();

        T d0 = x - base_node + (T)1;
        T z = (T)1.5 - d0;
        T z2 = z * z;
        w(0) = (T)0.5 * z2;

        T d1 = d0 - (T)1;
        w(1) = (T)0.75 - d1 * d1;

        T d2 = (T)1 - d1;
        T zz = (T)1.5 - d2;
        T zz2 = zz * zz;
        w(2) = (T)0.5 * zz2;

        dw(0) = -z;
        dw(1) = (T)-2*d1;
        dw(2) = zz;
    }

};