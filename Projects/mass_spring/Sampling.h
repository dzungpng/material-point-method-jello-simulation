#pragma once
#include "pcg32.h"

template<class T, int dim>
class Sampling {
public:
    using TV = Eigen::Matrix<T,dim,1>;
    /**
     * TODO: Stratefied sampling in a 3-dimensional grid.
     **/
    static std::vector<TV> stratefiedSampling(int &numSamples) {
        // The square root of the number of samples input
        int sqrtVal = (int) (std::sqrt((T) numSamples) + 0.5);
        // A number useful for scaling a square of size sqrtVal x sqrtVal to 1 x 1
        T invSqrtVal = 1.f / sqrtVal;

        std::vector<TV> samples;

        numSamples = sqrtVal * sqrtVal * sqrtVal;
        samples.resize(numSamples);

        // random number generator
        pcg32 rng;
        // rng.seed(42u, 54u);

        for (int i = 0; i < numSamples; i++) {
            TV sample; 

            int z = i % sqrtVal;
            int y = i % sqrtVal;    
            int x = i % sqrtVal;

            T randOffset_x = rng.nextFloat();
            T randOffset_y = rng.nextFloat();
            T randOffset_z = rng.nextFloat();

            T scaled_x = (T)x/(T)sqrtVal;
            T scaled_y = (T)y/(T)sqrtVal;
            T scaled_z = (T)z/(T)sqrtVal;

            float actualPlacement_x = scaled_x + (T)1/(T)sqrtVal * randOffset_x;
            float actualPlacement_y = scaled_y +  (T)1/(T)sqrtVal * randOffset_y;
            float actualPlacement_z = scaled_z + (T)1/(T)sqrtVal * randOffset_z;

            sample(0) = actualPlacement_x;
            sample(1) = actualPlacement_y;
            sample(2) = actualPlacement_z;

            samples.push_back(sample);
        }
        return samples;
    }

    /* 
    * compute 1D quadratic B spline weights 
    * x is assumed to be scaled in the index space (i.e., it is in a dx=1 grid)
    * w is 1x3 (row vector)
    */
    static void computeWeights1D(const T x, T &base_node, TV &wi) {
        base_node = floor(x-(T)0.5) + (T)1;
        wi.Zero();
        T d0 = x - base_node + (T)1;
        T z = (T)1.5 - d0;
        T z2 = z * z;
        wi(1) = (T)0.5 * z2;

        T d1 = d0 - (T)1;
        wi(2) = (T)0.75 - d1 * d1;

        T d2 = (T)1 - d1;
        T zz = (T)1.5 - d2;
        T zz2 = zz * zz;
        wi(3) = (T)0.5 * zz2;
    }

};