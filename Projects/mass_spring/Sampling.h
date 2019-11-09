#pragma once

template<class T, int dim>
class Sampling {
public:
    using TV = Eigen::Matrix<T,dim,1>;
    /**
     * Stratefied sampling in a 3-dimensional grid.
     **/
    static std::vector<TV> stratefiedSampling(int &numSamples) {
        // TODO
        // The square root of the number of samples input
        int sqrtVal = (int) (std::sqrt((float) numSamples) + 0.5);
        // A number useful for scaling a square of size sqrtVal x sqrtVal to 1 x 1
        float invSqrtVal = 1.f / sqrtVal;

        std::vector<TV> samples;

        numSamples = sqrtVal * sqrtVal * sqrtVal;
        samples.resize(numSamples);

        // random number generator

        return std::vector<TV>();
    }
};