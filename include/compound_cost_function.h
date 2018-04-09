#ifndef COMPOUND_COST_FUNCTION_H
#define COMPOUND_COST_FUNCTION_H


#include <ceres/ceres.h>

#include "transformtions.h"

/*!
 * \brief The InverseCompound3D class
 */
class Compound3D: public ceres::SizedCostFunction<6,6,6>{
public:

    /*!
     * \brief Constructor
     */
    Compound3D();

    ~Compound3D();

    /*!
     * \brief Evaluate
     * \param parameters
     * \param residuals
     * \param jacobians
     * \return
     */
    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    template<typename T>
    bool operator()(const T* const x1, const T* const x2, T* residuals) const
    {

        // Create the two transforms
        Eigen::Matrix<T,6,1> x1m;
        Eigen::Matrix<T,6,1> x2m;
        for (int i = 0;i<6;i++)
        {
            x1m(i,0) = x1[i];
            x2m(i,0) = x2[i];
        }

        // Compound
        Eigen::Matrix<T,6,1> x3m = compound3D<T>(x1m,x2m);

        // Retrieve result
        for (int i = 0;i<6;i++) residuals[i] = x3m[i];

        // Exit
        return true;
    }
};


#endif // COMPOUND_COST_FUNCTION_H
