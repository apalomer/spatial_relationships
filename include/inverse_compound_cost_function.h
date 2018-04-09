#ifndef INVERSE_COMPOUND_COST_FUNCTION_H
#define INVERSE_COMPOUND_COST_FUNCTION_H

#include <ceres/ceres.h>

#include "transformtions.h"

/*!
 * \brief The InverseCompound3D class
 */
class InverseCompound3D: public ceres::SizedCostFunction<6,6>{
public:

    /*!
     * \brief Constructor
     */
    InverseCompound3D();

    ~InverseCompound3D();

    /*!
     * \brief Evaluate
     * \param parameters
     * \param residuals
     * \param jacobians
     * \return
     */
    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    template<typename T>
    bool operator()(const T* const parameter, T* residuals) const
    {
        // Get transformation
        Eigen::Matrix<T,6,1> x;
        for (int i = 0;i<6;i++) x(i,0) = parameter[i];

        // Compute inversion
        Eigen::Matrix<T,6,1> xi = inverseCompound3D<T>(x);

        // Copy to residuals
        for (int i = 0;i<6;i++) residuals[i] = xi(i,0);

        // Exit
        return true;
    }
};

#endif // INVERSE_COMPOUND_COST_FUNCTION_H
