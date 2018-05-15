#ifndef QUATERNIONFROMRPY_H
#define QUATERNIONFROMRPY_H

#include <ceres/ceres.h>

#include "transformations.h"

/*!
 * \brief The Quaternion from RPY cost function
 * \callgraph
 * \callergraph
 */
class QuaternionFromRPY: public ceres::SizedCostFunction<4,3>
{
public:

    /*!
     * \brief Constructor
     * \callgraph
     * \callergraph
     */
    QuaternionFromRPY();

    ~QuaternionFromRPY();

    /*!
     * \brief Evaluate the cost function and its jacobians with the given parameters
     * \param parameters
     * \param[out] residuals
     * \param[out] jacobians
     * \return ture/false if the residuals and the jacobians have been properly evaluated.
     * \callgraph
     * \callergraph
     */
    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    /*!
     * \brief operator() evaluates the cost function for the given parameter.
     * \param parameter
     * \param[out] residuals xi = (-) parameter
     * \return ture/false if the residuals and the jacobians have been properly evaluated.
     * \callgraph
     * \callergraph
     */
    template<typename T>
    bool operator()(const T* const parameter, T* residuals) const
    {
        // Get quaternion
        Eigen::Matrix<T,3,1> rpy;
        for (int i = 0;i<3;i++) rpy(i,0) = parameter[i];

        // Compute inversion
        Eigen::Quaternion<T> q = quaternionFromRPY<T>(rpy);

        // Copy to residuals
        residuals[0] = q.w();
        residuals[1] = q.x();
        residuals[2] = q.y();
        residuals[3] = q.z();

        // Exit
        return true;
    }
};

#endif // QUATERNIONFROMRPY_H
