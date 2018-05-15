#ifndef QUATERNION_ANGLE_COST_FUNCTION_H
#define QUATERNION_ANGLE_COST_FUNCTION_H

#include <ceres/ceres.h>

#include "transformations.h"
#include "functions.h"

/*!
 * \brief The QuaternionAngle cost function class
 * \callgraph
 * \callergraph
 */
class QuaternionAngle: public ceres::SizedCostFunction<1,4,4>
{
public:

    /*!
     * \brief Constructor
     */
    QuaternionAngle();

    ~QuaternionAngle();

    /*!
     * \brief Evaluate the cost function with the given parameters
     * \param parameters
     * \param[out] residuals
     * \param[out] jacobians
     * \return ture/false if the residuals and the jacobians have been properly evaluated.
     * \callgraph
     * \callergraph
     */
    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;

    /*!
     * \brief operator() evaluates the cost function for the given parameters.
     * \param q1
     * \param q2
     * \param[out] residuals = quaternionAngle(q1,q2)
     * \return ture/false if the residuals and the jacobians have been properly evaluated.
     * \callgraph
     * \callergraph
     */
    template<typename T>
    bool operator()(const T* const q1, const T* const q2, T* residuals) const
    {
        Eigen::Quaternion<T> q_1(q1[0],q1[1],q1[2],q1[3]);
        Eigen::Quaternion<T> q_2(q2[0],q2[1],q2[2],q2[3]);
        residuals[0] = quaternionAngle<T>(q_1,q_2);
        return true;
    }
};

#endif // QUATERNION_ANGLE_COST_FUNCTION_H
