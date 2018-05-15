#include "quaternion_distance_cost_function.h"

QuaternionDistance::QuaternionDistance()
{
}
QuaternionDistance::~QuaternionDistance()
{

}

bool QuaternionDistance::Evaluate(const double * const *parameters, double *residuals, double **jacobians) const
{
    // Compute cost
    operator ()(parameters[0],parameters[1],residuals);

    // Compute Jacobians
    if (jacobians != NULL && (jacobians[1] != NULL || jacobians[0] != NULL))
    {
        // Get parameter
        Eigen::Quaterniond q1(parameters[0][0],parameters[0][1],parameters[0][2],parameters[0][3]);
        Eigen::Quaterniond q2(parameters[1][0],parameters[1][1],parameters[1][2],parameters[1][3]);

        // Compute Jacobians
        Eigen::Matrix<double,1,4> J1;
        Eigen::Matrix<double,1,4> J2;
        quaternionDistanceJacobian(q1,q2,J1,J2);
        // Assign to output
        if (jacobians[0] != NULL)
        {
            for (int i = 0;i<4;i++)
            {
                jacobians[0][i] = J1(0,i);
            }
        }
        if (jacobians[1] != NULL)
        {
            for (int i = 0;i<4;i++)
            {
                jacobians[1][i] = J2(0,i);
            }
        }
    }

    // Exit
    return true;
}
