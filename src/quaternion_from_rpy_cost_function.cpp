#include "quaternion_from_rpy_cost_function.h"

QuaternionFromRPY::QuaternionFromRPY()
{
}

QuaternionFromRPY::~QuaternionFromRPY()
{

}

bool QuaternionFromRPY::Evaluate(const double * const *parameters, double *residuals, double **jacobians) const
{
    // Compute cost
    operator ()(parameters[0],residuals);

    // Compute Jacobians
    if (jacobians != NULL && jacobians[0] != NULL)
    {

        Eigen::Matrix<double,3,1> rpy;
        for (int i = 0;i<3;i++)
        {
            rpy(i,0) = parameters[0][i];
        }

        Eigen::Matrix<double,4,3> J = quaternionFromRPYJacobian(rpy);
        for (int i = 0;i<12;i++)
        {
            jacobians[0][i] = J(i/3,i%3);
        }
    }
    return true;
}
