#include "compound_cost_function.h"

Compound3D::Compound3D()
{

}

Compound3D::~Compound3D()
{

}

bool Compound3D::Evaluate(const double * const *parameters, double *residuals, double **jacobians) const
{

    // Compute cost
    operator ()(parameters[0],parameters[1],residuals);

    // Compute Jacobians
    if (jacobians != NULL && (jacobians[1] != NULL || jacobians[0] != NULL)){

        // Get transformation
        Eigen::Matrix<double,6,1> x1;
        Eigen::Matrix<double,6,1> x2;
        for (int i = 0;i<6;i++)
        {
            x1(i,0) = parameters[0][i];
            x2(i,0) = parameters[1][i];
        }

        // Compute jacobian
        Eigen::Matrix<double,6,6> J1;
        Eigen::Matrix<double,6,6> J2;
        compound3DJacobian(x1,x2,J1,J2);

        // Assign to output
        if (jacobians[0] != NULL)
        {
            for (int i = 0;i<36;i++)
            {
                jacobians[0][i] = J1(i/6,i%6);
            }
        }
        if (jacobians[1] != NULL)
        {
            for (int i = 0;i<36;i++)
            {
                jacobians[1][i] = J2(i/6,i%6);
            }
        }
    }

    // Return
    return true;
}
