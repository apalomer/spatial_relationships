#include "inverse_compound_cost_function.h"

InverseCompound3D::InverseCompound3D()
{

}

InverseCompound3D::~InverseCompound3D()
{

}

bool InverseCompound3D::Evaluate(const double * const *parameters, double *residuals, double **jacobians) const
{

    // Compute cost
    operator ()(parameters[0],residuals);

    // Compute Jacobians
    if (jacobians != NULL && jacobians[0] != NULL){

        // Get transformation
        Eigen::Matrix<double,6,1> x;
        for (int i = 0;i<6;i++) x(i,0) = parameters[0][i];

        // Compute jacobian
        Eigen::Matrix<double,6,6> J = inverseCompound3DJacobian(x);

        // Assign to output
        for (int i = 0;i<36;i++)
        {
            jacobians[0][i] = J(i/6,i%6);
        }
    }

    // Return
    return true;
}
