#ifndef TESTER_H
#define TESTER_H

#include <ceres/ceres.h>
#include "functions.h"

/*!
 * Tester class for Compound3D and InverseCompound3D
 */
template<typename CostFunctionToProbe, int M = 0, int N0 = 0, int N1 = 0>
class Tester{
public:

    /*!
     * Display mode for the Tester class.
     */
    enum{
        DISPLAY_NONE,
        DISPLAY_ERROR,
        DISPLAY_ALL
    };

    /*!
     * \brief Constructor
     */
    Tester()
    {
        // Allocate residuals
        num_residuals = M;
        r_.resize(num_residuals,1);
        residuals_ = r_.data();

        // Check parameters
        if (N0 > 0)
        {
            parameters_size_.push_back(N0);
        }
        if (N1 > 0)
        {
            parameters_size_.push_back(N1);
        }

        // Allocate parameters and jacobians
        parameters_ = new double*[parameters_size_.size()];
        jacobians_automatic_ = new double*[parameters_size_.size()];
        jacobians_analytic_ = new double*[parameters_size_.size()];
        for (int i = 0;i<parameters_size_.size();i++)
        {
            // Parameters
            parameters_[i] = new double[parameters_size_[i]];

            // Jacobians
            ceres::Matrix J_autodiff;
            J_autodiff.resize(num_residuals,parameters_size_[i]);
            J_autodiff_.push_back(J_autodiff);
            jacobians_automatic_[i] = new double[num_residuals*parameters_size_[i]];
            ceres::Matrix J_analytic;
            J_analytic.resize(num_residuals,parameters_size_[i]);
            J_analytic_.push_back(J_analytic);
            jacobians_analytic_[i] = new double[num_residuals*parameters_size_[i]];
            ceres::Matrix e;
            e.resize(num_residuals,parameters_size_[i]);
            for (int j = 0;j<num_residuals;j++)
            {
                for (int k = 0 ;k<parameters_size_[i];k++)
                    e(j,k) = 0;
            }
            error_percentage_.push_back(e);
        }

        // Allocate cost functions
        cost_function_ =  new ceres::AutoDiffCostFunction<CostFunctionToProbe,M,N0,N1>(new CostFunctionToProbe);
        func_ = new CostFunctionToProbe;
    }

    ~Tester()
    {
        for (int i = 0;i<parameters_size_.size();i++)
        {
            delete[] parameters_[i];
            delete[] jacobians_analytic_[i];
            delete[] jacobians_automatic_[i];
        }
        delete[] jacobians_analytic_;
        delete[] jacobians_automatic_;
        delete[] parameters_;
        delete cost_function_;
        delete func_;
    }

    /*!
     * \brief test
     * \param iterations
     * \param max_error
     * \param display
     * \return
     */
    int test(int iterations, double max_error, int display = DISPLAY_ERROR)
    {
        int n_errors(0);
        for (int i = 0 ;i<iterations;i++)
        {
            // Create random x
            for (int j = 0;j<parameters_size_.size();j++)
            {
                Eigen::Matrix<double,6,1> aux;
                aux<<fRand(-100,100),fRand(-100,100),fRand(-100,100),fRand(-M_PI,M_PI),fRand(-M_PI,M_PI),fRand(-M_PI,M_PI);
                for (int k = 0;k<6;k++) parameters_[j][k] = aux(k,0);
            }

            // Evaluate automatic
            cost_function_->Evaluate(parameters_,residuals_,jacobians_automatic_);

            // Evaluate analytic
            func_->Evaluate(parameters_,residuals_,jacobians_analytic_);

            // Assign jacobians to matrices
            for(int j = 0;j<parameters_size_.size();j++)
            {
                for (int k = 0;k<num_residuals*parameters_size_[j];k++)
                {
                    J_analytic_[j](k/parameters_size_[j],k%parameters_size_[j]) = jacobians_analytic_[j][k];
                    J_autodiff_[j](k/parameters_size_[j],k%parameters_size_[j]) = jacobians_automatic_[j][k];
                }
            }

            // Compute difference
            std::vector<ceres::Matrix> err;
            std::vector<double> tej;
            bool berr(false);
            for (int j = 0;j<parameters_size_.size();j++)
            {
                ceres::Matrix e = J_analytic_[j] - J_autodiff_[j];
                err.push_back(e);
                double aux = e.lpNorm<Eigen::Infinity>();
                tej.push_back(aux);
                if (!berr)
                {
                    berr = aux > max_error;
                    if (berr) // only count the first
                    {
                        n_errors++;
                    }
                }
                if (berr)
                {
                    for (int rw = 0;rw<num_residuals;rw++)
                    {
                        for(int cl = 0;cl<parameters_size_[j];cl++)
                        {
                            if (fabs(e(rw,cl)) > max_error/10)
                                error_percentage_[j](rw,cl) = error_percentage_[j](rw,cl) + 1;
                        }
                    }
                }
            }

            if ((display == DISPLAY_ALL || (display == DISPLAY_ERROR && berr)) && display != DISPLAY_NONE)
            {
                std::cout<<"++++++++++++++++\nTest "<<i<<std::endl;
                std::cout<<"Evaluation:\n"<<r_.transpose()<<std::endl;
                for (int j = 0;j<parameters_size_.size();j++)
                {
                    std::cout<<"Analytic diff ("<<j<<"):\n"<<J_analytic_[j]<<std::endl;
                    std::cout<<"Automatic diff ("<<j<<"):\n"<<J_autodiff_[j]<<std::endl;
                    std::cout<<"Difference analytic-automatic ("<<j<<"):\n"<<err[j]<<std::endl;
                    std::cout<<"Jacobian error ("<<j<<"): "<<tej[j]<<std::endl;
                }
                std::cout<<"Test pass: "<<(!berr?"True":"False")<<std::endl;
            }
        }

        for (int j = 0;j<parameters_size_.size();j++)
        {
            for (int rw = 0;rw<num_residuals && n_errors > 0;rw++)
            {
                for(int cl = 0;cl<parameters_size_[j];cl++)
                {
                    error_percentage_[j](rw,cl) = error_percentage_[j](rw,cl)/n_errors;
                }
            }
            std::cout<<"Error per cell ("<<j<<"):\n"<<error_percentage_[j]<<std::endl;
        }

        return n_errors;
    }

protected:

    /*!
     * \brief parameters to be evaluated in the cost functions.
     */
    double** parameters_;

    /*!
     * \brief resutls of evaluating the cost function.
     */
    double* residuals_;

    /*!
     * \brief Jacobians computed using automatic differentiation from ceres.
     */
    double** jacobians_automatic_;

    /*!
     * \brief Jacobians computed usng analytic differentiation.
     */
    double** jacobians_analytic_;

    /*!
     * \brief Size of residual block.
     */
    int num_residuals;

    /*!
     * \brief Size of each parameter.
     */
    std::vector<int> parameters_size_;

    /*!
     * \brief Cost function to compute the automatic differentiation.
     */
    ceres::CostFunction* cost_function_;

    /*!
     * \brief Function to compute the analytic differentiation.
     */
    ceres::CostFunction* func_;

    /*!
     * \brief List of jacobians computed with automatic differentiation.
     */
    std::vector<ceres::Matrix> J_autodiff_;

    /*!
     * \brief List of jacobians computed with analytic differentiation.
     */
    std::vector<ceres::Matrix> J_analytic_;

    /*!
     * \brief Error for each element of the jacobians.
     */
    std::vector<ceres::Matrix> error_percentage_;

    /*!
     * \brief Residual expressed as a matrix (for easier display).
     */
    ceres::Matrix r_;

    /*!
     * \brief Parameters expressed as a matrix (for easier display).
     */
    ceres::Matrix p_;
};

#endif // TESTER_H
