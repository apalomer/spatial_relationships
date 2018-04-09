#include <ceres/gradient_checker.h>

#include "compound_cost_function.h"
#include "inverse_compound_cost_function.h"
#include "tester.h"

#define MAX_ERROR 0.0005
#define ITERATIONS 10
#define DISPLAY_LEVEL DISPLAY_NONE
#define COMPUTE_INVERSION false
#define COMPUTE_COMPOUNDING false

void usage(char** argv)
{
    std::cout<<"Usage: "<<argv[0]<<" OPTIONS\n\n";
    std::cout<<"OPTIONS:\n";
    std::cout<<"  -h    display this message\n";
    std::cout<<"  -c    test compounding (+) operation\n";
    std::cout<<"  -i    test inverse compounding (-) operation\n";
    std::cout<<"  -m    max error (default: "<<MAX_ERROR<<")\n";
    std::cout<<"  -d    display level (default: "<<DISPLAY_LEVEL<<", available: all ("<<DISPLAY_ALL<<"), errors ("<<DISPLAY_ERROR<<") and none ("<<DISPLAY_NONE<<"))\n";
    std::cout<<"  -n    number of tests (default: "<<ITERATIONS<<")\n";
}

int main(int argc, char** argv)
{
    // Parse inputs
    int n(ITERATIONS);
    int display_level(DISPLAY_LEVEL);
    bool compute_inversion(COMPUTE_INVERSION);
    bool compute_compounding(COMPUTE_COMPOUNDING);
    double max_error(MAX_ERROR);
    int i(1);
    while (i<argc)
    {
        std::string aux = argv[i];
        if (aux == "-c")
        {
            compute_compounding = true;
        }
        else if (aux == "-i")
        {
            compute_inversion = true;
        }
        else if (aux == "-m")
        {
            if (++i < argc)
            {
                n = atof(argv[i]);
            }
            else
            {
                std::cout<<"Error -m flag has to be followed by the max error to pass the test"<<std::endl;
                usage(argv);
                return -1;
            }
        }
        else if (aux == "-d")
        {
            if (++i < argc)
            {
                display_level = atoi(argv[i]);
                if (display_level != DISPLAY_ALL && display_level != DISPLAY_ERROR && display_level != DISPLAY_NONE)
                {
                    std::cout<<"Not recognized display option\n";
                    usage(argv);
                    return -1;
                }
            }
            else
            {
                std::cout<<"Error -d flag has to be followed by the display option"<<std::endl;
                usage(argv);
                return -1;
            }
        }
        else if (aux == "-n")
        {
            if (++i < argc)
            {
                n = atoi(argv[i]);
            }
            else
            {
                std::cout<<"Error -n flag has to be followed by the number of tests to do"<<std::endl;
                usage(argv);
                return -1;
            }
        }
        else if (aux == "-h")
        {
            usage(argv);
            return 0;
        }
        else
        {
            std::cout<<"Unrecognized option: "<<aux<<"\n";
            usage(argv);
            return -1;
        }
        i++;
    }

    // Test Jacobians
    if (compute_inversion)
    {
        std::cout<<"*****************************************************************"<<std::endl;
        std::cout<<"*                                                               *"<<std::endl;
        std::cout<<"*                Testing inverse compounding function           *"<<std::endl;
        std::cout<<"*                                                               *"<<std::endl;
        std::cout<<"*****************************************************************"<<std::endl;
        Tester<InverseCompound3D,6,6> ti;
        int te = ti.test(n,max_error,display_level);
        std::cout<<"Total errors: "<<te<<"/"<<n<<std::endl;
    }

    // Test Jacobians
    if (compute_compounding)
    {
        std::cout<<"*****************************************************************"<<std::endl;
        std::cout<<"*                                                               *"<<std::endl;
        std::cout<<"*                Testing compounding function                   *"<<std::endl;
        std::cout<<"*                                                               *"<<std::endl;
        std::cout<<"*****************************************************************"<<std::endl;
        Tester<Compound3D,6,6,6> tc;
        int te = tc.test(n,max_error,display_level);
        std::cout<<"Total errors: "<<te<<"/"<<n<<std::endl;
    }

    // Exit
    return 0;
}
