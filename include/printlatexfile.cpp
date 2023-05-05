#include <cstdlib>
#include "HDD_matrix.hpp"
#include "HDD_clusters.hpp"
#include "kernel_function.hpp"
#include "myHeaders.hpp"
#include "points_dt.hpp"
#include "LowRank_matrix.hpp"
#include "user_kernel.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    std::vector<ptsnD>* gridPoints;
    if(argc>1){
        numPoints = atoi(argv[1]);
        INTERACTION_TYPE_ALLOWED = atoi(argv[2]);
    }
    else
        numPoints = 10;
    N = pow(numPoints, NDIM);
    // Bounding box of the kernel
    Eigen::VectorXd X(NDIM), Y(NDIM);
    for (int i = 0; i < NDIM; ++i)
    {
        X(i) = -1;
        Y(i) = 1;
    }

    kernel_4d_test *ker = new kernel_4d_test();
    ker->get_points(gridPoints);
    std::cout << "Kernel formed" << std::endl;
    HODLRdD_matrix Kmat = HODLRdD_matrix(ker, gridPoints, X, Y);
    // Kmat.Assemble_matrix_operators();
    std::cout << "The size of K matrix " << Kmat.get_size() << std::endl;
    Kmat.print_matrix_latex();


    return 0;
}


