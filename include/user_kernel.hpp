#ifndef __user_kernel__
#define __user_kernel__

#include "myHeaders.hpp"
#include "points_dt.hpp"
#include "HDD_clusters.hpp"
#include "../integral_equation_4D/integral4d.hpp"

// Test Greens' function in 2d, 3d and 4d
class kernel_2d_test{
    double c = -1.0 / (2 * PI);
public:
    // location of particles in the domain
    std::vector<ptsnD> *gridPoints;
    kernel_2d_test()
    {
        gridPoints = new std::vector<ptsnD>;
        // Nodes in one dimension # Test considers a tensor grid in NDIM using loc_dir
        VectorXd loc_dir[NDIM];
        Eigen::VectorXd X(NDIM), Y(NDIM);
        for (int i = 0; i < NDIM; ++i)
        {
            X(i) = -1;
            Y(i) = 1;
        }
        int cts = 0;
        for (int i = 0; i < NDIM; ++i)
            loc_dir[i] = uniform_nodes(X(i), Y(i), numPoints);
        // Forming the tensor product grid of the points
        for (int i = 0; i < numPoints; i++)
        {
            for (int j = 0; j < numPoints; j++)
            {
                ptsnD temp;
                temp.x[0] = loc_dir[0][i];
                temp.x[1] = loc_dir[1][j];
                temp.id = cts++;
                gridPoints->push_back(temp);
            }
        }
    }
    void get_points(std::vector<ptsnD>*& src){
        src = this->gridPoints;
    }
    // The Green's function in 2D
    dtype_base Kernel_Fun(dtype_base x)
    {
            return c * log(x);
    }
    dtype_base getMatrixEntry(int i, int j)
    {
            double r;
            if (i == j)
                r = 0.0;
            else
                r = Kernel_Fun(nd_points::euclidean_distance(gridPoints->at(i), gridPoints->at(j))); 
            return r;
    }
    ~kernel_2d_test() {}
};

class kernel_3d_test{
    double c = 1.0 / (4 * PI);
public:
    std::vector<ptsnD> *gridPoints;
    kernel_3d_test()
    {
            // location of particles in the domain
            gridPoints = new std::vector<ptsnD>;
            // Nodes in one dimension # Test considers a tensor grid in NDIM using loc_dir
            VectorXd loc_dir[NDIM];
            Eigen::VectorXd X(NDIM), Y(NDIM);
            for (int i = 0; i < NDIM; ++i)
            {
                X(i) = -1;
                Y(i) = 1;
            }
            int cts = 0;
            for (int i = 0; i < NDIM; ++i)
                loc_dir[i] = uniform_nodes(X(i), Y(i), numPoints);
            // Forming the tensor product grid of the points
            for (int i = 0; i < numPoints; i++)
            {
                for (int j = 0; j < numPoints; j++)
                {
                    for (int k = 0; k < numPoints; k++){
                    ptsnD temp;
                    temp.x[0] = loc_dir[0][i];
                    temp.x[1] = loc_dir[1][j];
                    temp.x[2] = loc_dir[2][k];
                    temp.id = cts++;
                    gridPoints->push_back(temp);
                    }
                }
            }
    }
    void get_points(std::vector<ptsnD> *&src){
            src = this->gridPoints;
    }
    // The Green's function in 3D
    dtype_base Kernel_Fun(dtype_base x)
    {   
            return c/x;
    }
    dtype_base getMatrixEntry(int i, int j)
    {
            double r;
            if (i == j)
                r = 0.0;
            else
                r = Kernel_Fun(nd_points::euclidean_distance(gridPoints->at(i), gridPoints->at(j)));
            return r;
    }
    ~kernel_3d_test() {}
};

class kernel_4d_test{
    double c = -1.0 / (4 * PI * PI);
    double h = 1.0 / numPoints;
    double h4 = pow(h, 4.0);
    double kii = 0.0;
public:
    std::vector<ptsnD> *gridPoints;
    kernel_4d_test()
    {
            // location of particles in the domain
            double *a = new double[4];
            double *b = new double[4];
            a[0] = 0, a[1] = 0, a[2] = 0, a[3] = 0;
            b[0] = h * 0.5, b[1] = h * 0.5, b[2] = h * 0.5, b[3] = h * 0.5;
            kii = 1.0 + quadruple_integral(a, b); // Second kind
            kii /= h4;
            gridPoints = new std::vector<ptsnD>;
            // Nodes in one dimension # Test considers a tensor grid in NDIM using loc_dir
            VectorXd loc_dir[NDIM];
            Eigen::VectorXd X(NDIM), Y(NDIM);
            for (int i = 0; i < NDIM; ++i)
            {
                X(i) = -1;
                Y(i) = 1;
            }
            int cts = 0;
            for (int i = 0; i < NDIM; ++i)
                loc_dir[i] = uniform_nodes(X(i), Y(i), numPoints);
            // Forming the tensor product grid of the points
            for (int i = 0; i < numPoints; i++)
            {
                for (int j = 0; j < numPoints; j++)
                {
                    for (int k = 0; k < numPoints; k++){
                        for (int l = 0; l < numPoints; l++){
                            ptsnD temp;
                            temp.x[0] = loc_dir[0][i];
                            temp.x[1] = loc_dir[1][j];
                            temp.x[2] = loc_dir[2][k];
                            temp.x[3] = loc_dir[3][l];
                            temp.id = cts++;
                            gridPoints->push_back(temp);
                        }
                    }
                }
            }
    }
    void get_points(std::vector<ptsnD>*& src){
            src = this->gridPoints;
    }
    // The Green's function in 4D
    dtype_base Kernel_Fun(dtype_base x)
    {
            return c/(x*x);
    }
    dtype_base getMatrixEntry(int i, int j)
    {
            if (i == j)
                return kii;
            else
            {
                return Kernel_Fun(nd_points::euclidean_distance(gridPoints->at(i), gridPoints->at(j)));
            }
                
    }
    ~kernel_4d_test(){}
};


#endif