#ifndef __HDD_matrix__
#define __HDD_matrix__

#include "myHeaders.hpp"
#include "kernel_function.hpp"
#include "HDD_tree.hpp"

// TODO : Rectangular matrix structure ??
template <class Kernel>
class HODLRdD_matrix
{
    size_t N;
    Tree *HODLRdD_tree; // this is private to the matrix
    kernel_function<Kernel> *kernel_func;
    std::vector<ptsnD> *gridPoints;
    cluster *src;
public:
    // Constructor for the matrix data structure. Also, x1 and x2 is the bounding box of the cluster of points 
    HODLRdD_matrix(Kernel*& userkernel, std::vector<ptsnD>*& gPoints, Eigen::VectorXd x1, Eigen::VectorXd x2)
    {
        this->kernel_func = new kernel_function(userkernel);
        this->gridPoints = gPoints;
        N = gPoints->size();
        // Source cluster needs special treatment i.e., provide a guide list that maps the points
        src = new cluster(x1, x2, this->gridPoints);
        // Initialise the grid points to cluste
        // This initialises the grid points, this creates the hierarchical tre
        // Initialize the hierarchical 2^d tree
        HODLRdD_tree = new Tree(src, gPoints, kernel_func);
    }
    // Size of the matrix
    size_t get_size(){
        return N;
    }
    // Assemble the matrix operators 
    void Assemble_matrix_operators(){
        HODLRdD_tree->Initialise_tree();
    }
    // mat-vec product or mat-mat product
    Vec operator*(Vec x)
    {
        Vec b = HODLRdD_tree->matvec(x);
        return b;
    }

    // Assemble the operators and access to the Tree structure
    
    // solve using GMRES
    Vec solve(const Vec& b){
        Vec x = Vec::Zero(b.size());
        // Routine to compute the GMRES iterations
        return x;
    }
    // Destructor
    ~HODLRdD_matrix(){
    }
};
#endif