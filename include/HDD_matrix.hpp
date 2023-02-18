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
    public :
        // Constructor for the matrix data structure
        HODLRdD_matrix(Kernel *userkernel, std::vector<ptsnD> *gPoints)
        {
            this->kernel_func = new kernel_function(userkernel);
            this->gridPoints = gPoints;
            N = gPoints->size();
            // Initialise the grid points to cluste
            // This initialises the grid points, this creates the hierarchical tre
            // Initialize the hierarchical 2^d tree
            //HODLRdD_tree = new Tree();
    }
    // Size of the matrix
    size_t get_size(){
        return N;
    }
    // Assemble the matrix operators 
    void Assemble_matrix_operators();
    // mat-vec product
    Vec mat_mat_product(Vec& x){
        return HODLRdD_tree->matvec(x);
    } 
    Vec operator*(Vec x)
    {
        Vec b = mat_mat_product(x);
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