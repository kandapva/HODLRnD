#ifndef __HDD_matrix__
#define __HDD_matrix__

#include "myHeaders.hpp"

// TODO : Rectangular matrix structure ?? 
class HODLRdD_matrix{
size_t N;
Tree *HODLRdD_tree;
public:
// Constructor for the matrix data structure
    HODLRdD_matrix(){
        // Initialise the grid points to cluster

        
        // Initialize the hierarchical 2^d tree 
    }
    // Size of the matrix
    size_t get_size(){
        return N;
    }
    // mat-vec product
    
    // solve using GMRES
    // Destructor
    ~HODLRdD_matrix(){
    }
};
#endif