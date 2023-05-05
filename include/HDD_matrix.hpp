#ifndef __HDD_matrix__
#define __HDD_matrix__

#include "myHeaders.hpp"
#include "kernel_function.hpp"
#include "HDD_tree.hpp"
#include "GMRES.hpp"

// TODO : Rectangular matrix structure ??
template <class Kernel>
class HODLRdD_matrix
{
    size_t N;
    Tree<Kernel> *HODLRdD_tree; // this is private to the matrix
    kernel_function<Kernel> *kernel_func;
    std::vector<ptsnD> *gridPoints;
    cluster *src;
    double solve_time = 0.0;
public:
    // Constructor for the matrix data structure. Also, x1 and x2 is the bounding box of the cluster of points 
    HODLRdD_matrix(Kernel*& userkernel, std::vector<ptsnD>*& gPoints, Eigen::VectorXd x1, Eigen::VectorXd x2)
    {
        std::cout << "HMATRIX START" << std::endl;
        this->kernel_func = new kernel_function<Kernel>(userkernel);
        this->gridPoints = gPoints;
        N = gPoints->size();
        // Source cluster needs special treatment i.e., provide a guide list that maps the points
        src = new cluster(x1, x2, this->gridPoints);
        // Initialise the grid points to cluste
        // This initialises the grid points, this creates the hierarchical tre
        // Initialize the hierarchical 2^d tree
        HODLRdD_tree = new Tree<Kernel>(src, gPoints, kernel_func);
        std::cout<< "HMATRIX DONE" << std::endl; 
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
    Vec operator * (Vec x)
    {
        Vec b = HODLRdD_tree->mat_vec(x);
        return b;
    }

    // Assemble the operators and access to the Tree structure
    // solve using GMRES
    Vec solve(Vec& b){
        Vec x = Vec::Zero(b.size());
        x = b;
        HODLRdD_matrix<Kernel> *mat_obj;
        mat_obj = this;
        iterSolver<HODLRdD_matrix<Kernel> > solve_obj(100, N, eps_ACA);
        std::string gmres_history = "gmres_"+ std::to_string(N) +".txt";
        solve_obj.set_output_file(gmres_history);
        double start = omp_get_wtime();
        int k = solve_obj.GMRES(mat_obj, x, b);
        start = omp_get_wtime() - start;
        solve_time = start;
        std::cout << "GMRES took " << start << " with status " << k << std::endl;
        // Routine to compute the GMRES iterations
        return x;
    }
    void print_matrix_details(){
        std::cout << "============HODLRnD=================" << std::endl;
        HODLRdD_tree->print_tree_details();
        std::cout << "====================================" << std::endl;
    }

    void print_matrix_latex(){
        HODLRdD_tree->print_matrix_latex();
    }

    void print_matrix_statistics()
    {
        double n_FLOP = 0.0;
        size_t MAX_RANK = 0;
        double compression_ratio = 1.0/(double(N)*double(N));
        // 1 double is considered to consume 8 bytes of physical memory
        double memory_gb = 8 * pow(10,-9);
        HODLRdD_tree->get_stat(n_FLOP, MAX_RANK);
        compression_ratio *= n_FLOP;
        memory_gb *= n_FLOP;
        std::cout << getTimeStamp() << std::endl;
        std::cout << "Memory (in GB) : " << memory_gb << std::endl;
        std::cout << "Number of FLOP : " << n_FLOP << std::endl;
        std::cout << "Compression Ratio : " << compression_ratio << std::endl;
        std::cout << "Maximum rank across the Tree : " << MAX_RANK << std::endl;
        std::cout << "Time to Initialize (in s) : " << HODLRdD_tree->get_Init_time() << std::endl;
        std::cout << "Matrix-Vector time (in s) : " << HODLRdD_tree->get_mat_vec_time() << std::endl;
        std::cout << "Time to solution   (in s) : " << this->solve_time << std::endl;
    }
    // Destructor
    ~HODLRdD_matrix(){
    }
};
#endif