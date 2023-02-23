#include<iostream>
#include <Eigen/Dense>

#include "HDD_matrix.hpp"
#include "HDD_clusters.hpp"
#include "kernel_function.hpp"
#include "myHeaders.hpp"
#include "points_dt.hpp"
#include "LowRank_matrix.hpp"

using namespace std;

Eigen::VectorXd cheb_nodes(double a, double b, int n)
{
    Eigen::VectorXd X(n);
    double l, l1, param;
    l = 0.5 * (a + b);
    l1 = 0.5 * (b - a);
    for (int k = 0; k < n; k++)
    {
        param = (double)(k + 0.5) / n;
        X(k) = l - l1 * cos(param * 3.1412);
    }
    return X;
}

class userkernel{
    std::vector<ptsnD> gridPoints; // location of particles in the domain
public:
    userkernel()
    {
        VectorXd Xdir, Ydir;
        int numPoints = 10;
        // std::cout << "Cheb Points" << std::endl;
        Eigen::VectorXd X(2), Y(2);
        X(0) = -1;
        X(1) = -1;
        Y(0) = 1;
        Y(1) = 1;
        //gridPoints = new std::vector<ptsnD>; // location of particles in the domain
        int cts = 0;
        Xdir = cheb_nodes(X(0), Y(0), numPoints);
        //std::cout << Xdir << std::endl;
        Ydir = cheb_nodes(X(1), Y(1), numPoints);
        int N = numPoints * numPoints;
        for (size_t i = 0; i < numPoints; i++)
            for (size_t j = 0; j < numPoints; j++)
            {
                ptsnD temp;
                temp.x[0] = Ydir[i];
                temp.x[1] = Xdir[j];
                temp.id = cts++;
                gridPoints.push_back(temp);
            }
        }
    dtype_base getMatrixEntry(int i, int j){
        double r = 0.0;
        if(i != j){
        ptsnD a,b;
        a = gridPoints[i];
        b = gridPoints[j];
        r = nd_points::euclidean_distance(a, b);
        r = 1/r;
        }
        return r;
    }

    ~userkernel(){
    }
};


int main()
{
    VectorXd Xdir, Ydir;
    std::vector<ptsnD> *gridPoints; // location of particles in the domain
    int numPoints = 10;
    //std::cout << "Cheb Points" << std::endl;
    Eigen::VectorXd X(2), Y(2);
    X(0) = -1;
    X(1) = -1;
    Y(0) = 1;
    Y(1) = 1;
    gridPoints = new std::vector<ptsnD>; // location of particles in the domain
    int cts = 0;
    Xdir = cheb_nodes(X(0), Y(0), numPoints);
    //std::cout << Xdir << std::endl;
    Ydir = cheb_nodes(X(1), Y(1), numPoints);
    int N = numPoints * numPoints;
    for (size_t i = 0; i < numPoints; i++) 
        for (size_t j = 0; j < numPoints; j++)
        {
            ptsnD temp;
            temp.x[0] = Ydir[i];
            temp.x[1] = Xdir[j];
            temp.id = cts++;
            gridPoints->push_back(temp);
        }
    std::vector<size_t> v1;                  // vector with 100 ints.
    std::vector<size_t> v2;                  // vector with 100 ints.
    std::vector<size_t> v3;
    for(size_t i=0;i<50;i++)
        v1.push_back(i);
    for (size_t i = 50; i < 100; i++)
        v2.push_back(i);
    //std::iota(std::begin(v1), std::end(v1), 0); // Fill with 0, 1, ..., 99.
    
    //std::iota(std::begin(v2), std::end(v2), 50); // Fill with 0, 1, ..., 99.
    std::cout << "Grid_points A" << std::endl;
    cluster A(X,Y,gridPoints);
    A.add_points(v1);
    std::cout << "Grid_points B" << std::endl;
    cluster B(X, Y, gridPoints);
    B.add_points(v2);
    

    std::vector<cluster *> t;
    //A.print_cluster();
    // A.level_clustering(t);
    // for(int i=0;i<t.size();i++)
    //     t[i]->print_cluster();
    
    userkernel *add;
    add = new userkernel();
    kernel_function<userkernel> *kernelfunc = new kernel_function<userkernel>(add);
    HODLRdD_matrix Kmat = HODLRdD_matrix(add, gridPoints, X, Y);
    Kmat.Assemble_matrix_operators();
    //std::cout << "Entry " << add->getMatrixEntry(1,2) << std::endl;
    Vec x, b, bl;
    b = Vec::Zero(50);
    x = Vec::Ones(50, 1);
    Vec b1,b2,x_test;
    x_test = Vec::Ones(100, 1);
    std::cout << Kmat.get_size() << std::endl;
    b1 = Kmat * x_test;
    //std::cout << b1 << std::endl;
    //std::cout << Kmat.get_size() << std::endl;
    for (size_t i = 0; i < 100; i++)
        v3.push_back(i);
    b2 = kernelfunc->getMatrix(v3, v3) * x_test;
    //std::cout << b2 << std::endl;
    std::cout << "Relative Error.. hmatrix ... " << Vec_ops::relative_error(b2, b1) << std::endl;
    // std::cout << kernelfunc->getRow(1, v2) << endl;
    // std::cout << kernelfunc->getCol(1,v1) << endl;
    Mat K,L, R;

    K = kernelfunc->getMatrix(v1,v2);
    LowRankMat<userkernel> LR(kernelfunc, v1, v2, true);
    bl = LR * x;
    std::cout << K.rows() << "," << K.cols() << std::endl;
    std::cout << x.size() << std::endl;

    b = K*x;
    std::cout << "Rank of ACA.. " << LR.rank() << std::endl;
    std::cout << "Relative Error (ACA).." << Vec_ops::relative_error(b, bl) << std::endl;
    kernelfunc->ACA_FAST(L,R,0.000001,v1,v2);
    bl = L*(R.transpose()*x);
    std::cout << "Relative Error (ACA).."<< Vec_ops::relative_error(b,bl) << std::endl;
    return 0;
    }

/* points_dt Test case
******Headers******
#include "points_dt.hpp"
*******************
 // Verify the max_norm, 2-norm and the namespaces
    ptsnD A,B;
    A.x[0] = 0.5;
    A.x[1] = 0.5;
    B.x[0] = 0.75;
    B.x[1] = 0.75;

    cout << "Max norm .." << nd_points::max_norm_distance(A, B) << endl;
    cout << "2-norm .." << nd_points::euclidean_distance(A, B) << endl;

*/

    /*  Interaction type of two boxes provided the center of the clusters
        ptsnD c1,c2;
        c1.x[0] = 0.5;
        c1.x[1] = 0.5;
        c1.x[2] = 0.5;
        c2.x[0] = 2.5;
        c2.x[1] = 1.5;
        c2.x[2] = 1.5;
        int tmp = INTERACTION_type(c1, c2, 1.0);
        cout << "Interaction type .. "<< tmp << std::endl;

int INTERACTION_type(ptsnD A, ptsnD B, double L)
{
    int type_sharing = NDIM;
    // This ensures whether cluster |A|  |B| in the ith dimension image
    double dp =  nd_points::euclidean_distance(A, B) / L;
    dp *= dp;
    //std::cout << "Frac : "<< dp <<std::endl;
    if (double(NDIM)+1 > dp)
        type_sharing -= std::round(dp);
    else
        type_sharing = -1;
    return type_sharing;
}
    */