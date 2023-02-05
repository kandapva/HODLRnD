#include<iostream>
#include <eigen3/Eigen/Dense>

#include "HDD_clusters.hpp"

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

int main()
{
    VectorXd Xdir, Ydir;
    std::vector<ptsnD> *gridPoints; // location of particles in the domain

    std::cout << "Cheb Points" << std::endl;
    Eigen::VectorXd X(2), Y(2);
    X(0) = -1;
    X(1) = -1;
    Y(0) = 1;
    Y(1) = 1;
    gridPoints = new std::vector<ptsnD>; // location of particles in the domain
    int cts = 0;
    Xdir = cheb_nodes(X(0), Y(0), 6);
    std::cout << Xdir << std::endl;
    Ydir = cheb_nodes(X(1), Y(1), 6);
    for (size_t i = 0; i < 5; i++)
        for (size_t j = 0; j < 5; j++)
        {
            ptsnD temp;
            temp.x[0] = Ydir[i];
            temp.x[1] = Xdir[j];
            temp.id = cts++;
            gridPoints->push_back(temp);
        }
    std::vector<int> v(25);                  // vector with 100 ints.
    std::iota(std::begin(v), std::end(v), 0); // Fill with 0, 1, ..., 99.
    cluster A(X,Y,gridPoints);
    A.add_points(v);
    std::vector<cluster *> t;
    A.print_cluster();
    A.level_clustering(t);
    //for(int i=0;i<t.size();i++)
        //t[i]->print_cluster();
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