#ifndef __HDD_clusters__
#define __HDD_clusters__

#include "myHeaders.hpp"
#include "points_dt.hpp"

// Uniform clustering of points provided the bounding box, results in a set of clusters with N_max leaf points.
class cluster
{
    size_t N;
    std::vector<ptsnD> *gridPoints;
    double diam;
    double L;
    int level; // This is the level in the hierarchical clustering
    // Bounding box of the cluster! [x1[1], x2[1]] is the interval in dimension 1 and x1[1] < x2[1]
    // the bounding box here is provided by the user, TODO : If not provided the minimum possible hypercube that fits set of points needs to be done
    // The center of the cluster is not determined now, this we can do it while forming the Node of the tree.
public:
    size_t cluster_id;
    std::vector<size_t> index_of_points;
    double x1[NDIM], x2[NDIM];
    cluster(Eigen::VectorXd& x1_, Eigen::VectorXd& x2_,  std::vector<ptsnD>*& gPoints)
    {
        N = 0;
        this->gridPoints = gPoints;
        diam = 0.0;
        L = abs(x1_(0) - x2_(0));
        // TODO : Optimize copy
        for (int i = 0; i < NDIM; i++)
        {
            this->x1[i] = x1_(i);
            this->x2[i] = x2_(i);
            if (abs(this->x1[i] - this->x2[i]) > L)
                L = abs(this->x1[i] - this->x2[i]);
            diam += (x1_(i) - x2_(i)) * (x1_(i) - x2_(i));
            //std::cout << "[" << this->x1[i] << "," << this->x2[i] << "] "; 
        }
        //std::cout << std::endl;
        diam = sqrt(diam);
        cluster_id = 0;
    }
    void print_bounding_box(){
        std::cout << std::endl << "Bounding box : " ;
        for (int i = 0; i < NDIM; i++)
        {
            if(i != NDIM-1)
                std::cout << "[" << this->x1[i] << "," << this->x2[i] << "]x"; 
            else
                std::cout << "[" << this->x1[i] << "," << this->x2[i] << "]";
        }
        std::cout << std::endl;
    }
    void add_point(size_t a)
    {
        index_of_points.push_back(a);
        N++;
    }
    void add_points(std::vector<size_t> a)
    {
        index_of_points.insert(index_of_points.end(), a.begin(), a.end());
        N = index_of_points.size();
    }
    size_t get_cluster_size()
    {
        return N;
    }
    double get_diameter(){
        return diam;
    }
    void compute_cluster_center(ptsnD& c){
        // Compute cluster center
        //std::cout << "cluster center.."<<std::endl;
        for (int i = 0; i < NDIM; i++){
            //std::cout << "c.." << this->x1[i] << std::endl;
            //std::cout << "c.." << this->x2[i] << std::endl;
            c.x[i] = 0.5 * (this->x1[i] + this->x2[i]);
        }
        //std::cout << "cluster center computed.."<<std::endl;
    }
    void level_clustering(std::vector<cluster *>*& h_clusters){
        h_clusters->push_back(this);
        // Form 2^NDIM empty cluster by initialising through the bounding box of the parent
        for (int i = 0; i < NDIM; i++)
        {
            binary_clustering(h_clusters, i);
            //std::cout << i << " " << h_clusters->size() << std::endl;
        }
    }
    void binary_clustering(std::vector<cluster *>*& binary_clusters, int dim_i)
    {
        size_t n = binary_clusters->size();
        for (size_t i = 0; i < n; i++)
        {
            // Calculate boundary bdy
            double bdy = 0.0;
            double a = binary_clusters->at(0)->x1[dim_i];
            double b = binary_clusters->at(0)->x2[dim_i];
            bdy += 0.5 * (a + b);
            Eigen::VectorXd x1 = Eigen::VectorXd::Zero(NDIM);
            Eigen::VectorXd x2 = Eigen::VectorXd::Zero(NDIM);
            // Create the bounding box [x_1,x_2]^NDIM for the cluster
            for (int k = 0; k < NDIM; k++)
            {
                x1(k) = binary_clusters->at(0)->x1[k];
                x2(k) = binary_clusters->at(0)->x2[k];
            }
            // Create two clusters!! this subdivide the particles two ways TODO : HPC way<>
            x2(dim_i) = bdy;
            cluster *A, *B;
            A = new cluster(x1, x2, gridPoints);
            // Assigning cluster id = parent*2 in binary subdivision
            A->cluster_id = binary_clusters->at(0)->cluster_id * 2;
            x1(dim_i) = bdy;
            x2(dim_i) = b;
            B = new cluster(x1, x2, gridPoints);
            // Assigning cluster id = parent*2+1 in binary subdivision (this holds being the binary division of the domain)
            B->cluster_id = binary_clusters->at(0)->cluster_id * 2 + 1;

            // Here goes a while loop to subdivide the points in one dimension at the same time remove from source This ensure O(N) space for points at all times...
            for (size_t k = 0; k < binary_clusters->at(0)->get_cluster_size(); k++)
            {
                int point_i = binary_clusters->at(0)->index_of_points[k];
                // TODO : Use toggle to place the points inside the cluster A and B if it lies in the boundary
                if (gridPoints->at(point_i).x[dim_i] < bdy)
                    A->add_point(point_i);
                else
                    B->add_point(point_i);
            }
            binary_clusters->erase(binary_clusters->begin());
            binary_clusters->push_back(A);
            binary_clusters->push_back(B);
        }
    }
    double box_length(){
        return L;
    }
    friend int interaction_type(cluster*& A, cluster*& B);
    ~cluster(){
    }
};

int interaction_type(cluster*& A, cluster*& B){
    // int type_sharing = NDIM;
    // for(int i=0; i<NDIM; i++)
    //     if (abs(A->x2[i] - B->x2[i]) != 0 && abs(A->x1[i] - B->x1[i]) != 0) 
    //         type_sharing--;
    // return type_sharing;
    int type_sharing = NDIM;
    ptsnD a,b;
    A->compute_cluster_center(a);
    B->compute_cluster_center(b);
    double L = std::max(A->box_length(),B->box_length());
    // This ensures whether cluster |A|  |B| in the ith dimension image
    double dp = nd_points::euclidean_distance(a, b) / L;
    dp *= dp;
    //std::cout << "Frac : "<< dp <<std::endl;
    if (double(NDIM)+1 > dp) 
        type_sharing -= std::round(dp);
    else
        type_sharing = -1;
    return type_sharing;
}


// TODO : Find the bounding box
// TODO : Adaptive clustering at different levels
// TODO : Sort the points inside the cluster
#endif

/* Test case - Print points, binary clustering tested
            A->print_cluster();
            B->print_cluster();
    // Print point
    void print_cluster()
    {
        std::cout << "Cluster id - "<<cluster_id << std::endl;

        std::cout << "Bounding box " << std::endl;
        for (int j = 0; j < NDIM; j++)
            std::cout << "[" << x1[j]
                      << "," << x2[j] << "] ";
        std::cout << std::endl;
        for (int j = 0; j < NDIM; j++) std::cout << "x[" << j << "]"
                      << "     ";
        std::cout << N <<std::endl;
        for (int i = 0; i < N; i++)
        {
            std::cout << index_of_points[i] << "    ";
            for (int j = 0; j < NDIM; j++)
                std::cout << gridPoints->at(index_of_points[i]).x[j] << "     ";
            std::cout << std::endl;
        }
    }

*/