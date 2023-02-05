#ifndef Hclustering_HPP
#define Hclustering_HPP


#include "myHeaders.hpp"

const int NDIM = 2;

// Struture that holds the individual location of points
struct ptsnD{
    double x[NDIM];
    size_t id;  // User provides a unique identifier if not then kernel is defined by considering two points
};

// Eucledian distance between two points in any dimensions
double compute_distance(ptsnD& a,ptsnD& b){
    double sum = 0.0;
    for(int i=0;i< NDIM;i++){
        sum += (a.x[i] - b.x[i]) * (a.x[i] - b.x[i]);
    }
    // TODO : Update the above to reduce operation in openMP
    sum = sqrt(sum);
    return sum;
}

class cluster{
    size_t N;
    // Bounding box of the cluster! [x1[1], x2[1]] is the interval in dimension 1 and x1[1] < x2[1]
    // the bounding box here is provided by the user, TODO : If not provided the minimum possible hypercube that fits set of points needs to be done
    // The center of the cluster is not determined now, this we can do it while forming the Node of the tree.
public:
    double x1[NDIM], x2[NDIM];
    std::vector<ptsnD> gridPoints;
    size_t cluster_id;
    int level; // This is the level in the hierarchical clustering
    cluster(Eigen::VectorXd x1, Eigen::VectorXd x2)
    {
        N = 0;
        // TODO : Optimize copy
        for (int i = 0; i < NDIM; i++)
        {
            this->x1[i] = x1(i);
            this->x2[i] = x2(i);
        }
        }
        void add_point(ptsnD a){
            gridPoints.push_back(a);
            N++;
        }
        void add_points(std::vector<ptsnD>& gridPoints_){
            gridPoints.insert(gridPoints.end(), gridPoints_.begin(), gridPoints_.end());
            N = gridPoints.size();
        }
        void remove_point(size_t i){
            gridPoints.erase(gridPoints.begin() + i);
            N--;
        }

        void remove_all_points(){
            gridPoints.clear();
            N = 0;
        } 

        size_t get_cluster_size(){
            return N;
        }

        ~cluster(){
        }
};

// Uniform clustering of points provided the bounding box, results in a set of clusters with N_max leaf points. TODO : Find the bounding box
class Hclustering{
    int nLevel = 0 ;
    bool isLeaf = false;
    size_t N_max;
    size_t nClusters;
    bool isAdaptive = false; // TODO: Needs work Adaptive clustering
    // this class takes a complete set of points and provides 
    std::vector<cluster> h_clusters;
    public:
        double x1[NDIM], x2[NDIM]; // Bounding box of the root and is set at the initialization of the cluster
        Hclustering(cluster& src, size_t Nmx){
            N_max = Nmx;
            //  Updating the root domain in the clustering
            src.cluster_id = 0;
            src.level = 0;
            // Bounding box is set
            for(int i=0;i<NDIM; i++){
            this->x1[i] = src.x1[i];
            this->x2[i] = src.x2[i];
            }
            h_clusters.push_back(src); 
            // sanity check on whether cluster size is greater than N_max specified
            if(N_max < src.get_cluster_size())
                getHClustering();
        }
        // Call to this routine the first time 
        void getHClustering(){         
            while (!isLeaf){
                double N_buff = 0;
                nLevel++;
                size_t cluster_count = h_clusters.size();
                for (size_t l = 0; l < cluster_count; l++){
                    std::vector<cluster> binary_clusters;
                    binary_clusters.push_back(h_clusters[0]);
                    // Form 2^NDIM empty cluster by initialising through the bounding box of the parent
                    for(int i=0;i <NDIM;i++){
                        binary_clustering(binary_clusters, i);
                    }
                    // Remove the first in the cluster (this works like a queue, work added at back)
                    h_clusters.erase(h_clusters.begin());
                    // Adding the subdivided clusters
                    h_clusters.insert(h_clusters.end(), binary_clusters.begin(), binary_clusters.end());
                }
                // Check the top and add each and every point to the  appropriate clusters.
                for (size_t cl = 0; cl < h_clusters.size(); cl++)
                {
                    if (N_buff < h_clusters[cl].get_cluster_size())
                        N_buff = h_clusters[cl].get_cluster_size();
                }
                if (N_buff < N_max)
                    isLeaf = true;
            }
        }
        void binary_clustering(std::vector<cluster>& binary_clusters,int dim_i){
            size_t n = binary_clusters.size();
            for(size_t i=0;i<n;i++){
                // Calculate boundary bdy
                double bdy = 0.0;
                double a = binary_clusters[0].x1[dim_i];
                double b = binary_clusters[0].x2[dim_i];
                bdy += 0.5 * (a+b);
                Eigen::VectorXd x1 = Eigen::VectorXd::Zero(NDIM);
                Eigen::VectorXd x2 = Eigen::VectorXd::Zero(NDIM);
                // Create the bounding box [x_1,x_2]^NDIM for the cluster
                for(int k =0 ;k<NDIM;k++){
                    x1(k) = binary_clusters[0].x1[k];
                    x2(k) = binary_clusters[0].x2[k];
                }
                // Create two clusters!! this subdivide the particles two ways TODO : HPC way<>
                x2(dim_i) = bdy;
                cluster A(x1,x2);
                // Assigning cluster id = parent*2 in binary subdivision
                A.cluster_id  = binary_clusters[0].cluster_id*2; 
                x1(dim_i) = bdy;
                x2(dim_i) = b;
                cluster B(x1,x2);
                // Assigning cluster id = parent*2+1 in binary subdivision (this holds being the binary division of the domain)
                A.cluster_id = binary_clusters[0].cluster_id * 2 + 1;

                // Here goes a while loop to subdivide the points in one dimension at the same time remove from source This ensure O(N) space for points at all times...
                while (binary_clusters[0].gridPoints.size()!=0){
                    // TODO : Use toggle to place the points inside the cluster A and B if it lies in the boundary
                    if (binary_clusters[0].gridPoints[0].x[dim_i]<bdy){
                        A.add_point(binary_clusters[0].gridPoints[0]);
                        binary_clusters[0].remove_point(0);
                    }
                    else{
                        B.add_point(binary_clusters[0].gridPoints[0]);
                        binary_clusters[0].remove_point(0);
                    }
                }
                binary_clusters.erase(binary_clusters.begin());
            }
        }
        ~Hclustering(){
        }
};

// TODO : Adaptive clustering at different levels
// TODO : Sort the points inside the cluster
#endif

