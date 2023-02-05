#ifndef __HDD_tree__
#define __HDD_tree__

#include "myHeaders.hpp"
#include "HDD_node.hpp" 
#include "HDD_clusters.hpp"

// Hierarchical Tree structure 
class Tree{
size_t N;
std::vector<ptsnD> *gridPoints;
std::vector<std::vector<Node *>> obj_arr;
Node *root;
int level = 0;
cluster *src;
public:
// Constructor for the matrix data structure
    Tree(Eigen::VectorXd x1, Eigen::VectorXd x2, std::vector<ptsnD> *gPoints)
    {
        // Source luster needs special treatment i.e., provide a guide list that maps the points
        this->gridPoints = gPoints;
        src = new cluster(x1, x2, this->gridPoints);
        int Nlevel = gridPoints->size();
        std::vector<int> v(Nlevel);                  
        std::iota(std::begin(v), std::end(v), 0); 
        src->add_points(v);
        // create a root node
        root = new Node(src);
        obj_arr[level].push_back(root);

        // Perform clustering and form the tree
        while (Nlevel > Nmax)
        {
            obj_arr.resize(level + 1);
            // Form the next level through Hierarchical clustering 
            for(size_t i=0; i<obj_arr[level].size(); i++){
                std::vector<cluster *> t;
                obj_arr[level][i]->my_cluster->level_clustering(t);
                for(int j=0; j<t.size();j++){
                    Node *temp = new Node(t[j], obj_arr[level][i]);
                    obj_arr[level+1].push_back(temp);
                }
            }
            // Decide whether leaf level reached
            Nlevel = obj_arr[level+1][0]->n_particles; 
            for (size_t j = 0; j < obj_arr[level + 1].size(); j++) 
                if (obj_arr[level+1][j]->n_particles > Nlevel)
                    Nlevel = obj_arr[level + 1][j]->n_particles;
            level++;
        }
        // Mark level as leaf
    // Connection to the other nodes
    }
    // Initialise the matrix operators
    void Initialise_tree(){
        for(int i=0;i<=level;i++)
            for (size_t j = 0; j < obj_arr[i].size(); j++)
                obj_arr[i][j]->Initialize_node();
    }
    
    // mat-vec product
    Vec mat_vec(const Vec& x){
        Vec b = Vec::Zero(x.size());
        // Order the output vector 


        // perform node wise matrix vector 
        for(int i=0;i<=level;i++)
            for (size_t j = 0; j < obj_arr[i].size(); j++)
                obj_arr[i][j]->get_node_potential();
    }
    // Destructor
    ~Tree(){
    }
};
#endif