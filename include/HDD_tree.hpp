#ifndef __HDD_tree__
#define __HDD_tree__

#include "myHeaders.hpp"
#include "HDD_node.hpp" 
#include "HDD_clusters.hpp"

// Hierarchical Tree structure
template <class Kernel>
class Tree
{
    size_t N;
    std::vector<ptsnD> *gridPoints;
    std::vector<std::vector<Node<Kernel> *>> obj_arr;
    Node<Kernel> *root;
    int level = 0;
    cluster *src;
public:
// Constructor for the matrix data structure
    Tree(cluster*& src_, std::vector<ptsnD>*& gPoints, kernel_function<Kernel>*& usr_)
    {
        this->gridPoints = gPoints;
        this->src = src_;
        int Nlevel = gridPoints->size();
        std::vector<size_t> v(Nlevel);
        for (size_t i = 0; i < Nlevel; i++){
            v[i] = i;
        }
        //    std::iota(std::begin(v), std::end(v), size_t(0)); 
        src->add_points(v);
        // create a root node
        root = new Node(src, usr_);
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
                    Node<Kernel> *temp;
                    temp = new Node<Kernel>(t[j], obj_arr[level][i], usr_);
                    obj_arr[level+1].push_back(temp);
                }
            }
            // Decide whether leaf level reached
            Nlevel = obj_arr[level+1][0]->n_particles; 
            for (size_t j = 0; j < obj_arr[level + 1].size(); j++) {
                obj_arr[level + 1][j]->get_interaction_list();
                if (obj_arr[level + 1][j]->n_particles > Nlevel)
                    Nlevel = obj_arr[level + 1][j]->n_particles;
            }
                
            level++;
        }
        // Mark level as leaf
        for (size_t j = 0; j < obj_arr[level + 1].size(); j++)
            obj_arr[level][j]->isleaf = true;
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
        for (int i = 0; i <= level; i++)
            for (size_t j = 0; j < obj_arr[i].size(); j++)
                obj_arr[i][j]->set_node_charge(x);
        // perform node wise matrix vector
        for (int i = 0; i <= level; i++) 
            for (size_t j = 0; j < obj_arr[i].size(); j++)
                obj_arr[i][j]->get_node_potential();
        // Collect the output vector
        for (int i = 0; i <= level; i++)
            for (size_t j = 0; j < obj_arr[i].size(); j++)
                obj_arr[i][j]->collect_potential(b);
    }
    // Destructor
    ~Tree(){
    }
};
#endif