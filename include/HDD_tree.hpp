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
    std::vector<Node<Kernel> *> obj_lin_arr;
    Node<Kernel> *root;
    double MAT_VEC_TIME = 0.0;
    double INIT_TIME = 0.0;
    int level;
    cluster *src;
public:
// Constructor for the matrix data structure
    Tree(cluster*& src_, std::vector<ptsnD>*& gPoints, kernel_function<Kernel>*& usr_)
    {
        this->gridPoints = gPoints;
        this->src = src_;
        int Nlevel = gridPoints->size();
        //std::cout << Nlevel << std::endl;
        std::vector<size_t> v(Nlevel);
        for (int i = 0; i < Nlevel; i++){
            v[i] = i;
        }
        level = 0;
        //    std::iota(std::begin(v), std::end(v), size_t(0)); 
        src->add_points(v);
        // create a root node
        root = new Node<Kernel>(src, usr_);
        obj_arr.resize(1);
        obj_arr[level].push_back(root);
        level++;
        std::cout << "Root Formed" << std::endl;
        // Perform clustering and form the tree
        while (Nlevel > Nmax)
        {
            obj_arr.resize(level + 1);
            // Form the next level through Hierarchical clustering 
            for(size_t i=0; i<obj_arr[level-1].size(); i++){
                std::vector<cluster *> *t;
                t = new std::vector<cluster *>;
                obj_arr[level-1][i]->my_cluster->level_clustering(t);
                //std::cout << "Cluster recieved" << std::endl;
                for(size_t j=0; j<t->size();j++){
                    Node<Kernel> *temp;
                    temp = new Node<Kernel>(t->at(j), obj_arr[level-1][i], usr_);
                    obj_arr[level].push_back(temp);
                    obj_arr[level - 1][i]->Child.push_back(temp);
                    //std::cout<< "Node (" << level << "," << temp->get_id() << ") Formed" << std::endl;
                }
                delete t;
            }
            // Decide whether leaf level reached
            Nlevel = obj_arr[level][0]->n_particles; 
            for (size_t j = 0; j < obj_arr[level].size(); j++) {
                obj_arr[level][j]->get_interaction_list();
                if (obj_arr[level][j]->n_particles > (unsigned) Nlevel)
                    Nlevel = obj_arr[level][j]->n_particles;
            }
            level++;
        }
        std::cout << "Tree Formed with depth " << level-1 << std::endl;
        // Mark level as leaf
        for (size_t j = 0; j < obj_arr[level-1].size(); j++)
            obj_arr[level-1][j]->isleaf = true;
    }
    // Initialise the matrix operators
    void Initialise_tree(){
        double start = omp_get_wtime();
        for (int i = 0; i < level; i++)
            for (size_t j = 0; j < obj_arr[i].size(); j++)
                obj_arr[i][j]->Initialize_node();
        INIT_TIME = omp_get_wtime() - start;
        std::cout << "Matrix operators formed..." << std::endl;
    }
    
    // mat-vec product
    Vec mat_vec(Vec& x){
        Vec b = Vec::Zero(x.size());
        // Using openMP here requires a memeory race while reading
        for (int i = 0; i < level; i++)
            for (size_t j = 0; j < obj_arr[i].size(); j++)
                obj_arr[i][j]->set_node_charge(x);
        //std::cout << "x set" << std::endl;
        // perform node wise matrix vector
// #pragma omp parallel for collapse(2) schedule(dynamic)
        double start = omp_get_wtime();
        for (int i = 0; i < level; i++)
            for (size_t j = 0; j < obj_arr[i].size(); j++)
                obj_arr[i][j]->get_node_potential();
        MAT_VEC_TIME = omp_get_wtime() - start;
        // std::cout << "b compute" << std::endl;
        //  Collect the output vector
        for (int i = 0; i < level; i++) for (size_t j = 0; j < obj_arr[i].size(); j++)
            obj_arr[i][j]->collect_potential(b);
        // std::cout << "b calculate" << std::endl;
        // std::cout << "Mat-vec performed..." << std::endl;
        return b;
        }
        void get_stat(double &n_FLOP, size_t &MAX_RANK)
        {
            for (int i = 0; i < level; i++)
                for (size_t j = 0; j < obj_arr[i].size(); j++)
                {
                    n_FLOP += obj_arr[i][j]->compute_flop_count();
                    obj_arr[i][j]->find_max_rank(MAX_RANK);
                }
        }
        void print_tree_details()
        {
            std::cout << "Tree Depth:" << level << std::endl;
            std::cout << "______________________________" << std::endl;
            std::cout << std::endl;
            for (int i = 0; i < level; i++)
            {
                std::cout << "Level[" << i << "]" << std::endl;
                for (size_t j = 0; j < obj_arr[i].size(); j++)
                {
                    obj_arr[i][j]->print_node_details();
                    std::cout << std::endl;
                }
                std::cout << "______________________________" << std::endl;
            }
        }
        double get_mat_vec_time(){
            return MAT_VEC_TIME;
        }
        double get_Init_time()
        {
            return INIT_TIME;
        }
        // Destructor
        ~Tree()
        {
        }
    };
#endif