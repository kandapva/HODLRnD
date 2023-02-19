#ifndef __HDD_node__
#define __HDD_node__

#include "myHeaders.hpp"
#include "points_dt.hpp"
#include "HDD_clusters.hpp"
#include "kernel_function.hpp"

template <class Kernel>
class Node
{
    kernel_function<Kernel> *userkernel;
    size_t self_id;
    bool isleaf = false;
    bool isroot;
    // Matrix operators ... Full memory
    Mat P2P_self;
    Mat *P2P, *P2M, *L2P;
    Vec node_charge;
    Vec node_potential;
    // Matrix operators ... reduced memory
    // TODO : MEM efficient ACA
    bool is_admissible(Node *A){
        // Check the max norm between the center of two clusters
        double dist_btwn_cluster_center = max_norm_distance(this->cluster_center, A->cluster_center);
        if (dist_btwn_cluster_center < 1.5 * this->my_cluster->get_diameter())
            if (interaction_type(this, A) > INTERACTION_TYPE_ALLOWED)
                return false;
        return true;
    }

public:
    Node *parent;
    cluster *my_cluster;
    ptsnD cluster_center; // Computed post initialisation of all the nodes for the tree
    size_t n_particles;
    std::vector<Node *> Child;
    std::vector<Node *> my_neighbour_addr;
    int n_neighbours, n_intraction;
    std::vector<Node *> my_intr_list_addr;
    // Node can be initialised with just parent and child information
    Node(cluste *&src, kernel_function<Kernel>*& usr_)
    {
        this->my_cluster = src;
        this->userkernel = usr_;
        n_neighbours = 0;
        n_intraction = 0;
        n_particles = src->get_cluster_size();
        isroot = true;
        self_id = 0;
        this->my_cluster->compute_cluster_center(cluster_center);
        cluster_center.id = -1;
    }
    Node(cluster *&src, Node *&parent_, kernel_function<Kernel>*& usr_)
    {
        this->parent = parent_;
        this->userkernel = usr_;
        n_neighbours = 0;
        n_intraction = 0;
        n_particles = src->get_cluster_size();
        // Update child list of the parents
        parent->Child.push_back(this);
        isroot = false;
        self_id = src->cluster_id;
        this->my_cluster->compute_cluster_center(cluster_center);
        cluster_center.id = -1; // center points provided with id -1
    }
    void Initialize_node();
    void get_interaction_list();
    void get_node_potential();
    void set_node_charge(const Vec &b){
        node_charge = Vec::Zero(n_particles);
        for (size_t i = 0; i < n_particles; i++)
        {
            size_t tmp = my_cluster->index_of_points[i];
            node_charge(i) = b(tmp);
        }
    }
    void collect_potential(Vec& b){
        // Collects the potential
        // node_charge = Vec::Zero(n_particles);
        for (size_t i = 0; i < n_particles; i++)
        {
            size_t tmp = my_cluster->index_of_points[i];
            b(tmp) += node_potential(i);
        }
    }
};

// Routine to compute the interaction list and neighbor list
template<class Kernel>
void Node<Kernel>::Initialize_node()
{
    // Initialise the matrix operators - Full memory
    if (n_particles != 0)
    {
        node_charge = Vec::Zero(n_particles);
        L2P = new Mat[n_intraction];
        P2M = new Mat[n_intraction];
        for (int i = 0; i < n_intraction; i++)
            if (my_intr_list_addr[i]->n_particles != 0)
                userkernel->ACA_FAST(L2P[i], P2M[i], eps_ACA,
                                     my_cluster->index_of_points,
                                     my_intr_list_addr[i]->my_cluster->index_of_points);
        if (isleaf)
        {
            P2P = new Mat[n_neighbours];
            P2P_self = userkernel->getMatrix(my_cluster->index_of_points, my_cluster->index_of_points);
            for (int i = 0; i < n_neighbours; i++)
                if (my_neighbour_addr[i]->n_particles != 0)
                    P2P[i] = userkernel->getMatrix(my_cluster->index_of_points,
                                                   my_neighbour_addr[i]->my_cluster->index_of_points);
        }
    }
    //TODO : Initialise the matrix operator - Reduced memory
}
template<class Kernel>
void Node<Kernel>::get_interaction_list()
{
    if(!isroot){
        // The interaction list consists of nodes from two sources
        //      1) Siblings
        for (int i = 0; i < this->parent->Child.size(); i++){
            Node<Kernel> *tmp = this->parent->Child[i];
            if (this->self_id != tmp->self_id)
                if(this->is_admissible(tmp))
                    this->my_intr_list_addr.append(tmp);
                else
                    this->my_neighbour_addr.append(tmp);
        }
        //      2) Children of one's parent's neighbours
        for (int i = 0; i < this->parent->my_neighbour_addr.size(); i++){
            for (int j = 0; j < this->parent->my_neighbour_addr[i]->Child.size(); j++){
                Node<Kernel> *tmp = this->parent->my_neighbour_addr[i]->Child[j];
                if (this->is_admissible(tmp))
                    this->my_intr_list_addr.append(tmp);
                else
                    this->my_neighbour_addr.append(tmp);
            }   
        }
    }
    this->n_neighbours = this->my_neighbour_addr.size();
    this->n_intraction = this->my_intr_list_addr.size();
}

template<class Kernel>
void Node<Kernel>::get_node_potential(){
    // Routine for Full memory matvec
    if(isleaf){
        if(n_particles != 0){
            node_potential += P2P_self * node_charge;
            for (int i = 0; i < n_neighbours; i++)
                if (my_neighbour_addr[i]->n_particles != 0)
                    node_potential += P2P[i] * my_neighbour_addr[i]->node_charge;
        }
    }
    for (int i = 0; i < n_intraction; i++)
        if (my_intr_list_addr[i]->n_particles != 0)
            node_potential += (L2P[i] * (P2M[i].transpose() * my_neighbour_addr[i]->node_charge));
    // TODO : Routine Reduced memory matvec
}
#endif