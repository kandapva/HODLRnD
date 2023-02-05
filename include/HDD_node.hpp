#ifndef __HDD_node__
#define __HDD_node__

#include "myHeaders.hpp"
#include "points_dt.hpp"
#include "HDD_clusters.hpp"


class Node{
    size_t self_id;
    ptsnD cluster_center; // Computed post initialisation of all the nodes for the tree
    bool isleaf = false;
    bool isroot;
    // Matrix operators ... Full memory
    Mat P2P_self;
    Mat *P2P, *P2M, *L2P;
    Vec node_charge;
    Vec node_potential;
    // Matrix operators ... reduced memory
    // TODO : MEM efficient ACA

public:
    Node *parent;
    cluster *my_cluster;
    size_t n_particles;
    std::vector<Node *> Child;
    std::vector<Node *> my_neighbour_addr;
    int n_neighbours, n_intraction;
    std::vector<Node *> my_intr_list_addr;
    // Node can be initialised with just parent and child information
    Node(cluster* src){
        this->my_cluster = src;
        n_neighbours = 0;
        n_intraction = 0;
        n_particles = src->get_cluster_size();
        isroot = true;
        self_id = 0;
        this->my_cluster->compute_cluster_center(cluster_center);
        cluster_center.id = -1;
    }
    Node(cluster* src, Node *&parent_)
    {
        this->parent = parent_;
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
    void set_node_charge(Vec &b);
    void collect_potential(Vec &b);
};

// Routine to compute the interaction list and neighbor list
void Node::Initialize_node()
{
    // Initialise the matrix operators - Full memory

    //TODO : Initialise the matrix operator - Reduced memory
}

void Node::get_interaction_list(){

    n_neighbours = my_neighbour_addr.size();
    n_intraction = my_intr_list_addr.size();
}
void Node::set_node_charge(Vec& b){
    node_charge = Vec::Zero(n_particles);
    for (size_t i = 0; i < n_particles; i++){
        size_t tmp = my_cluster->index_of_points[i];
        node_charge(i) = b(tmp);
    }
}
void Node::get_node_potential(){
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

void Node::collect_potential(Vec &b)
{
    // Collects the potential
    node_charge = Vec::Zero(n_particles);
    for (size_t i = 0; i < n_particles; i++)
    {
        size_t tmp = my_cluster->index_of_points[i];
        b(tmp) = node_charge(i);
    }
}
#endif