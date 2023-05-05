#ifndef __HDD_node__
#define __HDD_node__

#include "myHeaders.hpp"
#include "points_dt.hpp"
#include "HDD_clusters.hpp"
#include "kernel_function.hpp"
#include "LowRank_matrix.hpp"

template <class Kernel>
class Node
{
    kernel_function<Kernel> *userkernel;
    size_t self_id;

    // Matrix operators ... Full memory
    Mat P2P_self;
    Mat *P2P, *P2M, *L2P;
    Vec node_charge;
    Vec node_potential;
    // Matrix operators ... reduced memory
    // TODO : MEM efficient ACA
    bool is_admissible(Node<Kernel>*& A){
        // Check the max norm between the center of two clusters
        double dist_btwn_cluster_center = nd_points::max_norm_distance(this->cluster_center, A->cluster_center);
        if (dist_btwn_cluster_center < 1.5 * this->my_cluster->box_length())
            if (interaction_type(this->my_cluster, A->my_cluster) > INTERACTION_TYPE_ALLOWED)
                return false;
        return true;
    }

public:
    Node<Kernel> *parent;
    cluster *my_cluster;
    ptsnD cluster_center; // Computed post initialisation of all the nodes for the tree
    size_t n_particles;
    bool isleaf = false;
    bool isroot;
    std::vector<Node<Kernel> *> Child;
    std::vector<Node<Kernel> *> my_neighbour_addr;
    int n_neighbours, n_intraction;
    std::vector<Node<Kernel> *> my_intr_list_addr;
    double my_flop_il = 0.0;
    size_t node_rank;
    // Node can be initialised with just parent and child information
    Node(cluster*& src, kernel_function<Kernel>*& usr_)
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
    Node(cluster*& src, Node*& parent_, kernel_function<Kernel>*& usr_)
    {
        this->my_cluster = src;
        this->parent = parent_;
        this->userkernel = usr_;
        n_neighbours = 0;
        n_intraction = 0;
        n_particles = src->get_cluster_size();
        //std::cout << "cluster size " << n_particles << std::endl;
        // Update child list of the parents
        //parent->Child.push_back(this);
        isroot = false;
        self_id = src->cluster_id;
        my_cluster->compute_cluster_center(cluster_center);
        cluster_center.id = -1; // center points provided with id -1
    }
    size_t get_id(){
        return this->self_id;
    }

    void Initialize_node();
    void find_max_rank(size_t &MAX_RANK);
    void update_matrix_node(std::vector<std::vector<int>>& mat_color, int lvl);
    double compute_flop_count();
    void get_interaction_list();
    void get_node_potential();
    void set_node_charge(const Vec &b){
        node_charge = Vec::Zero(n_particles);
        node_potential = Vec::Zero(n_particles);
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
    void print_node_details(){
        std::cout << "**************" << std::endl;
        std::cout << "Node [" << self_id << "] , Nunmber particles = " << n_particles << std::endl;
        my_cluster->print_bounding_box();
        std::cout << "Neighbours -";
        for (int i = 0; i < n_neighbours; i++)
            std::cout << " " << my_neighbour_addr[i]->self_id;
        std::cout << std::endl;
        std::cout << "Interaction List -";
        for (int i = 0; i < n_intraction; i++)
            std::cout << " " << this->my_intr_list_addr[i]->self_id;
        std::cout << std::endl;
        std::cout << "**************" << std::endl;
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
        // for (int i = 0; i < n_intraction; i++)
        //     if (my_intr_list_addr[i]->n_particles != 0)
        //         LR[i] = LowRankMat<Kernel>(userkernel, my_cluster->index_of_points,
        //                                        my_intr_list_addr[i]->my_cluster->index_of_points);
        // for (int i = 0; i < n_intraction; i++) 
        //     if (my_intr_list_addr[i]->n_particles != 0)
        //         userkernel->ACA_FAST(L2P[i], P2M[i], eps_ACA,
        //                                          my_cluster->index_of_points,
        //                                          my_intr_list_addr[i]->my_cluster->index_of_points);
        // if (isleaf)
        // {
        //     // P2P = new Mat[n_neighbours];
        //     // //P2P_self = userkernel->getMatrix(my_cluster->index_of_points, my_cluster->index_of_points);
        //     // for (int i = 0; i < n_neighbours; i++)
        //     //     if (my_neighbour_addr[i]->n_particles != 0)
        //     //         P2P[i] = userkernel->getMatrix(my_cluster->index_of_points,
        //     //                                        my_neighbour_addr[i]->my_cluster->index_of_points);
        // }
    }
    node_potential = Vec::Zero(n_particles);
    //TODO : Initialise the matrix operator - Reduced memory
}
template<class Kernel>
void Node<Kernel>::get_interaction_list()
{
    // std::cout << "Self--- " << self_id << std::endl;
    if(!isroot){
        // The interaction list consists of nodes from two sources
        //      1) Siblings
        for (size_t i = 0; i < this->parent->Child.size(); i++){
            Node<Kernel>* tmp;
            tmp = this->parent->Child[i];
            if (this->self_id != tmp->self_id){
                if(this->is_admissible(tmp))
                    this->my_intr_list_addr.push_back(tmp);
                else
                    this->my_neighbour_addr.push_back(tmp);
            }
        }
        //      2) Children of one's parent's neighbours
        for (size_t i = 0; i < this->parent->my_neighbour_addr.size(); i++){
            for (size_t j = 0; j < this->parent->my_neighbour_addr[i]->Child.size(); j++){
                Node<Kernel> *tmp;
                tmp = this->parent->my_neighbour_addr[i]->Child[j];
                if (this->is_admissible(tmp))
                    this->my_intr_list_addr.push_back(tmp);
                else
                    this->my_neighbour_addr.push_back(tmp);
            }   
        }
    }
    this->n_neighbours = this->my_neighbour_addr.size();
    this->n_intraction = this->my_intr_list_addr.size();
}

    template <class Kernel>
    double Node<Kernel>::compute_flop_count()
    {
        double my_flop_nn = 0.0;
        if(n_particles != 0){
            if(this->isleaf){
                my_flop_nn += n_particles;
                for (int i = 0; i < n_neighbours; i++){
                    if (my_neighbour_addr[i]->n_particles != 0)
                        my_flop_nn += my_neighbour_addr[i]->n_particles;
                }
                my_flop_nn *= n_particles;
            }
            // for (int i = 0; i < n_intraction; i++){
            //     if (my_intr_list_addr[i]->n_particles != 0){
            //         my_flop += LR[i].rank() * (n_particles + my_intr_list_addr[i]->n_particles);
            //     }
            // }
        }    
        return my_flop_nn;
    }

    // template <class Kernel>
    // void Node<Kernel>::find_max_rank(size_t& MAX_RANK){
    //     size_t tmp = MAX_RANK;
    //     if(n_particles != 0){
    //         for (int i = 0; i < n_intraction; i++)
    //         {
    //             if (my_intr_list_addr[i]->n_particles != 0)
    //             {
    //                 if ((size_t) LR[i].rank() > tmp)
    //                 {
    //                     tmp = LR[i].rank();
    //                 }
    //             }
    //         }
    //     }
    //     MAX_RANK = tmp;
    // }

    template <class Kernel>
    void Node<Kernel>::get_node_potential()
    {
    // Routine for Full memory matvec
    //std::cout << "Self ID " << self_id << std::endl;
    if(isleaf)
    {
        if(n_particles != 0)
        {
            //std::cout << "(" << P2P_self.rows() << "," << P2P_self.cols() <<") x ("<< node_charge.size() << "x1)" << std::endl;
            //node_potential += P2P_self * node_charge;
            for (size_t i = 0; i < n_particles; i++)
                    node_potential(i) += userkernel->getRow(my_cluster->index_of_points[i], my_cluster->index_of_points).dot(node_charge);
            //std::cout << "Neighbor id ";
            for (int i = 0; i < n_neighbours; i++)
            {
                //std::cout  << my_neighbour_addr[i]->self_id << " " << std::endl;
                if (my_neighbour_addr[i]->n_particles != 0)
                    //node_potential += P2P[i] * my_neighbour_addr[i]->node_charge;
                    for (size_t j = 0; j < n_particles; j++)
                        node_potential(j) += userkernel->getRow(my_cluster->index_of_points[j], my_neighbour_addr[i]->my_cluster->index_of_points).dot(my_neighbour_addr[i]->node_charge);
            }
            //std::cout << "Neighbor id ";
        }
    }
    //std::cout << std::endl;
    //std::cout << "Interaction  id " ;
    for (int i = 0; i < n_intraction; i++)
    {
        //std::cout << my_intr_list_addr[i]->self_id << " ";
        // if (my_intr_list_addr[i]->n_particles != 0)
        //     node_potential += (L2P[i] * (P2M[i].transpose() * my_intr_list_addr[i]->node_charge));
        if (my_intr_list_addr[i]->n_particles != 0){
            double tmp = omp_get_wtime();
            LowRankMat<Kernel> *LR = new LowRankMat<Kernel>[n_intraction];
            LR[i] = LowRankMat<Kernel>(userkernel, my_cluster->index_of_points,
                                               my_intr_list_addr[i]->my_cluster->index_of_points);
            init_time += omp_get_wtime() - tmp;
            node_potential += LR[i] * my_intr_list_addr[i]->node_charge;
            tmp = omp_get_wtime();
            this->my_flop_il += LR[i].rank() * (n_particles + my_intr_list_addr[i]->n_particles);
            if (this->node_rank < LR[i].rank())
                this->node_rank = LR[i].rank();
            meas_time += omp_get_wtime() - tmp;
        }
    }
    }
    //std::cout << std::endl;
    template<class Kernel>
    void Node<Kernel>::update_matrix_node(std::vector<std::vector<int> >& mat_color, int lvl)
    {
        // Create source list corresponding to leaf nodes, lvl is 4^{leaf_level - my_level}
        std::vector<size_t> row_indices;
        int drift = pow(pow(2,NDIM), lvl);
        int start_index = drift * int(this->self_id);
        for(int i=0; i< drift; i++)
            row_indices.push_back((size_t)start_index + i);

        // update the mat color based on the interaction list in the matrix
        for (int i = 0; i < n_intraction; i++)
        {
            // Create target list corresponding to leaf nodes
            std::vector<size_t> col_indices;      

            start_index = drift * int(my_intr_list_addr[i]->self_id);
            for (int k = 0; k < drift; k++)
                col_indices.push_back((size_t)start_index + k);

            // Make appropriate color changes
            for (int j = 0; j < (int)row_indices.size(); j++)
                for (int k = 0; k < (int)col_indices.size(); k++)
                {
                    mat_color[row_indices[j]][col_indices[k]] = 1;
                }
        }
    }

    // TODO : Routine Reduced memory matvec
#endif