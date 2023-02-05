#ifndef hodlrnd_node_HPP
#define hodlrnd_node_HPP


// LookUp Table for generating the Neighbours and interaction list
// Left Looking C ordering
const static int North[4] = {3, 2, -1, -1};
const static int South[4] = {-1, -1, 1, 0};
const static int East[4] = {1, -1, -1, 2};
const static int West[4] = {-1, 0, 3, -1};
const static int North_im[4] = {-1, -1, 1, 0};
const static int South_im[4] = {3, 2, -1, -1};
const static int East_im[4] = {-1, 0, 3, -1};
const static int West_im[4] = {1, -1, -1, 2};
const static int LR_region[4] = {2, 3, 0, 1};
const static double tol = pow(10, -10);
size_t MAX_RANK_H2 = 0;

class Node
{
    int north, south, east, west;

public:
    Node *parent;
    // Interactions Data
    points_dt *charges;
    Kernel *K_reorder;
    LowRank *lr;

    long int wload;
    long int lload;
    double fload;
    Vec node_charge;
    std::vector<Node *> Child;
    std::vector<size_t> crank;
    std::vector<Node *> my_neighbour_addr;
    std::vector<Node *> my_intr_list_addr;
    std::vector<int> *idx = new std::vector<int>;
    std::vector<int> *holdMyData = new std::vector<int>;
    int myproc;
    int myprocd;

    unsigned level_num, self_num, self_block;
    double x_start, y_start, x_end, y_end;
    double cx, cy;

    int charge_start, charge_end;
    int n_particles, n_neighbour, n_intraction;
    bool isleaf = false;
    bool isroot;
    std::vector<int> DLR, AFLR; // Diagonal Low Rank, Adjacent Full Rank Block
    Mat P2P_self;
    Mat *P2P, *P2M, *L2P;

    Node(double xs, double xe, double ys, double ye, points_dt *&charges) // Parent Node
    {
        level_num = 0;
        self_num = 0;
        north = -1;
        south = -1;
        east = -1;
        west = -1;
        x_start = xs;
        x_end = xe;
        y_start = ys;
        y_end = ye;
        isroot = true;
        cx = 0.5 * (xe + xs);
        cy = 0.5 * (ye + ys);
        this->charges = charges;
        myproc = 0;
        myprocd = 0;
        n_neighbour = 0;
        n_intraction = 0;
        wload = 0;
        std::cout << "Tol - " << tol << std::endl;
    }

    Node(unsigned level_num, unsigned self_num, Node *&obj, points_dt *&charges)
    {
        parent = obj;
        this->charges = charges;
        isroot = false;
        this->level_num = level_num;
        this->self_num = self_num;
        self_block = self_num % 4;
        myproc = 0;  // self_block%2;
        myprocd = 0; // self_block%2;
        get_domain();
        get_news();
        get_Low_rank_list();
    }
    // Gets the domain limits for a Node
    void get_domain();

    // Gets the Adjacent/Edge shared blocks
    void get_news();

    // Collects the indices of node for Low rank factorization
    void get_Low_rank_list();

    // Initalizes the HODLR2D operators
    void Initialize_node();
    void mat_vec(const Vec &x);
    void particle_to_child();
    void mark_leaf();
    // Utilities
    void set_max_rank();
    void print_self(string fname);

    void myload();
    double myload_full();
    ~Node()
    {
        delete idx;
    }
};

void Node::mark_leaf()
{
    this->isleaf = true;
    for (int i = 0; i < holdMyData->size(); i++)
        this->idx->push_back(holdMyData->at(i));
    n_particles = this->idx->size();
    delete holdMyData;
}

void Node::get_domain()
{
    if (self_block == 0)
    {
        x_start = parent->x_start;
        x_end = parent->cx;
        y_start = parent->y_start;
        y_end = parent->cy;
    }
    if (self_block == 1)
    {
        x_start = parent->cx;
        x_end = parent->x_end;
        y_start = parent->y_start;
        y_end = parent->cy;
    }
    if (self_block == 2)
    {
        x_start = parent->cx;
        x_end = parent->x_end;
        y_start = parent->cy;
        y_end = parent->y_end;
    }
    if (self_block == 3)
    {
        x_start = parent->x_start;
        x_end = parent->cx;
        y_start = parent->cy;
        y_end = parent->y_end;
    }
    cx = 0.5 * (x_end + x_start);
    cy = 0.5 * (y_end + y_start);
}

void Node::get_news()
{
    if (North[self_block] != -1)
        north = 4 * parent->self_num + North[self_block];
    else
    {
        if (parent->north != -1)
            north = 4 * parent->north + North_im[self_block];
        else
            north = -1;
    }
    if (South[self_block] != -1)
        south = 4 * parent->self_num + South[self_block];
    else
    {
        if (parent->south != -1)
            south = 4 * parent->south + South_im[self_block];
        else
            south = -1;
    }
    if (East[self_block] != -1)
        east = 4 * parent->self_num + East[self_block];
    else
    {
        if (parent->east != -1)
            east = 4 * parent->east + East_im[self_block];
        else
            east = -1;
    }
    if (West[self_block] != -1)
        west = 4 * parent->self_num + West[self_block];
    else
    {
        if (parent->west != -1)
            west = 4 * parent->west + West_im[self_block];
        else
            west = -1;
    }
    // Create the Vector with SWEN
    if (south != -1)
        AFLR.push_back(south);
    if (west != -1)
        AFLR.push_back(west);
    if (east != -1)
        AFLR.push_back(east);
    if (north != -1)
        AFLR.push_back(north);
    n_neighbour = AFLR.size();
}

void Node::get_Low_rank_list()
{
    DLR.push_back(parent->self_num * 4 + LR_region[self_block]);
    if (parent->south != -1)
        for (int i = 0; i < 4; i++)
            if (south != parent->south * 4 + i)
                DLR.push_back(parent->south * 4 + i);
    if (parent->west != -1)
        for (int i = 0; i < 4; i++)
            if (west != parent->west * 4 + i)
                DLR.push_back(parent->west * 4 + i);
    if (parent->east != -1)
        for (int i = 0; i < 4; i++)
            if (east != parent->east * 4 + i)
                DLR.push_back(parent->east * 4 + i);
    if (parent->north != -1)
        for (int i = 0; i < 4; i++)
            if (north != parent->north * 4 + i)
                DLR.push_back(parent->north * 4 + i);
    n_intraction = DLR.size();
}

void Node::Initialize_node()
{
    if (n_particles != 0)
    {
        node_charge = Vec::Zero(n_particles);
        L2P = new Mat[n_intraction];
        P2M = new Mat[n_intraction];
        for (int i = 0; i < n_intraction; i++)
            if (my_intr_list_addr[i]->n_particles != 0)
                lr->getFactorization(L2P[i], P2M[i], tol, charge_start,
                                     my_intr_list_addr[i]->charge_start,
                                     n_particles,
                                     my_intr_list_addr[i]->n_particles);
        if (isleaf)
        {
            P2P = new Mat[n_neighbour];
            P2P_self = K_reorder->getMatrix(charge_start, charge_start, n_particles, n_particles);
            for (int i = 0; i < n_neighbour; i++)
                if (my_neighbour_addr[i]->n_particles != 0)
                    P2P[i] = K_reorder->getMatrix(charge_start,
                                                  my_neighbour_addr[i]->charge_start,
                                                  n_particles,
                                                  my_neighbour_addr[i]->n_particles);
        }
    }
}

void Node::set_max_rank()
{
    if (n_particles != 0)
        for (int i = 0; i < n_intraction; i++)
        {
            if (my_intr_list_addr[i]->n_particles != 0)
            {
                crank.push_back(L2P[i].cols());
                if (MAX_RANK_H2 < L2P[i].cols())
                    MAX_RANK_H2 = L2P[i].cols();
            }
        }
}
void Node::mat_vec(const Vec &x) //,Vec& b)
{
    if (n_particles != 0)
    {
        for (int i = 0; i < n_intraction; i++)
            if (my_intr_list_addr[i]->n_particles != 0)
                node_charge +=
                    (L2P[i] * (P2M[i].transpose() * x.segment(my_intr_list_addr[i]->charge_start, my_intr_list_addr[i]->n_particles)));
        if (isleaf)
        {
            node_charge += P2P_self * x.segment(charge_start, n_particles);
            for (int i = 0; i < n_neighbour; i++)
                if (my_neighbour_addr[i]->n_particles != 0)
                    node_charge += P2P[i] * x.segment(my_neighbour_addr[i]->charge_start, my_neighbour_addr[i]->n_particles);
        }
    }
}





#endif