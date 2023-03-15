# HODLRnD
HODLR in higher dimensions. 

Edit the following hyperparameters based on your kernel and application. 
NDIM -> The dimension space of the particles
Nmax -> The maximum number of particles at the leaf level (This decides the depth of the hierarchical tree)
eps_ACA -> Accuracy in compression by ACA
SYS_SIZE -> Above this size the matrix operations are performed in an memory efficient way
INTERACTION_TYPE_ALLOWED -> This represents d' (to decide on which weak admissibility criteria)

The following options are allowed in deciding the hierarchical tree and/or the admissibility criteria
In Two dimenstions:
* d' = -1 -> H-matrix with standard admissibility criteria $\eta = \sqrt{2}$.
* d' = 0 -> Interaction list includes Vertex sharing neighbours, HODLR2D matrix structure.
* d' = 1 -> Interaction list includes all clusters that are not self. HODLR with quad tree 
In Three dimensions:
* d' = -1 -> H-matrix with standard admissibility criteria $\eta = \sqrt{3}$.
* d' = 0 -> Interaction list includes Vertex sharing neighbours, HODLR3D matrix structure.
* d' = 1 -> Interaction list includes Vertex sharing and edge sharing neighbours
* d' = 2 -> Interaction list includes all clusters that are not self. HODLR with oct tree  
In Higher dimensions:
* d' = -1 -> H-matrix with standard admissibility criteria $\eta = \sqrt{NDIM}$.
* d' = 0 -> Interaction list includes Vertex sharing neighbours, HODLRnD matrix structure.
* ....
* d' = NDIM-1 -> Interaction list includes all clusters that are not self. HODLR with $2^{NDIM}$ tree.  

## Guidelines for Pull request:
* While providing a pull request, please ensure that the executable does not gets included.
* One possible workaround would be to compile with extensions .app, .exe which is ignored. 
