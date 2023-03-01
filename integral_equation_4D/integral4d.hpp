#ifndef Integral2D_HPP
#define Integral2D_HPP
#include <iostream>
#include <cmath>
#include "../include/myHeaders.hpp"
#include "Gauss_Legendre_Nodes_and_Weights.hpp"

int n_gauss_nodes = numPoints;


double Integrand(double x,double y, double z, double w)
{
    return 1/(x*x + y*y + z*z + w*w);
}

// a = [x_low,y_low] and b = [x_upper,y_upper]
double quadruple_integral(double*& a, double*& b){
    double result = 0.0;
    double tx, ty, tz, tw, cx, cy, cz, cw, Lx, Ly, Lz, Lw;
    double *nodes_x,*weights_x;
    Gauss_Legendre_Nodes_and_Weights(n_gauss_nodes, nodes_x, weights_x);
    double *nodes_y,*weights_y;
    Gauss_Legendre_Nodes_and_Weights(n_gauss_nodes, nodes_y, weights_y);
    double *nodes_z,*weights_z;
    Gauss_Legendre_Nodes_and_Weights(n_gauss_nodes, nodes_z, weights_z);
    double *nodes_w,*weights_w;
    Gauss_Legendre_Nodes_and_Weights(n_gauss_nodes, nodes_w, weights_w);
    // Shift the nodes from [-1,1] to that interval
    cx = 0.5 * (b[0] - a[0]);
    cy = 0.5 * (b[1] - a[1]);
    cz = 0.5 * (b[2] - a[2]);
    cw = 0.5 * (b[3] - a[3]);

    Lx = 0.5 * (b[0] + a[0]);   // half length
    Ly = 0.5 * (b[1] + a[1]);
    Lz = 0.5 * (b[2] + a[2]);
    Lw = 0.5 * (b[3] + a[3]);

    // Gauss-Legendre Quadrature formula
    for(int i=0; i<n_gauss_nodes; i++)
    {
        tx = cx*nodes_x[i] + Lx;
        for(int j=0; j<n_gauss_nodes; j++)
        {
            ty = cy*nodes_y[j] + Ly;
            for(int k=0; k<n_gauss_nodes; k++)
            {
                tz = cz*nodes_z[k] + Lz;
                for(int l=0; l<n_gauss_nodes; l++)
                {
                    tw = cw*nodes_w[l] + Lw;
                    result += weights_x[i] * weights_y[j]* weights_z[k] * weights_w[l] * Integrand(tx,ty,tz,tw);
                }
            }
        }
    }
    // Final result due to change of variables
    result *= (cx*cy*cz*cw);
    //std::cout << "The quadruple integral is " << result << std::endl;
    return result;
}
#endif