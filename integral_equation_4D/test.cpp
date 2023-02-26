#include "integral4d.hpp"
int main(){
    double* a = new double[4];
    double* b = new double[4];

    a[0] = -1, a[1] = -1, a[2] = -1, a[3] = -1;
    b[0] = 1, b[1] = 1, b[2] = 1, b[3] = 1;
    std::cout << quadruple_integral(a,b) << std::endl;
}