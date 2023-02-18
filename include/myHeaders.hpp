#ifndef myHeaders_HPP
#define myHeaders_HPP

#include <iostream>
#include <chrono>
#include <ctime>
#include <set>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <complex>
#include <iterator>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <list>
#include <algorithm>
#include <numeric>
#include <stack>

#define _USE_MATH_DEFINES
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using std::string;
using namespace Eigen;

// HODLR matrix parameters
const int NDIM = 2;
const int Nmax = 500;
// The admissibility is based on the max norm of the center
const double eta_admissible = 1;
const double eps_aca = pow(10,-8);

int mod(int a, int b){
    return ((a % b + b) % b);
}

#ifdef USE_FLOAT
    using dtype = float;
    using dtype_base = float;
    using Mat = Eigen::MatrixXf;
    using Vec = Eigen::VectorXf;
    #define Calc_dist(x1, y1, x2, y2) float(sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)))
#endif

#ifdef USE_DOUBLE
    using dtype = double;
    using dtype_base = double;
    using Mat = Eigen::MatrixXd;
    using Vec = Eigen::VectorXd;
    #define abs_(x) ((x < 0) ? (-x) : (x))
    #define conj_(x) (x)
    #define Calc_dist(x1, y1, x2, y2) double(sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)))
#endif

#ifdef USE_COMPLEX32
    using dtype = std::complex<float>;
    using dtype_base = float;
    using Mat = Eigen::MatrixXcf;
    using Vec = Eigen::VectorXcf;
    const std::complex<float> I(0.0, 1.0);
#endif

#ifdef USE_COMPLEX64
    using dtype = std::complex<double>;
    using dtype_base = double;
    using Mat = Eigen::MatrixXcd;
    using Vec = Eigen::VectorXcd;
    #define abs_(x) (std::abs((x)))
    #define conj_(x) (std::conj((x)))
    const std::complex<double> I(0.0, 1.0);
    #define Calc_dist(x1, y1, x2, y2) float(sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)))
#endif

std::string inline getTimeStamp()
{
    std::string s;
    s = "\n" + string(__DATE__) + "," + string(__TIME__) + "\n";
    return s;
}

namespace Vec_ops
{
    dtype_base inline relative_error(Vec &X_ori, Vec &X_comp)
    {
        dtype_base result = 0.0;
        result = (X_ori - X_comp).norm() / X_ori.norm();
        return result;
    }
    dtype inline dot_product(const Vec &x, const Vec &y)
    {
        dtype tmp = dtype(0.0);
        for (size_t i = 0; i < x.size(); i++)
            tmp += x(i) * y(i);
        return tmp;
    }
}

#ifdef USE_CEREAL
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
namespace storedata
{
    Eigen::VectorXd inline load_vec(std::string fname)
    {
        Eigen::VectorXd X;
        std::vector<dtype_base> K;
        std::ifstream infile(fname, std::ios::binary);
        cereal::PortableBinaryInputArchive iarchive(infile); // Create an input archive
        iarchive(K);                                         // Read the data from the archive
        X = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(K.data(), K.size());
        return X;
    }

    void inline save_vec(std::string fname, Eigen::VectorXd X)
    {
        std::vector<dtype_base> K(X.data(), X.data() + X.size());
        std::ofstream outfile(fname, std::ios::binary);
        cereal::PortableBinaryOutputArchive oarchive(outfile); // Create an output archive
        oarchive(CEREAL_NVP(K));                               // Write the data to the archive
    }
}
#endif

#endif

