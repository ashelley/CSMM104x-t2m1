#ifndef DEBUG_HELPERS_H
#define DEBUG_HELPERS_H

#include <iostream>
#include <Eigen/Core>

#define DEBUG_MODE 1


static Eigen::IOFormat MATRIX_CLEANFORMAT(4, 0, ", ", "\n", "[", "]");

static inline 
void DEBUGPrintVector(VectorXs& v) {
    #if DEBUG_MODE
    std::cout << v.format(MATRIX_CLEANFORMAT) << std::endl;
    #endif
}

static inline
void DEBUGPrintVector(const VectorXs& v) {
    #if DEBUG_MODE
    std::cout << v.format(MATRIX_CLEANFORMAT) << std::endl;
    #endif
}

#endif