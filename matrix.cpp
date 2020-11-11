#include <vector>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include "matrix.h"

namespace matrix {

namespace {

const double DOUBLE_GAP = 1e-12;

bool is_match(double a, double b) { //instead of ==
    return (fabs(a - b) < DOUBLE_GAP);
}

//it is support function only for determinant()
template <typename T> T sum_for_determinant(std::vector<T>& L_str,
                                            std::vector<T>& U_str,
                                            size_t n) {
    assert(n <= L_str.size());
    assert(n <= U_str.size());
    assert(n >= 0);

    T sum = static_cast<T>(0);
    for(size_t i = 0; i < n; ++i) {
        sum += L_str[i] * U_str[i];
    }
    return sum;
}

//it is support function only for determinant()
template <typename T> void L_step(const Square_Matrix<T>& matrix,
                                  std::vector<std::vector<T>>& L,
                                  std::vector<std::vector<T>>& U,
                                  std::vector<T>& main_diag,
                                  size_t i) {
    L[i - 1][0] = matrix(i, 0) / main_diag[0];

    for(size_t j = 1; j < i; ++j) {
        L[i - 1][j] = (matrix(i, j) - sum_for_determinant<T>(L[i - 1], U[j - 1], j)) / main_diag[j];
    }
}

//it is support function only for determinant()
template <typename T> void U_step(const Square_Matrix<T>& matrix,
                                  std::vector<std::vector<T>>& L,
                                  std::vector<std::vector<T>>& U,
                                  size_t i) {

    for(size_t j = i; j < (matrix.size() - 1); ++j) {
        U[j][i] = matrix(i, j + 1) - sum_for_determinant<T>(L[i - 1], U[j], i);
    }
}

template<typename T>
 T universal_determinant(const Square_Matrix<T> input_matrix)  {
    if(input_matrix.size() == 0) {
        return static_cast<T>(0);
    }

    Square_Matrix<T> matrix(input_matrix);
    const size_t matr_size = matrix.size();
    std::vector<std::vector<T>> L;
    std::vector<std::vector<T>> U; //for LU decomposition it is transposed U (for speed)
    std::vector<T> main_diag; //of U

    //L and transposed U store only elems under their main diagonals
    //main_diag - it's U main diagonal
    //L main diagonal should have only "1" elements and not stores
    L.resize(matr_size - 1);
    for(size_t i = 0; i < (matr_size - 1); ++i) {
        L[i].resize(i + 1);
    }
    U.resize(matr_size - 1);
    for(size_t i = 0; i < (matr_size - 1); ++i) {
        U[i].resize(i + 1);
    }
    main_diag.resize(matr_size);


    //It isn't standart LU decompose method, because it should work correctly with some main minor = 0
    //
    //L and U are counted line by line and synchronously by steps
    //At each step we count (in this order):
    //1) 1 ROW of L
    //2) 1 diagonal elem
    //3) CHECK
    //4) 1 COLUMN of TRANSPOSED U (1 ROW of U)
    //
    //CHECK: if we have 0 in U main diag, there are two cases:
    //
    //1. we can fix it by adding other matrix rows under current row
    //   trying to make main minor nondegen
    //
    //   In this case we sholud repeat points 1), 2) and 3) after each attempt
    //
    //2. if (1) way can't fix main minor so the matrix is degen
    //   determinant = 0
    //
    //First step should be different with others (for convenience)
    //point 1) isn't exists there (but it isn't all difference)
    size_t str_n = 1;
    while(is_match(matrix(0, 0), static_cast<T>(0))) {
        if(str_n == matr_size) {
            return static_cast<T>(0);
        }
        matrix.add_row_to_row(str_n, 0);
        ++str_n;
    }

    main_diag[0] = matrix(0, 0);
    for(size_t j = 0; j < matr_size - 1; ++j) {
        U[j][0] = matrix(0, j + 1);
    }


    //other steps after the first
    for(size_t i = 1; i < matr_size; ++i) {
        size_t str_n;

        //building L
        L_step<T>(matrix, L, U, main_diag, i);

        //building main_diag, checking two special cases
        main_diag[i] = matrix(i, i) - sum_for_determinant<T>(L[i - 1], U[i - 1], i);

        //check
        str_n = i + 1;
        while(is_match(main_diag[i], static_cast<T>(0))) {
            if(str_n == matr_size) { //unfixing case
                return static_cast<T>(0);
            }

            //fix matrix
            matrix.add_row_to_row(str_n, i);

            //rebuild L row and main_diag elem
            L_step<T>(matrix, L, U, main_diag, i);
            main_diag[i] = matrix(i, i) - sum_for_determinant<T>(L[i - 1], U[i - 1], i);

            ++str_n;
        }


        //building U
        U_step<T>(matrix, L, U, i);

    }

    //calculating determinant
    T det = main_diag[0];
    for(size_t i = 1; i < matr_size; ++i) {
        det *= main_diag[i];
    }
    return det;
}

} //ending namespace

template<>
double Square_Matrix<double>::determinant() const {
    return universal_determinant(*this);
}

template<>
int Square_Matrix<int>::determinant() const {
    Square_Matrix<double> matr_buf(*this);
    return lround(universal_determinant(matr_buf));
}

template<>
long long int Square_Matrix<long long int>::determinant() const {
    Square_Matrix<double> matr_buf(*this);
    return llround(universal_determinant(matr_buf));
}

//... other types if will be need

} //ending namespace matrix {}
