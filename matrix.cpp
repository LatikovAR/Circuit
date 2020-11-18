#include <vector>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <utility>
#include <list>

#include "matrix.h"

namespace matrix {

namespace {

template <typename T>
class Universal_Determinant {
private:
    Square_Matrix<T> matrix_;
    const size_t size_;
    std::vector<std::vector<T>> L_;
    std::vector<std::vector<T>> U_; //for LU decomposition it is transposed U (for speed)
    std::vector<T> main_diag_;
    T determinant = static_cast<T>(0);

    static constexpr double DOUBLE_GAP = 1e-12;

    static bool is_match(double a, double b) { //instead of ==
        return (fabs(a - b) < DOUBLE_GAP);
    }

    T sum_for_determinant(size_t L_str_it, size_t U_str_it, size_t n);

    void L_step(size_t i);

    void U_step(size_t i);

public:
    Universal_Determinant(const Square_Matrix<T> input_matrix);

    T result() const { return determinant; }
};


template<typename T>
class Linear_Equations_System {
private:
    Matrix<T> matrix_;
    std::vector<T> column_;

    static constexpr double DOUBLE_GAP = 1e-12;

    static bool is_match(double a, double b) { //instead of ==
        return (fabs(a - b) < DOUBLE_GAP);
    }
public:
    Linear_Equations_System(Matrix<T> matrix, std::vector<T> column):
        matrix_(matrix), column_(column) {
        assert((matrix.column_size() == column.size()) && "invalid data sizes");
    }

    //returns any of the solutions if it isn't one
    //if no solution returning result undefined and bool = false
    //class object condition undefined after this method
    std::pair<std::vector<T>, bool> solve();


};

} // namespace

//-------------------------------------Determinant methods------------------------------------

template<>
double Square_Matrix<double>::determinant() const {
    return Universal_Determinant<double>(*this).result();
}

template<>
int Square_Matrix<int>::determinant() const {
    Square_Matrix<double> matr_buf(*this);
    return lround(Universal_Determinant<double>(matr_buf).result());
}

template<>
long long int Square_Matrix<long long int>::determinant() const {
    Square_Matrix<double> matr_buf(*this);
    return llround(Universal_Determinant<double>(matr_buf).result());
}

//... other types if will be need

namespace {

template <typename T>
T Universal_Determinant<T>::sum_for_determinant(size_t L_str_it, size_t U_str_it, size_t n) {
    assert(n <= L_[L_str_it].size());
    assert(n <= U_[U_str_it].size());
    assert(n >= 0);

    T sum = static_cast<T>(0);
    for(size_t i = 0; i < n; ++i) {
        sum += L_[L_str_it][i] * U_[U_str_it][i];
    }
    return sum;
}


template <typename T>
void Universal_Determinant<T>::L_step(size_t i) {
    L_[i - 1][0] = matrix_(i, 0) / main_diag_[0];

    for(size_t j = 1; j < i; ++j) {
        L_[i - 1][j] = (matrix_(i, j) - sum_for_determinant(i - 1, j - 1, j)) / main_diag_[j];
    }
}


template <typename T>
void Universal_Determinant<T>::U_step(size_t i) {

    for(size_t j = i; j < (matrix_.size() - 1); ++j) {
        U_[j][i] = matrix_(i, j + 1) - sum_for_determinant(i - 1, j, i);
    }
}

template<typename T>
Universal_Determinant<T>::Universal_Determinant(const Square_Matrix<T> input_matrix):
    matrix_(input_matrix),
    size_(matrix_.size())
{
    if(input_matrix.size() == 0) {
        return;
    }

    //L and transposed U store only elems under their main diagonals
    //main_diag - it's U main diagonal
    //L main diagonal should have only "1" elements and not stores
    L_.resize(size_ - 1);
    for(size_t i = 0; i < (size_ - 1); ++i) {
        L_[i].resize(i + 1);
    }
    U_.resize(size_ - 1);
    for(size_t i = 0; i < (size_ - 1); ++i) {
        U_[i].resize(i + 1);
    }
    main_diag_.resize(size_);


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
    while(is_match(matrix_(0, 0), static_cast<T>(0))) {
        if(str_n == size_) {
            return;
        }
        matrix_.add_row_to_row(str_n, 0);
        ++str_n;
    }

    main_diag_[0] = matrix_(0, 0);
    for(size_t j = 0; j < size_ - 1; ++j) {
        U_[j][0] = matrix_(0, j + 1);
    }


    //other steps after the first
    for(size_t i = 1; i < size_; ++i) {
        size_t str_n;

        //building L
        L_step(i);

        //building main_diag, checking two special cases
        main_diag_[i] = matrix_(i, i) - sum_for_determinant(i - 1, i - 1, i);

        //check
        str_n = i + 1;
        while(is_match(main_diag_[i], static_cast<T>(0))) {
            if(str_n == size_) { //unfixing case
                return;
            }

            //fix matrix
            matrix_.add_row_to_row(str_n, i);

            //rebuild L row and main_diag elem
            L_step(i);
            main_diag_[i] = matrix_(i, i) - sum_for_determinant(i - 1, i - 1, i);

            ++str_n;
        }


        //building U
        U_step(i);

    }

    //calculating determinant
    determinant = main_diag_[0];
    for(size_t i = 1; i < size_; ++i) {
        determinant *= main_diag_[i];
    }
}

} //namespace

//-------------------------------------Linear_Equations_System methods---------------------------

template<typename T1, typename T2>
std::vector<T2> vector_conversion(const std::vector<T1>& inp_vector) {
    std::vector<T2> res_vector;
    res_vector.resize(inp_vector.size());
    for(size_t i = 0; i < inp_vector.size(); ++i) {
        res_vector[i] = static_cast<T2>(inp_vector[i]);
    }
    return res_vector;
}

std::pair<std::vector<double>, bool> solve_linear_equations(Square_Matrix<double> matrix,
                                                            std::vector<double> column) {
    return Linear_Equations_System<double>(Matrix<double>(matrix), column).solve();
}

std::pair<std::vector<double>, bool> solve_linear_equations(Square_Matrix<int> matrix,
                                                            std::vector<int> column) {
    return Linear_Equations_System<double>(Matrix<double>(Square_Matrix<double>(matrix)),
                                           vector_conversion<int, double>(column)).solve();
}

std::pair<std::vector<double>, bool> solve_linear_equations(Square_Matrix<long long> matrix,
                                                            std::vector<long long> column) {
    return Linear_Equations_System<double>(Matrix<double>(Square_Matrix<double>(matrix)),
                                           vector_conversion<long long, double>(column)).solve();
}

std::pair<std::vector<double>, bool> solve_linear_equations(Matrix<double> matrix,
                                                            std::vector<double> column) {
    return Linear_Equations_System<double>(matrix, column).solve();
}

std::pair<std::vector<double>, bool> solve_linear_equations(Matrix<int> matrix,
                                                            std::vector<int> column) {
    return Linear_Equations_System<double>(Matrix<double>(matrix),
                                           vector_conversion<int, double>(column)).solve();
}

std::pair<std::vector<double>, bool> solve_linear_equations(Matrix<long long> matrix,
                                                            std::vector<long long> column) {
    return Linear_Equations_System<double>(Matrix<double>(matrix),
                                           vector_conversion<long long, double>(column)).solve();
}

namespace {

//uses Gauss method
//solution will get with two passes:
//
//1) straight pass - building an up triangle matrix with "1" on the main diag
//
//2) reverse pass - building an identity matrix with solution in the boung column
template<typename T>
std::pair<std::vector<T>, bool> Linear_Equations_System<T>::solve() {
    std::vector<T> answer;
    answer.resize(matrix_.row_size());
    std::list<bool> trace; //for going back the same way as straight pass

    //first pass (straight)

    // number of already calculated rows/columns
    size_t row_built = 0;
    size_t column_built = 0;

    while((row_built < matrix_.column_size()) && (column_built < matrix_.row_size())) {

        size_t suitable_row_num = row_built; //with row[column built] != 0
        while((suitable_row_num < matrix_.column_size()) &&
              is_match(matrix_(suitable_row_num, column_built), 0)) {
            ++suitable_row_num;
        }

        if(suitable_row_num == matrix_.column_size()) {
            ++column_built;
            trace.push_back(false);
            continue;
        }

        if(suitable_row_num != row_built) {
            matrix_.swap_rows(row_built, suitable_row_num);
            std::swap(column_[row_built], column_[suitable_row_num]);
        }


        T norm_koef = matrix_(row_built, column_built);
        assert(!is_match(norm_koef, 0));
        norm_koef = 1 / norm_koef;
        matrix_.mult_row_on_number(row_built, norm_koef);
        column_[row_built] *= norm_koef;

        for(size_t i = row_built + 1; i < matrix_.column_size(); ++i) {
            T k = matrix_(i, column_built);
            matrix_.sub_row_to_row(row_built, i, k);
            column_[i] -= (column_[row_built] * k);
        }

        ++column_built;
        ++row_built;
        trace.push_back(true);
    }

    //down rows (0, 0 ..... 0 | ?)
    //if ? != 0 - no solution
    for(size_t i = row_built; i < matrix_.column_size(); ++i) {
        if(!is_match(column_[i], 0)) {
            return std::pair<std::vector<T>, bool>(answer, false);
        }
    }

    //second pass (reverse)

    --column_built;
    --row_built;
    size_t i = matrix_.row_size();
    while(i > 0) {
        --i;
        if(i > column_built) {
            answer[i] = 0.0;
            continue;
        }

        if(trace.back() == true) {
            for(size_t j = 0; j < row_built; ++j) {
                T k = matrix_(j, column_built);
                matrix_.sub_row_to_row(row_built, j, k);
                column_[j] -= (column_[row_built] * k);
            }
            answer[i] = column_[row_built];
            --column_built;
            --row_built;
            trace.pop_back();
        }

        else {
            answer[i] = 0.0;
            --column_built;
            trace.pop_back();
        }
    }

    return std::pair<std::vector<T>, bool>(answer, true);
}

} //namespace

} //ending namespace matrix
