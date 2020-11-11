#include <cassert>
#include <vector>
#include <iostream>
#include <ctime>

#include "matrix.h"

using namespace matrix;

namespace matrix_tests {

void unit_test0() {
    Square_Matrix<double> m({{2.1, 0.5},
                             {1.8, 5.4}});

    m.transpose();
    Square_Matrix<double> m_new({{2.1, 1.8},
                                 {0.5, 5.4}});

    assert(m_new == m);
    std::cout << "Unit_test0 complete.\n";
}

void unit_test1() {
    Square_Matrix<double> m({});

    m.transpose();
    Square_Matrix<double> m_new({});

    assert(m_new == m);
    std::cout << "Unit_test1 complete.\n";
}

void unit_test2() {
    Square_Matrix<int> m({{2, 0},
                          {1, 5}});

    assert(m(0, 0) == 2);
    std::cout << "Unit_test2 complete.\n";
}

void unit_test3() {
    Square_Matrix<double> m({{1.0, 0.0},
                             {0.0, 1.0}});

    double det = m.determinant();

    assert(det == 1.0);
    std::cout << "Unit_test3 complete.\n";
}

void unit_test4() {
    Square_Matrix<int> m({{1, 1},
                          {1, 1}});

    int det = m.determinant();

    assert(det == 0);
    std::cout << "Unit_test4 complete.\n";
}

void unit_test5() {
    Square_Matrix<int> m1({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m2({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m3({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m1_res({{3, 3},
                               {3, 3}});

    Square_Matrix<int> m2_res({{2, 2},
                               {2, 2}});

    Square_Matrix<int> m3_res({{1, 1},
                               {1, 1}});

    m1 += m2 += m3;

    assert(m3 == m3_res);
    assert(m2 == m2_res);
    assert(m1 == m1_res);
    std::cout << "Unit_test5 complete.\n";
}

void unit_test6() {
    Square_Matrix<int> m1({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m2({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m3({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m1_res({{-1, -1},
                               {-1, -1}});

    Square_Matrix<int> m2_res({{2, 2},
                               {2, 2}});

    Square_Matrix<int> m3_res({{1, 1},
                               {1, 1}});

    m1 -= m2 *= m3;

    assert(m3 == m3_res);
    assert(m2 == m2_res);
    assert(m1 == m1_res);
    std::cout << "Unit_test6 complete.\n";
}

void unit_test7() {
    Square_Matrix<int> m1({{1, 0},
                           {0, 1}});

    Square_Matrix<int> m2({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m3 = m1 * m2;

    Square_Matrix<int> m1_res({{1, 0},
                               {0, 1}});

    Square_Matrix<int> m2_res({{1, 1},
                               {1, 1}});

    Square_Matrix<int> m3_res({{1, 1},
                               {1, 1}});

    assert(m3 == m3_res);
    assert(m2 == m2_res);
    assert(m1 == m1_res);
    std::cout << "Unit_test7 complete.\n";
}

void unit_test8() {
    Square_Matrix<int> m({{1, 1, 1},
                          {1, 1, 1},
                          {1, 1, 1}});

    int det = m.determinant();

    assert(det == 0);
    std::cout << "Unit_test8 complete.\n";
}

void unit_test9() {
    Square_Matrix<int> m({{2, 0, 1},
                          {1, 0, 1},
                          {1, 1, 1}});

    int det = m.determinant();

    assert(det == -1);
    std::cout << "Unit_test9 complete.\n";
}

void unit_test10() {
    Square_Matrix<int> m({{0, 0, 1, 1},
                          {0, 0, 1, 2},
                          {1, 1, 1, 1},
                          {1, 2, 1, 1}});

    int det = m.determinant();

    assert(det == 1);
    std::cout << "Unit_test10 complete.\n";
}

namespace  {

Square_Matrix<long long int> make_add_matrix(const Square_Matrix<long long int>& matrix,
                                             size_t e_row,
                                             size_t e_column) {
    assert(e_row < matrix.size());
    assert(e_column < matrix.size());
    std::vector<std::vector<long long int>> new_rows;
    new_rows.resize(matrix.size() - 1);
    for(size_t i = 0; i < (matrix.size() - 1); ++i) {
        new_rows[i].resize(matrix.size() - 1);
    }

    size_t i_new = 0, j_new = 0;
    for(size_t i = 0; i < matrix.size(); ++i) {
        if(i != e_row) {
            j_new = 0;
            for(size_t j = 0; j < matrix.size(); ++j) {
                if(j != e_column) {
                    new_rows[i_new][j_new] = matrix(i, j);
                    ++j_new;
                }
            }
            ++i_new;
        }
    }

    return Square_Matrix<long long int>(new_rows);
}

//naive determinant
//write with templates too lazily
long long int bad_determinant(const Square_Matrix<long long int>& matrix) {
    if(matrix.size() == 0) {
        return 0;
    }
    if(matrix.size() == 1) {
        return matrix(0, 0);
    }

    long long int det = 0;
    for(size_t i = 0; i < matrix.size(); ++i) {
        Square_Matrix<long long int> add_matrix = make_add_matrix(matrix, 0, i);
        long long int minor = bad_determinant(add_matrix);

        if(i % 2 == 1) {
            minor *= -1;
        }

        det += (minor * matrix(0, i));
    }

    return det;
}

Square_Matrix<long long int> matrix_generator(size_t matrix_size) {
    std::vector<std::vector<long long int>> rows;
    rows.resize(matrix_size);
    for(size_t i = 0; i < (matrix_size); ++i) {
        rows[i].resize(matrix_size);
    }

    for(size_t i = 0; i < matrix_size; ++i) {
        for(size_t j = 0; j < matrix_size; ++j) {
            rows[i][j] = (rand() % 10);
        }
    }

    return Square_Matrix<long long int>(rows);
}
} //ending namespace

void long_test(int n_tests = 1000) {
    srand(0);
    size_t sum_bad_time = 0, sum_norm_time = 0;

    for(int i = 0; i < n_tests; ++i) {
        Square_Matrix<long long int>m = matrix_generator(7);

        long start_time = clock();
        long long int bad_det = bad_determinant(m);
        long bad_time = clock() - start_time;

        start_time = clock();
        long long int det = m.determinant();
        long norm_time = clock() - start_time;

        assert(det == bad_det);
        std::cout << "Bad det = " << bad_det << std::endl;
        std::cout << "Det = " << det << std::endl;
        std::cout << "Bad time = " << bad_time << std::endl;
        std::cout << "Normal time = " << norm_time << std::endl;

        sum_bad_time += static_cast<size_t>(bad_time);
        sum_norm_time += static_cast<size_t>(norm_time);
    }

    std::cout << "Sum bad time = " << sum_bad_time << std::endl;
    std::cout << "Sum normal time = " << sum_norm_time << std::endl;
    std::cout << "Long test complete\n";
}

void unit_test11() {
    Matrix<int> m1({{1, 1},
                    {1, 2},
                    {1, 1}});

    Matrix<int> m2({{1, 1, 1},
                    {1, 1, 1}});

    m1 = m1 * m2;

    Matrix<int> m1_new({{2, 2, 2},
                        {3, 3, 3},
                        {2, 2, 2}});
    assert(m1 == m1_new);

    std::cout << "Unit_test11 complete.\n";
}

void unit_test12() {
    Matrix<int> m1({{1, 1},
                    {1, 2},
                    {1, 1}});

    m1.transpose();

    Matrix<int> m1_new({{1, 1, 1},
                        {1, 2, 1}});
    assert(m1 == m1_new);

    std::cout << "Unit_test12 complete.\n";
}

}
