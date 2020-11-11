#pragma once

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <utility>
#include <cmath>
#include <new>
#include <typeinfo>

namespace matrix {

template <typename T> class Abstract_Matrix {
protected:
    static constexpr double DOUBLE_GAP = 1e-12;

    static bool is_match(double lhs, double rhs) { //instead of ==
        return (fabs(lhs - rhs) < DOUBLE_GAP);
    }

    static bool is_match(int lhs, int rhs) {
        return (lhs == rhs);
    }

    static bool is_match(long long int lhs, long long int rhs) {
        return (lhs == rhs);
    }

public:
    virtual ~Abstract_Matrix() {} //why =0 don't works?
    virtual T operator()(size_t row_i, size_t column_i) const = 0;
    virtual bool operator==(const Abstract_Matrix<T>& rhs) const = 0;
    virtual void print() const = 0;
};



template <typename T>
class Square_Matrix final : public Abstract_Matrix<T> {
private:
    size_t size_ = 0;
    T *data_ = nullptr;

    using Abstract_Matrix<T>::is_match;

public:
    Square_Matrix(const T* inp_data, size_t inp_size);

    Square_Matrix(const std::vector<std::vector<T>>& input_rows);

    ~Square_Matrix() override;

    Square_Matrix(const Square_Matrix& matr): Square_Matrix(matr.data_, matr.size_) {}

    Square_Matrix& operator= (const Square_Matrix& matr)&;

    template<typename U> Square_Matrix(const Square_Matrix<U>& matr);

    void transpose() const;

    void add_row_to_row(size_t src_num, size_t dst_num) const;

    void sub_row_to_row(size_t src_num, size_t dst_num) const;

    T operator()(size_t row_i, size_t column_i) const override;

    size_t size() const { return size_; }

    bool operator==(const Abstract_Matrix<T>& inp_rhs) const override;

    Square_Matrix& operator+=(const Square_Matrix& rhs)&;

    Square_Matrix& operator-=(const Square_Matrix<T>& rhs)&;

    Square_Matrix& operator*=(const Square_Matrix<T>& rhs)&;

    void print() const override;

    T determinant() const;
};


template <typename T> class Matrix final : public Abstract_Matrix<T> {
private:
    size_t column_size_ = 0;
    size_t row_size_ = 0;
    T *data_ = nullptr;

    using Abstract_Matrix<T>::is_match;

public:
    //matrix elems should be contained in inp_data row by row
    Matrix(const T* inp_data, size_t inp_column_size, size_t inp_row_size);

    Matrix(const std::vector<std::vector<T>>& input_rows);

    ~Matrix() override;

    Matrix(const Matrix& matr): Matrix(matr.data_, matr.column_size_, matr.row_size_) {}

    Matrix& operator= (const Matrix& matr)&;

    template<typename U> Matrix(const Matrix<U>& matr);

    void transpose();

    void add_row_to_row(size_t src_num, size_t dst_num) const;

    void sub_row_to_row(size_t src_num, size_t dst_num) const;

    T operator()(size_t row_i, size_t column_i) const override;

    size_t row_size() const { return row_size_; }
    size_t column_size() const { return column_size_; }

    bool operator==(const Abstract_Matrix<T>& inp_rhs) const override;

    Matrix& operator+=(const Matrix& rhs)&;

    Matrix& operator-=(const Matrix<T>& rhs)&;

    Matrix& operator*=(const Matrix<T>& rhs)&;

    void print() const override;
};


//----------------------------methods for Square_Matrix--------------------------------
template<typename T>
Square_Matrix<T>::Square_Matrix(const T* inp_data, size_t inp_size):
    size_(inp_size)
{
    if(size_ > 0) data_ = new T[size_ * size_];

    for(size_t i = 0; i < size_ * size_; ++i) {
        data_[i] = inp_data[i];
    }
}

template<typename T>
Square_Matrix<T>::Square_Matrix(const std::vector<std::vector<T>>& input_rows):
    size_(input_rows.size())
{
    for(const std::vector<T> &row : input_rows) {
        assert((row.size() == size_) && ("Invalid matrix size"));
    }

    if(size_ > 0) data_ = new T[size_ * size_];

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i * size_ + j] = input_rows[i][j];
        }
    }
}

template<typename T>
Square_Matrix<T>::~Square_Matrix() {
    if(size_ > 0) delete [] data_;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator= (const Square_Matrix<T>& matr)& {
    if(this != &matr) {
        if((size_ > 0)) delete [] data_;

        size_ = matr.size();

        if((size_ > 0)) data_ = new T[size_ * size_];

        for(size_t i = 0; i < size_; ++i) {
            for(size_t j = 0; j < size_; ++j) {
                data_[i * size_ + j] = matr(i, j);
            }
        }
    }

    return *this;
}

template<typename T>
template<typename U>
Square_Matrix<T>::Square_Matrix(const Square_Matrix<U>& matr):
    size_(matr.size())
{
    if(size_ > 0) data_ = new T[size_ * size_];
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i * size_ + j] = static_cast<T>(matr(i, j));
        }
    }
}

template<typename T>
void Square_Matrix<T>::transpose() const {
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < i; ++j) {
            std::swap(data_[i * size_ + j], data_[j * size_ + i]);
        }
    }
}

template<typename T>
void Square_Matrix<T>::add_row_to_row(size_t src_num, size_t dst_num) const {
    assert((src_num < size_) && "invalid row number");
    assert((dst_num < size_) && "invalid row number");

    for(size_t i = 0; i < size_; ++i) {
        data_[dst_num * size_ + i] += data_[src_num * size_ + i];
    }
}

template<typename T>
void Square_Matrix<T>::sub_row_to_row(size_t src_num, size_t dst_num) const {
    assert((src_num < size_) && "invalid row number");
    assert((dst_num < size_) && "invalid row number");

    for(size_t i = 0; i < size_; ++i) {
        data_[dst_num * size_ + i] -= data_[src_num * size_ + i];
    }
}

template<typename T>
T Square_Matrix<T>::operator()(size_t row_i, size_t column_i) const {
    assert((row_i < size_) && "invalid row");
    assert((column_i < size_) && "invalid column");
    return data_[row_i * size_ + column_i];
}

template<typename T>
bool Square_Matrix<T>::operator==(const Abstract_Matrix<T>& inp_rhs) const {
    if(typeid(*this) != typeid(inp_rhs)) return false; //IDE?

    const Square_Matrix<T>& rhs = static_cast<const Square_Matrix&>(inp_rhs);

    if(size_ != rhs.size()) return false;

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            if(!is_match(data_[i * size_ + j], rhs(i, j))) return false;
        }
    }

    return true;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator+=(const Square_Matrix<T>& rhs)& {
    assert((size_ == rhs.size()) && "different matrix sizes");
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i * size_ + j] += rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator-=(const Square_Matrix<T>& rhs)& {
    assert((size_ == rhs.size()) && "different matrix sizes");
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i * size_ + j] -= rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator*=(const Square_Matrix<T>& rhs)& { //naive algorithm (almost)
    assert((size_ == rhs.size()) && "different matrix sizes");
    if(size_ == 0) return *this;

    T *new_data = new T[size_ * size_];

    Square_Matrix new_rhs(rhs);
    new_rhs.transpose(); //for cache friendly

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            new_data[i * size_ + j] = 0;
            for(size_t k = 0; k < size_; ++k) {
                new_data[i * size_ + j] += (data_[i * size_ + k] * new_rhs(j, k));
            }
        }
    }

    delete [] data_;
    data_ = new_data;

    return *this;
}

template<typename T>
void Square_Matrix<T>::print() const {
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            std::cout << data_[i * size_ + j] << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T>
Square_Matrix<T> operator+(const Square_Matrix<T>& lhs, const Square_Matrix<T>& rhs) {
    Square_Matrix<T> tmp(lhs);
    tmp += rhs;
    return tmp;
}

template<typename T>
Square_Matrix<T> operator-(const Square_Matrix<T>& lhs, const Square_Matrix<T>& rhs) {
    Square_Matrix<T> tmp(lhs);
    tmp -= rhs;
    return tmp;
}

template<typename T>
Square_Matrix<T> operator*(const Square_Matrix<T>& lhs, const Square_Matrix<T>& rhs) {
    Square_Matrix<T> tmp(lhs);
    tmp *= rhs;
    return tmp;
}

template<>
double Square_Matrix<double>::determinant() const;

template<>
int Square_Matrix<int>::determinant() const;

template<>
long long int Square_Matrix<long long int>::determinant() const;



//----------------------------methods for Matrix----------------------------------------

template<typename T>
Matrix<T>::Matrix(const T* inp_data, size_t inp_column_size, size_t inp_row_size):
    column_size_(inp_column_size),
    row_size_(inp_row_size)
{
    if((row_size_ > 0) && (column_size_ > 0)) data_ = new T[row_size_ * column_size_];

    for(size_t i = 0; i < row_size_ * column_size_; ++i) {
        data_[i] = inp_data[i];
    }
}

template<typename T>
Matrix<T>::Matrix(const std::vector<std::vector<T>>& input_rows):
    column_size_(input_rows.size()),
    row_size_((input_rows.size() > 0) ? input_rows[0].size() : 0)
{
    for(const std::vector<T> &row : input_rows) {
        assert((row.size() == row_size_) && ("Invalid matrix size"));
    }

    if((column_size_ > 0) && (row_size_ > 0)) data_ = new T[row_size_ * column_size_];

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i * row_size_ + j] = input_rows[i][j];
        }
    }
}

template<typename T>
Matrix<T>::~Matrix() {
    if((column_size_ > 0) && (row_size_ > 0)) delete [] data_;
}

template<typename T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& matr)& {
    if(this != &matr) {
        if((column_size_ > 0) && (row_size_ > 0)) delete [] data_;

        column_size_ = matr.column_size();
        row_size_ = matr.row_size();

        if((column_size_ > 0) && (row_size_ > 0)) data_ = new T[row_size_ * column_size_];

        for(size_t i = 0; i < column_size_; ++i) {
            for(size_t j = 0; j < row_size_; ++j) {
                data_[i * row_size_ + j] = matr(i, j);
            }
        }
    }

    return *this;
}

template<typename T>
template<typename U>
Matrix<T>::Matrix(const Matrix<U>& matr):
    column_size_(matr.column_size()),
    row_size_(matr.row_size())
{
    if((column_size_ > 0) && (row_size_ > 0)) data_ = new T[row_size_ * column_size_];
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i * row_size_ + j] = static_cast<T>(matr(i, j));
        }
    }
}

template<typename T>
void Matrix<T>::transpose() {
    T* new_data = nullptr;
    if((column_size_ > 0) && (row_size_ > 0)) new_data = new T[row_size_ * column_size_];

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            new_data[j * column_size_ + i] = data_[i * row_size_ + j];
        }
    }

    if((column_size_ > 0) && (row_size_ > 0)) delete [] data_;

    data_ = new_data;
    std::swap(row_size_, column_size_);
}

template<typename T>
void Matrix<T>::add_row_to_row(size_t src_num, size_t dst_num) const {
    assert((src_num < column_size_) && "invalid row number");
    assert((dst_num < column_size_) && "invalid row number");

    for(size_t i = 0; i < row_size_; ++i) {
        data_[dst_num * row_size_ + i] += data_[src_num * row_size_ + i];
    }
}

template<typename T>
void Matrix<T>::sub_row_to_row(size_t src_num, size_t dst_num) const {
    assert((src_num < row_size_) && "invalid row number");
    assert((dst_num < row_size_) && "invalid row number");

    for(size_t i = 0; i < row_size_; ++i) {
        data_[dst_num * row_size_ + i] -= data_[src_num * row_size_ + i];
    }
}

template<typename T>
T Matrix<T>::operator()(size_t row_i, size_t column_i) const {
    assert((row_i < column_size_) && "invalid row");
    assert((column_i < row_size_) && "invalid column");
    return data_[row_i * row_size_ + column_i];
}

template<typename T>
bool Matrix<T>::operator==(const Abstract_Matrix<T>& inp_rhs) const {
    if(typeid(*this) != typeid(inp_rhs)) return false;

    const Matrix<T>& rhs = static_cast<const Matrix&>(inp_rhs);

    if((row_size_ != rhs.row_size()) || (column_size_ != rhs.column_size())) return false;

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            if(!is_match(data_[i * row_size_ + j], rhs(i, j))) return false;
        }
    }

    return true;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs)& {
    assert((row_size_ == rhs.row_size()) && "different matrix sizes");
    assert((column_size_ == rhs.column_size()) && "different matrix sizes");
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i * row_size_ + j] += rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs)& {
    assert((row_size_ == rhs.row_size()) && "different matrix sizes");
    assert((column_size_ == rhs.column_size()) && "different matrix sizes");
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i * row_size_ + j] -= rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs)& { //naive algorithm (almost)
    assert((row_size_ == rhs.column_size()) && "invalid matrix sizes");

    if((row_size_ == 0) || (column_size_ == 0)) return *this;

    if(rhs.row_size() == 0) {
        delete [] data_;
        data_ = nullptr;
        row_size_ = 0;
        return *this;
    }

    T *new_data = new T[column_size_ * rhs.row_size()];

    Matrix new_rhs(rhs);
    new_rhs.transpose(); //for cache friendly

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < new_rhs.column_size(); ++j) {
            new_data[i * new_rhs.column_size() + j] = 0;
            for(size_t k = 0; k < row_size_; ++k) {
                new_data[i * new_rhs.column_size() + j] += (data_[i * row_size_ + k] * new_rhs(j, k));
            }
        }
    }

    delete [] data_;
    data_ = new_data;
    row_size_ = rhs.row_size();

    return *this;
}

template<typename T>
void Matrix<T>::print() const {
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            std::cout << data_[i * row_size_ + j] << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T>
Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    Matrix<T> tmp(lhs);
    tmp += rhs;
    return tmp;
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    Matrix<T> tmp(lhs);
    tmp -= rhs;
    return tmp;
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    Matrix<T> tmp(lhs);
    tmp *= rhs;
    return tmp;
}

}
