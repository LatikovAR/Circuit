#pragma once

#include <cstdlib>
#include <new>
#include <vector>

namespace dyn_arr {

template <typename T>
class dynamic_array final {
private:
    const size_t size_;
    T* data_;
public:
    dynamic_array(size_t size): size_(size) {
        data_ = new T[size];
    }

    dynamic_array(const std::vector<T>& data): size_(data.size()) {
        data_ = new T[size_];

        for(size_t i = 0; i < size_; ++i) {
            data_[i] = data[i];
        }
    }

    ~dynamic_array() { delete [] data_; }

    dynamic_array(const dynamic_array& rhs) = delete;

    dynamic_array& operator=(const dynamic_array&) = delete;

    size_t size() const { return size_; }

    T& operator[](size_t i) { return data_[i]; } //undefined if i >= size_

    const T& operator[](size_t i) const { return data_[i]; } //undefined if i >= size_
};

}
