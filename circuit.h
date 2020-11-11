#pragma once

#include <vector>

#include "matrix.h"

namespace circuit {

struct Node final {
private:
    size_t begin_;
    size_t end_;
    double R_;
    double U_;
    double I_;

    static bool is_prev_input_ok;
public:
    Node(size_t begin, size_t end, double R, double U):
        begin_(begin), end_(end),
        R_(R), U_(U), I_(0.0) {}

    static Node input_node();

    void print() const;
};


class Circuit final {
private:
    std::vector<Node> nodes;
};



}
