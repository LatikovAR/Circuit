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
    static bool is_eof;
public:
    Node(size_t begin, size_t end, double R, double U):
        begin_(begin), end_(end),
        R_(R), U_(U), I_(0.0) {}

    static Node input_node();

    static std::vector<Node> input_nodes();

    void print() const;

    static void print_nodes(const std::vector<Node>& nodes);

    size_t begin() { return begin_; }
    size_t end() { return end_; }
    double R() { return R_; }
    double U() { return U_; }
    double I() { return I_; }

    void set_I(double I) { I_ = I; }
};


class Circuit final {
private:
    std::vector<Node> nodes_;
public:
    Circuit(std::vector<Node> nodes);
};



}
