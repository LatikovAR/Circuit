#pragma once

#include <vector>
#include <utility>

#include "matrix.h"

namespace circuit {

struct Edge final {
private:
    size_t begin_;
    size_t end_;
    double R_;
    double U_;
    double I_;

    static bool is_prev_input_ok;
    static bool is_eof;
public:
    Edge(size_t begin, size_t end, double R, double U):
        begin_(begin), end_(end),
        R_(R), U_(U), I_(0.0) {}

    static Edge input_edge();

    static std::vector<Edge> input_edges();

    void print() const;

    static void print_edges(const std::vector<Edge>& edges);

    size_t begin() const { return begin_; }
    size_t end() const { return end_; }
    double R() const  { return R_; }
    double U() const { return U_; }
    double I() const { return I_; }

    void set_I(double I) { I_ = I; }
};


struct Vertex final {
private:
    std::vector<std::pair<Vertex*, const Edge*>> nodes_;
public:
    void add_node(Vertex* connected_vert, const Edge* edge) {
        nodes_.push_back(std::pair<Vertex*, const Edge*>(connected_vert, edge));
    }

    size_t nodes_num() const { return nodes_.size(); }

    Vertex* next_vertex(size_t i) const {
        if(i >= nodes_.size()) return nullptr;
        return nodes_[i].first;
    }

    const Edge* edge(size_t i) const {
        if(i >= nodes_.size()) return nullptr;
        return nodes_[i].second;
    }
};


class Circuit final {
private:
    std::vector<Edge> edges_;
    std::vector<Vertex> vertices_;

    void build_circuit_graph();
public:
    Circuit(const std::vector<Edge>& edges);
    void print_vertices_all() const;
};



}
