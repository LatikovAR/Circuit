#pragma once

#include <vector>
#include <utility>

#include "matrix.h"
#include "dynamic_array.h"

namespace circuit {

struct Edge_Info final {
private:
    size_t begin_;
    size_t end_;
    double R_;
    double U_;
    double I_;
    bool solved_ = false; //true when I_ will be defined

    //for correct input of the edges
    static bool is_prev_input_ok;
    static bool is_eof;
public:
    Edge_Info() {}

    Edge_Info(size_t begin, size_t end, double R, double U):
        begin_(begin), end_(end),
        R_(R), U_(U) {}

    static Edge_Info input_edge_info();

    static std::vector<Edge_Info> input_edges_info();

    void print() const;

    static void print_edges_info(const std::vector<Edge_Info>& edges_info);

    size_t begin() const { return begin_; }
    size_t end() const { return end_; }
    double R() const  { return R_; }
    double U() const { return U_; }
    double I() const { return I_; } //undefined if solved_ == false
    bool is_solved() const { return solved_; }

    void set_I(double I) { I_ = I; solved_ = true; }
};



enum Edge_Condition { IN_CYCLE, OUT_OF_CYCLE, UNDEFINED };

class Vertex; //declaration for Edge

struct Edge final {
    Vertex* vertex1;
    Vertex* vertex2;
    Edge_Info* edge_info;
    Edge_Condition condition;

    Edge(): vertex1(nullptr), vertex2(nullptr), edge_info(nullptr) {}

    Edge(Vertex* vertex1_inp, Vertex* vertex2_inp, Edge_Info* edge_info_inp):
        vertex1(vertex1_inp),
        vertex2(vertex2_inp),
        edge_info(edge_info_inp),
        condition(UNDEFINED) {}
};



//This class describes a graph vertex
class Vertex final {
private:
    std::vector<Edge*> edges_; //edges from this vertex

    //counter of edges from graph cycles
    //after adding all edges supposed in_cycle for this counter before check
    size_t num_edges_undefined_ = 0;
public:
    bool visited = false; //for graph traversal

    //adding edges after starting solving process is UB
    void add_edge(Edge* edge) {
        edges_.push_back(edge);
        if (edge->condition == UNDEFINED) ++num_edges_undefined_;
    }

    size_t edges_num() const { return edges_.size(); }

    //returns nullptr if edge_num invalid
    const Edge* edge(size_t edge_num) const {
        if(edge_num >= edges_.size()) return nullptr;
        return edges_[edge_num];
    }

    size_t num_edges_undefined() const { return num_edges_undefined_; }

    const Edge* find_undefined_edge() const;

    //returns another Vertex of this lone edge (nullptr if can't be done)
    Vertex* define_lone_edge_as_out_of_cycle();
};


class Circuit final {
private:

    //Dynamic_array is non-standart container that used for Circuit data
    //Circuit data containers should have some specifications:
    //
    //1) They mustn't allow allocations after building, because circuit graph has to store
    //some pointers to its data. So non-const vector and its analogs can't be used
    //
    //2) They must allow to change their elements (with constant position in memory)
    //So const vector and its analogs can't be used
    //
    //3) A size of the containers is computing during runtime and unknown at the
    //compilation moment. So array can't be used.

    dyn_arr::dynamic_array<Edge_Info> edges_info_; //circuit information about edges
    dyn_arr::dynamic_array<Vertex> vertices_; //graph vertices contains here after build_circuit_graph()
    dyn_arr::dynamic_array<Edge> edges_; //graph edges contains here after build_circuit_graph()

    void build_circuit_graph();

    void check_elems_beyond_cycles();

    void find_all_currents();

public:
    Circuit(const std::vector<Edge_Info>& edges_info_);

    void print_vertices_all() const;
    void print_edges_all() const;
};



}
