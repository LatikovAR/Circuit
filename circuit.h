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
    size_t number_;
    size_t second_number_; //numeration for solving equations, without already defined edges

    //for correct input of the edges
    static bool is_prev_input_ok;
    static bool is_eof;
public:

    Edge_Info() {}

    Edge_Info(size_t begin, size_t end, double R, double U, size_t number):
        begin_(begin), end_(end),
        R_(R), U_(U),  number_(number) {}

    static Edge_Info input_edge_info(size_t edge_number);

    static std::vector<Edge_Info> input_edges_info();

    void print() const;

    static void print_edges_info(const std::vector<Edge_Info>& edges_info);

    size_t begin() const { return begin_; }
    size_t end() const { return end_; }
    double R() const  { return R_; }
    double U() const { return U_; }
    double I() const { return I_; } //undefined if solved_ == false
    bool is_solved() const { return solved_; }
    size_t number() const { return  number_; }
    size_t second_number() const { return  second_number_; }

    void set_I(double I) { I_ = I; solved_ = true; }
    void set_second_number(size_t n) { second_number_ = n; }
};




class Circuit final {
private:
    enum Condition { IN_CYCLE, OUT_OF_CYCLE, UNDEFINED };

    class Vertex; //declaration for Edge

    struct Edge final {
        Vertex* vertex1;
        Vertex* vertex2;
        Edge_Info* edge_info;
        Condition condition;

        Edge() {}

        Edge(Vertex* vertex1_inp, Vertex* vertex2_inp, Edge_Info* edge_info_inp):
            vertex1(vertex1_inp),
            vertex2(vertex2_inp),
            edge_info(edge_info_inp),
            condition(UNDEFINED) {}

        //returns nullptr if edge is cycle
        Vertex* next_vertex(const Vertex* cur_vertex);
    };


    //This class describes a graph vertex
    class Vertex final {
    private:
        std::vector<Edge*> edges_; //edges from this vertex

        //counter of edges from graph cycles
        //after adding all edges supposed in_cycle for this counter before check
        size_t num_edges_undefined_ = 0;

        size_t number_;

        //for find cycle
        static void add_cycle_to_data(std::pair<std::vector<std::pair<Vertex*, size_t>>, std::vector<Edge*>>& trace,
                                      size_t begin_iterator,
                                      std::vector<std::vector<std::pair<Vertex*, Edge*>>>& cycles_data,
                                      Edge* last_edge);
    public:
        bool visited = false; //for graph traversal
        Condition condition = UNDEFINED;

        //adding edges after starting solving process is UB
        void add_edge(Edge* edge) {
            edges_.push_back(edge);
            if (edge->condition == UNDEFINED) ++num_edges_undefined_;
        }

        size_t edges_num() const { return edges_.size(); }

        size_t number() const { return number_; }

        void define_number(size_t num) { number_ = num; }

        //returns nullptr if edge_num invalid
        const Edge* edge(size_t edge_num) const {
            if(edge_num >= edges_.size()) return nullptr;
            return edges_[edge_num];
        }

        size_t num_edges_undefined() const { return num_edges_undefined_; }

        const Edge* find_undefined_edge() const;

        //returns another Vertex of this lone edge (nullptr if can't be done)
        Vertex* define_lone_edge_as_out_of_cycle();

        //0 - ok, 1 - vertex can be in cycle
        int define_vertex_outside_any_cycle_as_visited();

        static void find_cycle(std::pair<std::vector<std::pair<Vertex*, size_t>>, std::vector<Edge*>>& trace,
                               std::vector<std::vector<std::pair<Vertex*, Edge*>>>& all_cycles);

    };

    struct comp_for_build_circuit_graph;

    bool validity_ = true;

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

    std::vector<std::vector<std::pair<Vertex*, Edge*>>> all_cycles_; //in circuit graph

    void build_circuit_graph();

    void check_elems_beyond_cycles();

    //this method find all "linearly independent" cycles
    //"linearly independent" cycles - that cycles, which can't be maked from edges of the other cycles
    void find_cycles();

    std::pair<std::vector<double>, bool> make_and_solve_linear_cicruit_equations();

    void find_all_currents();

    //also set I in this edges to 0
    void define_all_undefined_elems_as_out_of_cycle();

    //returns number of undefined edges
    size_t set_second_numbers_in_edge_info();
public:
    Circuit(const std::vector<Edge_Info>& edges_info_);

    void print_vertices_all() const; //for debug
    void print_edges_all() const; //for debug

    //for debug, and also check cycles validity
    void print_cycles_all() const;

    bool validity() const { return validity_; }

    //Undefined if validity() = false
    void print_circuit() const ;
};

} //namespace circuit
