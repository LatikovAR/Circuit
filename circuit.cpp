#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <set>
#include <cassert>

#include "circuit.h"
#include "dynamic_array.h"

namespace circuit {

//------------------------------------Edge_Info methods-------------------------------------
bool Edge_Info::is_prev_input_ok = true;
bool Edge_Info::is_eof = false;

Edge_Info Edge_Info::input_edge_info() {
    size_t begin = 0;
    size_t end = 0;
    double R = 0;
    double U = 0;

    std::string buf = {};
    std::string trash = {};

    is_prev_input_ok = true;
    is_eof = false;

    std::getline(std::cin, buf, '\n');
    if(std::cin.eof()) {
        is_eof = true;
        if(buf == "") { //sometimes it can be useful
            is_prev_input_ok = false;
            return Edge_Info(0, 0, 0.0, 0.0);
        }
    }
    if(std::cin.fail() && !(std::cin.eof())) {
        std::cout << "Warning: invalid str input format\n";
        is_prev_input_ok = false;
        return Edge_Info(0, 0, 0.0, 0.0);
    }

    std::stringstream s_buf;

    s_buf << buf;

    s_buf >> begin;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(0, 0, 0.0, 0.0);
    }

    s_buf >> trash;
    if(!(s_buf.good()) || (trash != "--")) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, 0, 0.0, 0.0);
    }

    s_buf >> end;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, 0, 0.0, 0.0);
    }

    s_buf >> trash;
    if(!(s_buf.good()) || (trash != ",")) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, end, 0.0, 0.0);
    }

    s_buf >> R;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, end, 0.0, 0.0);
    }

    s_buf >> trash;
    if((s_buf.eof() && s_buf) && (trash == ";")) {
        return Edge_Info(begin, end, R, 0.0);
    }

    if(!(s_buf.good()) || (trash != ";")) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, end, R, 0.0);
    }

    s_buf >> U;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, end, R, 0.0);
    }

    s_buf >> trash;
    if((s_buf.eof() && s_buf) && (trash == "V")) {
        return Edge_Info(begin, end, R, U);
    }

    std::cout << "Warning: invalid input format\n";
    is_prev_input_ok = false;
    return Edge_Info(begin, end, R, U);
}

std::vector<Edge_Info> Edge_Info::input_edges_info() {
    std::vector<Edge_Info> edges_info;

    while(is_eof == false) {
        Edge_Info edge_info = input_edge_info();
        if(is_prev_input_ok) {
            edges_info.push_back(edge_info);
        }
    }

    return edges_info;
}

void Edge_Info::print() const {
    std::cout << "Edge: " << begin_ << " -- " << end_ << std::endl;
    std::cout << "R = " << R_ << std::endl;
    std::cout << "U = " << U_ << std::endl;
    if(is_solved())
        std::cout << "I = " << I_ << "\n";
    else std::cout << "I = UNDEFINED\n";
    std::cout << std::endl;
}

void Edge_Info::print_edges_info(const std::vector<Edge_Info> &edges_info) {
    for(const Edge_Info& edge_info : edges_info) {
        edge_info.print();
    }
}



//---------------------------------------Vertex methods--------------------------------------

const Edge* Vertex::find_undefined_edge() const {
    for(const Edge* edge: edges_) {
        if(edge->condition == UNDEFINED) return edge;
    }
    return nullptr;
}


Vertex* Vertex::define_lone_edge_as_out_of_cycle() {
    if(num_edges_undefined_ != 1) return nullptr;

    for(Edge* edge: edges_) {
        if(edge->condition == UNDEFINED) {
            edge->condition = OUT_OF_CYCLE;
            num_edges_undefined_--;

            //if edge is out of cycle, it can't have any current
            edge->edge_info->set_I(0.0);

            if(this != edge->vertex1) {
                assert((edge->vertex1->num_edges_undefined_ > 0) &&
                       "Invalid graph: num_edges_undefined counter is invalid");

                (edge->vertex1->num_edges_undefined_)--;

                return edge->vertex1;
            }
            else {
                assert((this != edge->vertex2) &&
                       "Invalid graph: lone undefined edge can't be cycle");
                assert((edge->vertex2->num_edges_undefined_ > 0) &&
                       "Invalid graph: num_edges_undefined counter is invalid");

                (edge->vertex2->num_edges_undefined_)--;

                return edge->vertex2;
            }
        }
    }

    assert(0 && "Invalid graph: num_edges_undefined counter is invalid");
    return nullptr;
}



//---------------------------------------Circuit methods-------------------------------------

//for Circuit constructor
namespace  {
size_t num_of_vertices(const std::vector<Edge_Info>& edges_info) {
    std::set<size_t> vertices_nums;
    for(const Edge_Info& edge_info : edges_info) {
        vertices_nums.insert(edge_info.begin());
        vertices_nums.insert(edge_info.end());
    }
    return vertices_nums.size();
}
}

Circuit::Circuit(const std::vector<Edge_Info>& edges_info):
    edges_info_(edges_info),
    vertices_(num_of_vertices(edges_info)),
    edges_(edges_info.size())
{
    build_circuit_graph();
    print_vertices_all();
    find_all_currents();
    print_vertices_all();
}



//for set in the build_circuit_graph()
namespace  {
struct comp {
    bool operator() (const std::pair<size_t, Vertex*>& lhs,
                     const std::pair<size_t, Vertex*>& rhs) const {
        //std::cout << 1 << std::endl; - for my interest =)
        return (lhs.first < rhs.first);
    }
};
}

void Circuit::build_circuit_graph() {
    //getting all numbers of existing vertices
    std::set<size_t> vertices_nums;
    for(size_t i = 0; i < edges_info_.size(); ++i) {
        vertices_nums.insert(edges_info_[i].begin());
        vertices_nums.insert(edges_info_[i].end());
    }

    assert(vertices_.size() == vertices_nums.size());

    //building set of vertices_nums connected with vector of vertices
    std::set<std::pair<size_t, Vertex*>, comp> vertices_connected_nums;
    size_t i = 0;
    for(const size_t& elem : vertices_nums) {
        vertices_connected_nums.insert(std::pair<size_t, Vertex*>(elem, &(vertices_[i])));
        ++i;
    }



    //pushing edges and vertices in circuit
    for(size_t i = 0; i < edges_info_.size(); ++i) {
        Edge_Info& edge_info = edges_info_[i];
        //getting pointers to connected vertices
        auto iter1 = vertices_connected_nums.find(std::pair<size_t, Vertex*>(edge_info.begin(), nullptr));
        auto iter2 = vertices_connected_nums.find(std::pair<size_t, Vertex*>(edge_info.end(), nullptr));
        Vertex* vert1 = iter1->second;
        Vertex* vert2 = iter2->second;

        //creating edge
        edges_[i] = Edge(vert1, vert2, &edge_info);

        //connecting vertices with edge
        vert1->add_edge(&(edges_[i]));
        vert2->add_edge(&(edges_[i]));
    }
}

//for debug
void Circuit::print_vertices_all() const {
    for(size_t i = 0; i < vertices_.size(); ++i) {
        for(size_t j = 0; j < vertices_[i].edges_num(); ++j) {
            assert(vertices_[i].edge(j) != nullptr);
            vertices_[i].edge(j)->edge_info->print();
        }
        std::cout << std::endl;
    }
}

//for debug
void Circuit::print_edges_all() const {
    for(size_t i = 0; i < edges_.size(); ++i) {
        edges_[i].edge_info->print();
    }
}


//it's simple algorithm which separates only vertices with <=1 edge and their edges
//of course it isn't all possible cases of non-cycle elems
//but this step is only for optimisation and can be deleted
//
//all separated vertices will be marked as visited
//all separated edges will be marked as OUT_OF_CYCLE
void Circuit::check_elems_beyond_cycles() {
    for(size_t i = 0; i < vertices_.size(); ++i) {

        if(vertices_[i].visited == false) { //unnecessary to check already separated vertex

            if(vertices_[i].num_edges_undefined() == 1) {

                //if we define edge of this vertex as OUT_OF_CYCLE
                //we should go and check another vertex, connected with this edge

                Vertex *next_vert = &(vertices_[i]);

                do {
                    next_vert = next_vert->define_lone_edge_as_out_of_cycle();
                    assert(next_vert != nullptr);
                } while(next_vert->num_edges_undefined() == 1);
            }
        }
    }
}

void Circuit::find_all_currents() {
    check_elems_beyond_cycles(); //optional
}

}
