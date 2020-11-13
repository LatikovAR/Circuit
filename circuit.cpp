#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <set>
#include <cassert>

#include "circuit.h"

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
        if(buf == "") {
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
    std::vector<Edge_Info> edges;

    while(is_eof == false) {
        Edge_Info edge = input_edge_info();
        if(is_prev_input_ok) {
            edges.push_back(edge);
        }
    }

    return edges;
}

void Edge_Info::print() const {
    std::cout << "Edge: " << begin_ << " -- " << end_ << std::endl;
    std::cout << "R = " << R_ << std::endl;
    std::cout << "U = " << U_ << std::endl;
    std::cout << std::endl;
}

void Edge_Info::print_edges_info(const std::vector<Edge_Info> &edges_info) {
    for(const Edge_Info& edge_info : edges_info) {
        edge_info.print();
    }
}




//---------------------------------------Circuit methods-------------------------------------

Circuit::Circuit(const std::vector<Edge_Info>& edges_info):
    edges_info_(edges_info)
{
    build_circuit_graph();
    print_vertices_all();
    print_edges_all();
    find_all_currents();
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
    for(const Edge_Info& edge_info : edges_info_) {
        vertices_nums.insert(edge_info.begin());
        vertices_nums.insert(edge_info.end());
    }

    //initialising vector of vertices
    vertices_.resize(vertices_nums.size());

    //building set of vertices_nums connected with vector of vertices
    std::set<std::pair<size_t, Vertex*>, comp> vertices_connected_nums;
    size_t i = 0;
    for(const size_t& elem : vertices_nums) {
        vertices_connected_nums.insert(std::pair<size_t, Vertex*>(elem, &(vertices_[i])));
        i++;
    }

    //pushing edges and vertices in circuit
    for(Edge_Info& edge_info : edges_info_) {
        //getting pointers to connected vertices
        auto iter1 = vertices_connected_nums.find(std::pair<size_t, Vertex*>(edge_info.begin(), nullptr));
        auto iter2 = vertices_connected_nums.find(std::pair<size_t, Vertex*>(edge_info.end(), nullptr));
        Vertex* vert1 = iter1->second;
        Vertex* vert2 = iter2->second;

        //creating edge


        edges_.push_back(Edge(vert1, vert2, &edge_info));

        //connecting vertices with edge
        vert1->add_edge(&(edges_[edges_.size() - 1]));
        vert2->add_edge(&(edges_[edges_.size() - 1]));
    }
}

//for debug
void Circuit::print_vertices_all() const {
    for(const Vertex& vertex: vertices_) {
        for(size_t i = 0; i < vertex.edges_num(); ++i) {
            assert(vertex.edge(i) != nullptr);
            std::cout << vertex.edge(i)->edge_info->begin() << " -- " <<
                         vertex.edge(i)->edge_info->end() << ", R = " <<
                         vertex.edge(i)->edge_info->R() << "; U = " <<
                         vertex.edge(i)->edge_info->U() << "; I = ";

            if(vertex.edge(i)->edge_info->is_solved())
                std::cout << vertex.edge(i)->edge_info->I() << ";\n";

            else std::cout << "UNDEFINED\n";
        }
        std::cout << std::endl;
    }
}

//for debug
void Circuit::print_edges_all() const {
    for(const Edge& edge: edges_) {
        std::cout << edge.edge_info->begin() << " -- " <<
                     edge.edge_info->end() << ", R = " <<
                     edge.edge_info->R() << "; U = " <<
                     edge.edge_info->U() << "; I = ";

        if(edge.edge_info->is_solved())
            std::cout << edge.edge_info->I() << ";\n";

        else std::cout << "UNDEFINED\n";
    }
}


//it's simple algorithm which separates only vertices with <=1 edge and their edges
//of course it isn't all possible cases of non-cycle elems
//but this step is only for optimisation and can be deleted
//
//all separated vertices will be marked as visited
//all separated edges will be marked as OUT_OF_CYCLE
void Circuit::check_elems_beyond_cycles() {
    for(Vertex& vertex : vertices_) {

        if(vertex.visited == false) { //unnecessary to check already separated vertex

            if(vertex.edges_num_in_cycle_or_undefined() == 1) {

            }
        }
    }
}

void Circuit::find_all_currents() {
    check_elems_beyond_cycles(); //optional
}

}
