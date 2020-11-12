#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <set>

#include "circuit.h"

namespace circuit {

//------------------------------------Node methods-------------------------------------
bool Edge::is_prev_input_ok = true;
bool Edge::is_eof = false;

Edge Edge::input_edge() {
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
    }
    if(std::cin.fail() && !(std::cin.eof())) {
        std::cout << "Warning: invalid str input format\n";
        is_prev_input_ok = false;
        return Edge(0, 0, 0.0, 0.0);
    }

    std::stringstream s_buf;

    s_buf << buf;

    s_buf >> begin;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge(0, 0, 0.0, 0.0);
    }

    s_buf >> trash;
    if(!(s_buf.good()) || (trash != "--")) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge(begin, 0, 0.0, 0.0);
    }

    s_buf >> end;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge(begin, 0, 0.0, 0.0);
    }

    s_buf >> trash;
    if(!(s_buf.good()) || (trash != ",")) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge(begin, end, 0.0, 0.0);
    }

    s_buf >> R;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge(begin, end, 0.0, 0.0);
    }

    s_buf >> trash;
    if((s_buf.eof() && s_buf) && (trash == ";")) {
        return Edge(begin, end, R, 0.0);
    }

    if(!(s_buf.good()) || (trash != ";")) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge(begin, end, R, 0.0);
    }

    s_buf >> U;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge(begin, end, R, 0.0);
    }

    s_buf >> trash;
    if((s_buf.eof() && s_buf) && (trash == "V")) {
        return Edge(begin, end, R, U);
    }

    std::cout << "Warning: invalid input format\n";
    is_prev_input_ok = false;
    return Edge(begin, end, R, U);
}

std::vector<Edge> Edge::input_edges() {
    std::vector<Edge> nodes;

    while(is_eof == false) {
        Edge node = input_edge();
        if(is_prev_input_ok) {
            nodes.push_back(node);
        }
    }

    return nodes;
}

void Edge::print() const {
    std::cout << "Node: " << begin_ << " -- " << end_ << std::endl;
    std::cout << "R = " << R_ << std::endl;
    std::cout << "U = " << U_ << std::endl;
    std::cout << std::endl;
}

void Edge::print_edges(const std::vector<Edge>& edges) {
    for(const Edge& edge : edges) {
        edge.print();
    }
}



//---------------------------------------Circuit methods-------------------------------------

Circuit::Circuit(const std::vector<Edge>& edges):
    edges_(edges)
{
    build_circuit_graph();
    print_vertices_all();
}

//for build_circuit_graph
namespace  {
struct comp {
    bool operator() (const std::pair<size_t, Vertex*>& lhs,
                     const std::pair<size_t, Vertex*>& rhs) const {
        return (lhs.first < rhs.first);
    }
};
}

void Circuit::build_circuit_graph() {
    //getting all numbers of existing vertices
    std::set<size_t> vertices_nums;
    for(const Edge& edge : edges_) {
        vertices_nums.insert(edge.begin());
        vertices_nums.insert(edge.end());
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

    //pushing nodes in vector of vertices
    //using default comparator for pair
    for(const Edge& edge : edges_) {
        auto iter_begin = vertices_connected_nums.find(std::pair<size_t, Vertex*>(edge.begin(), nullptr));
        auto iter_end = vertices_connected_nums.find(std::pair<size_t, Vertex*>(edge.end(), nullptr));
        (*iter_begin).second->add_node((*iter_end).second, &edge);
        (*iter_end).second->add_node((*iter_begin).second, &edge);
    }
}

void Circuit::print_vertices_all() const {
    for(const Vertex& vertex: vertices_) {
        for(size_t i = 0; i < vertex.nodes_num(); ++i) {
            std::cout << vertex.edge(i)->begin() << " -- " <<
                         vertex.edge(i)->end() << ", R = " <<
                         vertex.edge(i)->R() << "; U = " <<
                         vertex.edge(i)->U() << "; I = " <<
                         vertex.edge(i)->I() << ";\n";
        }
    }
}

}
