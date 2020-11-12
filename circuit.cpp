#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "circuit.h"

namespace circuit {

//------------------------------------Node methods-------------------------------------
bool Node::is_prev_input_ok = true;
bool Node::is_eof = false;

Node Node::input_node() {
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
    if(std::cin.fail()) {
        std::cout << "Warning: invalid str input format\n";
        is_prev_input_ok = false;
        return Node(0, 0, 0.0, 0.0);
    }

    std::stringstream s_buf;

    s_buf << buf;

    s_buf >> begin;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Node(0, 0, 0.0, 0.0);
    }

    s_buf >> trash;
    if(!(s_buf.good()) || (trash != "--")) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Node(begin, 0, 0.0, 0.0);
    }

    s_buf >> end;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Node(begin, 0, 0.0, 0.0);
    }

    s_buf >> trash;
    if(!(s_buf.good()) || (trash != ",")) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Node(begin, end, 0.0, 0.0);
    }

    s_buf >> R;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Node(begin, end, 0.0, 0.0);
    }

    s_buf >> trash;
    if((s_buf.eof() && s_buf) && (trash == ";")) {
        return Node(begin, end, R, 0.0);
    }

    if(!(s_buf.good()) || (trash != ";")) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Node(begin, end, R, 0.0);
    }

    s_buf >> U;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Node(begin, end, R, 0.0);
    }

    s_buf >> trash;
    if((s_buf.eof() && s_buf) && (trash == "V")) {
        return Node(begin, end, R, U);
    }

    std::cout << "Warning: invalid input format\n";
    is_prev_input_ok = false;
    return Node(begin, end, R, U);
}

std::vector<Node> Node::input_nodes() {
    std::vector<Node> nodes;

    while(is_eof == false) {
        Node node = input_node();
        if(is_prev_input_ok) {
            nodes.push_back(node);
        }
    }

    return nodes;
}

void Node::print() const {
    std::cout << "Node: " << begin_ << " -- " << end_ << std::endl;
    std::cout << "R = " << R_ << std::endl;
    std::cout << "U = " << U_ << std::endl;
    std::cout << std::endl;
}

void Node::print_nodes(const std::vector<Node>& nodes) {
    for(const Node& node : nodes) {
        node.print();
    }
}



//---------------------------------------Circuit methods-------------------------------------

Circuit::Circuit(std::vector<Node> nodes):
    nodes_(nodes)
{

}

}
