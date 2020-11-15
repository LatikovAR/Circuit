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




//---------------------------------------Edge methods----------------------------------------

Vertex* Edge::next_vertex(const Vertex* cur_vertex) {
    if(vertex1 != cur_vertex) return vertex1;
    if(vertex2 != cur_vertex) return vertex2;
    return nullptr;
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

            //we shouldn't check this vertices any more
            visited = true;

            Vertex* next_vertex = edge->next_vertex(this);

            assert((next_vertex != nullptr) &&
                   "Invalid graph: lone undefined edge can't be cycle");

            assert((next_vertex->num_edges_undefined_ > 0) &&
                   "Invalid graph: num_edges_undefined counter is invalid");

            (next_vertex->num_edges_undefined_)--;

            return next_vertex;
        }
    }

    assert(0 && "Invalid graph: num_edges_undefined counter is invalid");
    return nullptr;
}


int Vertex::define_vertex_outside_any_cycle_as_visited() {
    for(const Edge* edge : edges_) {
        if(edge->condition != OUT_OF_CYCLE) return -1;
    }

    visited = true;

    return 0;
}


namespace  { //for Vertex::find_cycle()
void add_cycle_to_data(std::pair<std::vector<std::pair<Vertex*, size_t>>, std::vector<Edge*>>& trace,
                       size_t begin_iterator,
                       std::vector<std::vector<std::pair<Vertex*, Edge*>>>& cycles_data,
                       Edge* last_edge) {
    std::vector<std::pair<Vertex*, Edge*>> cycle;

    for(size_t i = begin_iterator; i < (trace.first.size() - 1); ++i) {
        //pushing vertex and edge in cycle

        Vertex* pushed_vertex = trace.first[i].first;
        Edge* pushed_edge = trace.second[i];
        cycle.push_back(std::pair<Vertex*, Edge*>(pushed_vertex, pushed_edge));

        //marking edge before pushing vertex as IN_CYCLE
        trace.second[i]->condition = IN_CYCLE;
    }

    //last cycle element
    size_t pushed_vertex_iterator = trace.first.size() - 1;
    Vertex* pushed_vertex = trace.first[pushed_vertex_iterator].first;
    cycle.push_back(std::pair<Vertex*, Edge*>(pushed_vertex, last_edge));
    last_edge->condition = IN_CYCLE;

    cycles_data.push_back(cycle);
}
}


void Vertex::find_cycle(std::pair<std::vector<std::pair<Vertex*, size_t>>, std::vector<Edge*>>& trace,
                        std::vector<std::vector<std::pair<Vertex*, Edge*>>>& all_cycles) {
    assert(trace.first.size() > 0);

    while(1) {
        Vertex* cur_vertex = trace.first[trace.first.size() - 1].first;
        size_t edge_iterator = trace.first[trace.first.size() - 1].second;
        //I can't use size_t& edge_iterator because vector can make reallocations
        //so it's changing will be written by really strange way
        //but this places will be marked as //edge_iterator++ =)

        //end condition
        if((trace.first.size() == 1) && (cur_vertex->edges_num() == edge_iterator)) break;

        //if no more edges for current vertex, we should go back
        if(cur_vertex->edges_num() == edge_iterator) {
            trace.first.pop_back();
            trace.second.pop_back();
            continue;
        }

        Edge* cur_edge = cur_vertex->edges_[edge_iterator];

        //we mustn't go back to previous vertex with the same way
        if((trace.second.size() > 0) && (cur_edge == trace.second[trace.second.size() - 1])) {
            (trace.first[trace.first.size() - 1].second)++; //edge_iterator++;
            continue;
        }


        //if we know, that this edge already checked
        if(cur_edge->condition != UNDEFINED) {
            (trace.first[trace.first.size() - 1].second)++; //edge_iterator++;
            continue;
        }

        Vertex* next_vertex = cur_edge->next_vertex(cur_vertex);

        if(next_vertex->visited == true) {

            //searching next_vertex in trace
            size_t i = trace.first.size();
            while(i > 0) {
                --i;
                if(trace.first[i].first == next_vertex) {
                    add_cycle_to_data(trace, i, all_cycles, cur_edge);
                }
            }

            (trace.first[trace.first.size() - 1].second)++; //edge_iterator++;
            continue;
        }

        else {
            //going to the next vertex
            trace.second.push_back(cur_edge);
            (trace.first[trace.first.size() - 1].second)++; //edge_iterator++;
            trace.first.push_back(std::pair<Vertex*, size_t>(next_vertex, 0));
            next_vertex->visited = true;
            continue;
        }
    }
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
    std::cout << "All vertices:\n";
    for(size_t i = 0; i < vertices_.size(); ++i) {
        std::cout << "Vertex:\n";
        std::cout << "visited = " << vertices_[i].visited << std::endl;
        for(size_t j = 0; j < vertices_[i].edges_num(); ++j) {
            assert(vertices_[i].edge(j) != nullptr);
            vertices_[i].edge(j)->edge_info->print();
        }
        std::cout << std::endl;
    }
}

//for debug
void Circuit::print_edges_all() const {
    std::cout << "All edges:\n";
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

            if(vertices_[i].num_edges_undefined() == 0) {
                std::cout << "Warning: it looks like algorithm invalidation";

                int res = vertices_[i].define_vertex_outside_any_cycle_as_visited();
                assert((res == 0) && "Invalid num_edges_undefined() counter");
            }

            if(vertices_[i].num_edges_undefined() == 1) {

                //if we define edge of this vertex as OUT_OF_CYCLE
                //we should go and check another vertex, connected with this edge

                Vertex *next_vert = &(vertices_[i]);

                do {
                    next_vert = next_vert->define_lone_edge_as_out_of_cycle();
                    assert(next_vert != nullptr);
                } while(next_vert->num_edges_undefined() == 1);

                //if another vertex, connected with this edge, hasn't any more undefined edges
                if(next_vert->num_edges_undefined() == 0) {
                    int res = next_vert->define_vertex_outside_any_cycle_as_visited();
                    assert((res == 0) && "Invalid num_edges_undefined() counter");
                }
            }
        }
    }
}



void Circuit::find_cycles(std::vector<std::vector<std::pair<Vertex*, Edge*>>>& all_cycles) {
    all_cycles.clear();

    for(size_t i = 0; i < vertices_.size(); ++i) {
        if(vertices_[i].visited == false) {
            vertices_[i].visited = true;

            //this storage contain trace of the graph traversal in vertices and edges
            //for all vertices also contains iterators of edge,
            //which was selected in this traversal
            std::pair<std::vector<std::pair<Vertex*, size_t>>, std::vector<Edge*>> trace;
            trace.first.push_back(std::pair<Vertex*, size_t>(&(vertices_[i]), 0));
            vertices_[i].find_cycle(trace, all_cycles);
        }
    }
}


void Circuit::print_cycles_all(const std::vector<std::vector<std::pair<Vertex*, Edge*>>>& all_cycles) const {
    std::cout << "Cycles:\n";
    for(const std::vector<std::pair<Vertex*, Edge*>>& cycle : all_cycles) {
        if(cycle.size() == 0) continue;

        Edge* cur_edge = cycle[0].second;
        Edge* next_edge;

        if(cycle.size() == 1) {
            std::cout << cur_edge->edge_info->begin() << " -- " << cur_edge->edge_info->end();
        }
        else {
            size_t next_vert_num;
            size_t start_vert_num;
            next_vert_num = cur_edge->edge_info->begin();

            next_edge = cycle[1].second;
            if((next_edge->edge_info->begin() != next_vert_num) &&
               (next_edge->edge_info->end() != next_vert_num)) {
                std::cout << next_vert_num << " -- ";
                start_vert_num = next_vert_num;
                next_vert_num = cur_edge->edge_info->end();
            }
            else {
                std::cout << cur_edge->edge_info->end() << " -- ";
                start_vert_num = cur_edge->edge_info->end();
            }
            assert(((next_edge->edge_info->begin() == next_vert_num) ||
                    (next_edge->edge_info->end() == next_vert_num)) && "invalid cycle");

            for(size_t i = 1; i < cycle.size() - 1; ++i) {
                std::cout << next_vert_num << " -- ";
                cur_edge = next_edge;
                next_edge = cycle[i + 1].second;

                if(next_vert_num == cur_edge->edge_info->begin()) {
                    next_vert_num = cur_edge->edge_info->end();
                }
                else {
                    assert((next_vert_num == cur_edge->edge_info->end()) && "invalid cycle");
                    next_vert_num = cur_edge->edge_info->begin();
                }
            }

            std::cout << next_vert_num << " -- ";
            cur_edge = next_edge;
            if(next_vert_num == cur_edge->edge_info->begin()) {
                next_vert_num = cur_edge->edge_info->end();
            }
            else {
                assert((next_vert_num == cur_edge->edge_info->end()) && "invalid cycle");
                next_vert_num = cur_edge->edge_info->begin();
            }
            std::cout << next_vert_num;

            assert((next_vert_num == start_vert_num) && "invalid cycle");
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


void Circuit::find_all_currents() {
    check_elems_beyond_cycles(); //optional

    //all cycles will be contained here like a sequence of vertices and edges
    //vertex[0] -- edge[0] -- vertex[1] -- edge[1] --......
    //begin and end vertex/edge no matter
    std::vector<std::vector<std::pair<Vertex*, Edge*>>> all_cycles;

    //this method find all "linearly independent" cycles
    //"linearly independent" cycles - that cycles, which can't be maked from edges of the other cycles
    find_cycles(all_cycles);

    //we do this, because all edges in cycle already defined as IN_CYCLE
    //it also set I in this edges to 0
    define_all_undefined_edges_as_out_of_cycle();

    print_cycles_all(all_cycles);
}

void Circuit::define_all_undefined_edges_as_out_of_cycle() {
    for(size_t i = 0; i < edges_.size(); ++i) {
        if(edges_[i].condition == UNDEFINED) {
            edges_[i].condition = OUT_OF_CYCLE;
            edges_[i].edge_info->set_I(0.0);
         }
    }
}

}
