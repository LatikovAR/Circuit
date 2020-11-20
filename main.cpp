#include <iostream>

#include "matrix.h"

#define WORK

#ifdef WORK
#include "circuit.h"
#endif

#ifdef MATRIX_TEST
#include "matrix_test.h"
using namespace matrix_tests;
#endif

int main() {
#ifdef MATRIX_TEST
    matrix_test();
#endif
#ifdef WORK
    std::vector<circuit::Edge_Info> edges_info = circuit::Edge_Info::input_edges_info();
    circuit::Circuit circuit(edges_info);

    if(circuit.validity() == true) {
        circuit.print_circuit();
    }
    else {
        std::cout << "Input circuit can't be solved.\n";
        std::cout << "Maybe you set some R <= 0.0.\n";
    }
#endif
    return 0;
}
