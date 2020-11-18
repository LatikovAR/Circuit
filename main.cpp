#include <iostream>

#include "matrix.h"
#include "circuit.h"


#ifdef MATRIX_TEST
#include "matrix_test.h"
using namespace matrix_tests;
#endif

int main() {
#ifdef MATRIX_TEST
    //matrix_test();
    unit_test13();
    unit_test14();
    unit_test15();
#endif

    std::vector<circuit::Edge_Info> edges_info = circuit::Edge_Info::input_edges_info();
    circuit::Circuit circ(edges_info);
    return 0;
}
