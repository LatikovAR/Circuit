#include <iostream>

#include "matrix.h"
#include "circuit.h"

#define MATRIX_TEST

#ifdef MATRIX_TEST
#include "matrix_test.h"
using namespace matrix_tests;
#endif

int main() {
#ifdef MATRIX_TEST
    unit_test0();
    unit_test1();
    unit_test2();
    unit_test3();
    unit_test4();
    unit_test5();
    unit_test6();
    unit_test7();
    unit_test8();
    unit_test9();
    unit_test10();
    unit_test11();
    unit_test12();
    long_test(10);
#endif

    //circuit::Node n = circuit::Node::input_node();
    //n.print();
    return 0;
}
