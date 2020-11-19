CXX = g++
CXXFLAGS = 

WORK_SRC =  main.cpp matrix.cpp circuit.cpp
MATRIX_TEST_SRC =  main.cpp matrix.cpp matrix_test.cpp

work: $(WORK_SRC)
	$(CXX) $(CXXFLAGS) $? -DWORK -o prog
test: $(MATRIX_TEST_SRC)
	$(CXX) $(CXXFLAGS) $? -DMATRIX_TEST -o prog
