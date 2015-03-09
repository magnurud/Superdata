#include <iostream>
#include "matrix.hh"

using namespace std;


int main(int argc, char** argv) {
#ifdef COMP_GNU
    cout << "!!! GNU !!!" << endl;
#endif
#ifdef COMP_INTEL
    cout << "!!! INTEL !!!" << endl;
#endif

    size_t cols = 4;
    size_t rows = cols*cols;
    Matrix<double, ColMajor> A(rows, cols);
    int a = 0;
    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) {
            A(i, j) = a++;
        }
    }
    A.print();
    A.mkTransposed();
    cout << endl;
    A.print();




    return 0;
}
