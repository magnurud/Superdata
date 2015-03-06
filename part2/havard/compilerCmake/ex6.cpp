#include <iostream>
using std::cout;
using std::endl;

int main(){
    cout << "Running" << endl;
#ifdef HAVE_TE
    cout << "HAVE_TE" << endl;
#endif 
#ifdef COMP_GNU
    cout << "GNU" << endl;
#endif 
#ifdef COMP_INTEL
    cout << "INTEL" << endl;
#endif 

    return 0;
}

