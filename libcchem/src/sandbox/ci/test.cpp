#include "runtime.hpp"

int main(int args, char **argv) {

    extern void fci(size_t no, size_t na, size_t nb);

    fci(16, 8, 8);

}
