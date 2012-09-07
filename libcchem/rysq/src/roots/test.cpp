#include "roots/roots.hpp"

int main(int args, const char **argv) {
    rysq_roots_initialize();
#ifdef INTEGER64
    int64_t N = 12;
#else
    int32_t N = 12;
#endif
    double X = 0;
    double t[N], w[N];
    rysq_roots_(&N, &X, t, w);

    printf("N = %i, X = %e\n", N, X);

    for (int i = 0; i < N; ++i) {
	printf("t = %e, w = %e\n", t[i], w[i]);
    }

}
