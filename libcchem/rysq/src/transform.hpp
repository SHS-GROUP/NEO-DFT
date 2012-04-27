#ifndef _RYSQ_TRANSFORM_HPP_
#define _RYSQ_TRANSFORM_HPP_


template<class bra, int N, class Transform>
static void rysq::Quadrature::eval(const Shell::Impl &c, const Shell::Impl &d,
				   int K, const double* const *C, const double *Ints2d,
				   Transform &transform) {
     
    typedef double T;
    static const int ldN = (sizeof(T) > 8) ? N : N + N%2;
    static const int Nij = ldN*(bra::A::L+1)*(bra::B::L+1);
    const int dim2d = Nij*(c.L+1)*(d.L+1);

    const T *Ix = &Ints2d[0*dim2d];
    const T *Iy = &Ints2d[1*dim2d];
    const T *Iz = &Ints2d[2*dim2d];

    int spk = (c.type < 0);
    int spl = (d.type < 0);

    for(int l = d.first, kl = 0; l <= d.last; ++l) {
	const int lx = (c.L+1)*LX[l];
	const int ly = (c.L+1)*LY[l];
	const int lz = (c.L+1)*LZ[l];

	const int lsp = (spl && l) << spk;

	for(int k = c.first; k <= c.last; ++k, ++kl) {
	    const double *Ckl = C[(spk && k) + lsp];

	    const int klx = Nij*(lx + LX[k]);
	    const int kly = Nij*(ly + LY[k]);
	    const int klz = Nij*(lz + LZ[k]);

	    static const int ni = bra::A::size;
	    static const int nj = bra::B::size;
	    double I[ni*nj] __attribute__((aligned(16))) = { 0.0 };
	    
	    int flags = 0;
	    double screen = 0.0;
	    kernels::quadrature<bra,N>(flags, screen, K, Ckl, dim2d,
				       &Ix[klx], &Iy[kly], &Iz[klz], 1.0, I);

	    // 	    fock::eval<ni,nj>(k - c.first, l - d.first, kl,
	    // 			      D[0], D[1], D[2], D[3], D[4], D[5],
	    // 			      F[0], F[1], F[2], F[3], F[4], F[5], I);

	    transform(k - c.first, l - d.first, kl, I);

	}
    }

}


#endif /* _RYSQ_TRANSFORM_HPP_ */
