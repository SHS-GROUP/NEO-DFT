
void hf::fock::kernel(const Basis &basis,
		      const BlockMatrix &D, BlockMatrix &F, const double *I,
		      int a, int b, int c, int d) {

    int ni = basis(a).size();
    int nj = basis(b).size();
    int nk = basis(c).size();
    int nl = basis(d).size();

    typedef BlockMatrix::Index Index;

    const double *Dij = D[Index(a,b)];
    const double *Dkl = D[Index(c,d)];
    const double *Dik = D[Index(a,c)];
    const double *Dil = D[Index(a,d)];
    const double *Djl = D[Index(b,d)];
    const double *Djk = D[Index(b,c)];

    double *Fij = F[Index(a,b)];
    double *Fkl = F[Index(c,d)];
    double *Fik = F[Index(a,c)];
    double *Fil = F[Index(a,d)];
    double *Fjl = F[Index(b,d)];
    double *Fjk = F[Index(b,c)];

    for (int l = 0, kl = 0, ijkl = 0; l < nl; ++l) {

        for (int k = 0; k < nk; ++k, ++kl) {
	    int jk = k*nj;
	    int jl = l*nj;

            for (int j = 0, ij = 0; j < nj; ++j, ++jk, ++jl) {
		int ik = k*ni;
		int il = l*ni;

                for (int i = 0; i < ni; ++i, ++ij, ++ik, ++il, ++ijkl) {
                    double valueI = I[ijkl];
		    //if(fabs(valueI) < 1e-10) continue;

                    Fij[ij] += Dkl[kl] * valueI*4.0;
		    Fkl[kl] += Dij[ij] * valueI*4.0;
		    Fik[ik] -= Djl[jl] * valueI;
		    Fil[il] -= Djk[jk] * valueI;
		    Fjl[jl] -= Dik[ik] * valueI;
		    Fjk[jk] -= Dil[il] * valueI;
                }
            }

        }
    }
    //std::cout << F << std::endl;
}
