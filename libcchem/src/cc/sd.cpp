#include "core/wavefunction.hpp"
#include "basis/basis.hpp"
#include "integral/eri.hpp"
#include "cc/cc.hpp"
#include "cc/tensor.hpp"
#include "blas.hpp"
#include "foreach.hpp"
#include "utility/timer.hpp"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/adaptor.hpp>

namespace cc {
namespace sd {

    template<class A, typename R>
    void copy(const double *shell, A &a,
	      const R &p, const R &q, const R &r, const R &s) {
	for (size_t l = 0; l < s.size(); ++l) {
	    for (size_t k = r.start(); k < r.stop(); ++k) {
		for (size_t j = q.start(); j < q.stop(); ++j) {
		    for (size_t i = p.start(); i < p.stop(); ++i) {
	    		    // printf("zzz %i %i %i %i %e\n", l,j,k,i,
	    		    // 	   *shell);
			a[l][j][k][i] = *shell++;
		    }
		}
	    }
	}
    }

    void vvvv(const Wavefunction &wf, Map<const Array*> t,
	      Array &u2, Array &c) {

	utility::timer timer;

	namespace ublas = boost::numeric::ublas;
	using ublas::make_matrix;
	using ublas::column_major;

	const Basis &basis = wf.basis();
	int N = basis.size();
	int no = wf.active().size();
	int nv = wf.virtuals().size();
	//Matrix Co = wf.C(wf.active());   // C(i,p)
	Matrix Cv = wf.C(wf.virtuals()); // C(p,a)

	for (int a = 0; a < nv; ++a) {
	    for (int i = 0; i < N; ++i) {
		std::cout <<"hhh " << i << " " << a << " " << Cv(i,a) << std::endl;
	    }
	}

	struct {
	    std::vector<double> data[3];
	    double* operator[](int i) { return &data[i][0]; }
	} buffer;
	buffer.data[0].resize(no*no*N*basis.max().size());
	buffer.data[1].resize(N*N*no);
	buffer.data[2].resize(std::max<int>(nv*nv*no,
					    (N*N*N)*basis.max().size()));


	std::cout <<  "no " << no << ", nv " << nv << std::endl;

	// c(p,q,i,j) = 0.5*(t(i,j,a,b) + t(i,a)*t(j,b))*C(a,p)*C(b,q)
	{
	    Matrix t1(no,nv);
	    Matrix pb(N,nv), ab(nv,nv);
	    {
		size_t begin[] = { 0, 0 };
		size_t end[] =   { no, nv };
 		t("i,a")->get(t1.data().begin(), begin, end);
	    }
	    for (int j = 0; j < no; ++j) {
		size_t begin[] = { 0, 0, 0, j };
		size_t end[] = { nv, nv, no, j+1 };
		t("a,b,i,j")->get(buffer[2], begin, end);
		for (int i = 0; i < no; ++i) {
		    BOOST_AUTO(tij,
			       make_matrix<column_major>(nv, nv, buffer[2]+i*nv*nv));
		    BOOST_AUTO(cij,
			       make_matrix<column_major>(N, N, buffer[1]+i*N*N));
		    ab = tij + ublas::outer_prod(ublas::row(t1,i),
		    				 ublas::row(t1,j));
		    blas::gemm(0.5, Cv, ab, 0, pb); // c(p,b)
		    blas::gemm(1, pb, ublas::trans(Cv), 0, cij); // c(p,r)
		}
		end[0] = N;
		end[1] = N;
		c.put(buffer[1], begin, end);
	    }
	}


	integral::Screening screening(basis.shells().size());
	integral::Eri eri(basis, &screening);	

	foreach (const Basis::Shell &S, basis) {
	    //std::fill(buffer[2], buffer[2]+S.size()*S.stop()*N*N, 0);
	    boost::multi_array_ref<double,4>
		G(buffer[2], boost::extents[S.size()][S.stop()][N][N]);
	    foreach (const Basis::Shell &Q, basis.range(&S+1)) {

		utility::timer t;

	    	foreach (const Basis::Block &R, basis.blocks()) {
	    	    foreach (const Basis::Block &P, basis.blocks()) {
	    		eri(P, Q, R, S);
	    		size_t size = (P.shell().size()*R.shell().size()*
	    			       Q.data()->size()*S.data()->size());
	    		for (size_t i = 0; i < eri.quartets().size(); ++i) {
	    		    BOOST_AUTO(const &q, eri.quartets()[i]);
	    		    copy(eri.data() + i*size, G,
	    			 basis[q[0]], basis[q[1]],
	    			 basis[q[2]], basis[q[3]]);
	    		}
	    	    }
	    	}

		// std::cout << "shells (" << S.index() << "," << Q.index() << "): "
		// 	  << t
		// 	  << ", weight: " << ((S.data()->K()*Q.data()->K())*
		// 	   		      (S.data()->L()*Q.data()->L()+1))
		// 	  << ", size: " << (S.data()->size()*Q.data()->size())
		// 	  << std::endl;


	    }

	    // for (int s = 0, qs = 0; s < S.size(); ++s) {
	    // 	for (int q = 0; q < S.stop(); ++q, ++qs) {
	    // 	    for (int j = 0, ij = 0; j < N; ++j) {
	    // 		for (int i = 0; i < N; ++i, ++ij) {
	    // 		    printf("xxx %i %i %i %i %e\n", s,q,j,i,
	    // 			   G[s][q][j][i]);
	    // 		}
	    // 	    }
	    // 	}
	    // }

	    size_t M = S.size()*S.stop();
	    BOOST_AUTO(t, make_matrix<column_major>(no*no, M, buffer[0]));
	    BOOST_AUTO(g, make_matrix<column_major>(N*N, M, G.data()));

	    // c(p,r,i,j)'*g(p,r,q,s) -> u(i,j,q,s) 
	    for (int j = 0; j < no; ++j) {
		BOOST_AUTO(cj, make_matrix<column_major>(N*N, no, buffer[1]));
		BOOST_AUTO(tj, ublas::subrange(t, j*no, (j+1)*no, 0, t.size2())); 
		size_t begin[] = { 0, 0, 0, j };
		size_t end[] = { N, N, no, j+1 };
		c.get(cj.data().begin(), begin, end);
		blas::gemm(1, cj, g, 0, tj);
	    }
	    size_t begin[] = { 0, 0, 0, S.start() };
	    size_t end[] = { no, no, S.stop(), S.stop() };
	    u2.put(t.data().begin(), begin, end);

	    
	    //printf("shell stops %i\n", S.stop());

	    for (int s = S.start(), qs = 0; s < S.stop(); ++s) {
	    	for (int q = 0; q < S.stop(); ++q, ++qs) {
	    	    for (int j = 0, ij = 0; j < no; ++j) {
	    		for (int i = 0; i < no; ++i, ++ij) {
	    		    printf("xxx %i %i %i %i %e\n", q,s,i,j, t(ij,qs));
	    		}
	    	    }
	    	}
	    }
	
	}


	// (i,j,q,s) -> (i,j,a,s)
	for (int s = 0; s < N; ++s) {
	    BOOST_AUTO(t, make_matrix<column_major>(no*no, nv, buffer[0]));
	    BOOST_AUTO(u, make_matrix<column_major>(no*no, N, buffer[1]));
	    // (i,j,q,s)
	    if (s) {
		u2.get(u.data().begin(),
		       u2.index(0,0,0,s),
		       u2.index(no,no,s,s+1));
	    }
	    // (j,i,s,q)
	    {
		u2.get(u.data().begin()+(s*no*no),
		       u2.index(0,0,s,s),
		       u2.index(no,no,s+1,N));
	    }
	    // (i,j,q,s) = (j,i,s,q)
	    for (int k = s; k < N; ++k) {
		BOOST_AUTO(ui, make_matrix<column_major>(no, no, &u(0,k)));
		ui.assign(Matrix(ublas::trans(ui)));
	    }
	    blas::gemm(1, u, Cv, 0, t); 
	    u2.put(t.data().begin(),
	    	   u2.index(0,0,0,s),
	    	   u2.index(no,no,nv,s+1));
	}

	// (i,j,a,s) -> (i,j,a,b)
	for (int a = 0; a < nv; ++a) {
	    BOOST_AUTO(t, make_matrix<column_major>(no*no, nv, buffer[0]));
	    BOOST_AUTO(u, make_matrix<column_major>(no*no, N, buffer[1]));
	    size_t begin[] = { 0, 0, a, 0 };
	    size_t end[] = { no, no, a+1, N };
	    u2.get(u.data().begin(), begin, end);
	    blas::gemm(1, u, Cv, 0, t);
	    end[3] = nv;
	    u2.put(t.data().begin(), begin, end);
	}

	//throw;

	std::cout << "CCSD: v(ab,ef)*c(ef,ij): " << timer << std::endl;

    }
    

}
}
