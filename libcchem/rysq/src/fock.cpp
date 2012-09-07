#include "rysq-fock.hpp"

#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/to_tuple.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/lambda/lambda.hpp>
#include "externals/boost/numeric/ublas/storage_adaptors.hpp"

#include "foreach.hpp"

#include "kernel/new.hpp"
#include "fock-transform.hpp"
#include "transpose.hpp"

using namespace rysq;

typedef kernel::Eri<> kernel_type;

Fock::Fock(const Quartet<Shell> &quartet) {

    transpose_ = Transpose(quartet[0].size < quartet[1].size,
 			   quartet[2].size < quartet[3].size,
 			   quartet.bra().size() < quartet.ket().size());
    kernel_ = kernel::new_<fock::Transform>(transpose_(quartet));
    // std::cout << transpose_(quartet) << std::endl;
}

Fock::~Fock() { delete static_cast<kernel_type*>(kernel_); }

// void Fock::operator()(const Quartet<Center> &centers,
// 		      density_set D, fock_set F,
// 		      const Parameters &parameters) {
//     // (*pImpl)(centers, D, F, parameters);
// }

template<typename T, size_t N, typename U, class F>
boost::array<T,N> make_array(U (&input)[N], F unary) {
    boost::array<T,N> array;
    std::transform(input, input + N, array.begin(), unary);
    return array;
}

template<class I, class O, class F>
void transform2(I begin,I end, O output, F binary) {
    std::transform(begin, end, output, output, binary);
}

namespace ublas =  boost::numeric::ublas;

template<typename T>
struct kernel_data {
    typedef typename boost::remove_pointer<T>::type value_type;
    typedef std::pair<int,int> index_type;
    // matrix types
    typedef ublas::column_major layout;
    typedef ublas::const_array_ref_adaptor<double> adaptor;
    typedef ublas::matrix<double, layout, adaptor> matrix;
    struct block : matrix {
	void check_range(size_t size1, size_t size2) const {
	    if (size1*size2 != this->size1()*this->size2()) {
		std::cout << size1 << " " << this->size1() << ", "
			  << size2 << " " << this->size2() << std::endl;
		throw std::range_error("");
	    }
	}
	void swap(size_t size1, size_t size2, const double *data = 0) {
	    if (data) check_range(size1, size2);
	    matrix tmp(size1, size2, adaptor(size1*size2, data));
	    matrix::swap(tmp);
	}
	void assign(size_t size1, size_t size2, const double *data) {
	    check_range(size1, size2);
	    if (data == this->data().begin()) return;
	    swap(size1, size2, data);
	    using namespace boost::lambda;
	    for_each(_2 = _1);
	}
	template<class F>
	void for_each(const F &f) {
	    double *__restrict__ data_ = const_cast<double*>(data().begin());
    	    double *__restrict__ buffer_ = const_cast<double*>(buffer);
	    size_t m = size1(), n = size2();
	    if (trans) {
		for (size_t j = 0; j < n; ++j) {
		    for (size_t i = 0; i < m; ++i) {
			f(*data_++, buffer_[j+i*n]);
		    }
		}
	    }
	    else {
		for (size_t i = 0; i < m*n; ++i) {
		    f(*data_++, *buffer_++);
		}
	    }
	}
	T buffer;
	bool trans;
        index_type index;
    };
    typedef block* iterator;
    typedef const block* const_iterator;
    kernel_data() {
	end_ = data_;
	for (int i = 0; i < 4; ++i) {
	    std::fill(map_[i], map_[i] + 4, (block*)NULL);
	}
    }
    iterator begin() { return data_; }
    iterator end() { return end_; }
    block& operator[](index_type ij) {
	block* &b = map_[ij.first][ij.second];
	if (!b) {
	    b = end_++;
	    b->index = ij;
	}
	return *b;
    }
    size_t size(){ return data_.size(); }
    operator boost::array<T,6>() {
	boost::array<T,6> array;
	for (size_t i = 0; i < array.size(); ++i) {
	    //index_type index(index_list::at(i)[0], index_list::at(i)[1]);
	    array[i] = data_[i].buffer;//(*this)[index].buffer;
	}
	return array;
    }
private:
    block data_[16];
    block* map_[4][4];
    block *end_;
};

void Fock::operator()(const std::vector<Center> &centers,
		      const std::vector<Int4> &quartets,
		      density_matrix_set D, fock_matrix_set F,
		      const Parameters &parameters) {

    typedef kernel_data<const double*> density_data;
    typedef kernel_data<double*> fock_data;

    kernel_type &kernel = *static_cast<kernel_type*>(kernel_);
    density_data density_;
    fock_data fock_;

    std::vector<double> buffer;
    size_t pos = 0;

    {
	size_t size = 0;
	for (size_t k = 0; k < index_list::size(); ++k) {
	    size_t i = index_list::at(k)[0];
	    size_t j = index_list::at(k)[1];

	    std::pair<int,int> ij(transpose_.index[i], transpose_.index[j]);
	    bool trans = (ij.first > ij.second);
	    if (trans) std::swap(ij.first, ij.second);

	    density_data::block &density_block = density_[ij];
	    fock_data::block &fock_block = fock_[ij];
	    density_block.trans = trans;
	    fock_block.trans = trans;

	    i = ij.first;
	    j = ij.second;

	    fock_block.swap(F.get(i,j)->block().size1(),
			    F.get(i,j)->block().size2());
	    density_block.swap(D.get(i,j)->block().size1(),
			       D.get(i,j)->block().size2());
	    size += fock_block.data().size();
	    size += density_block.data().size();
	}

	buffer.resize(size);
	size_t i = 0;
	foreach (density_data::block &block, density_) {
	    block.buffer = &buffer.at(i);
	    i += block.data().size();
	}
	pos = i;
	foreach (fock_data::block &block, fock_) {
	    block.buffer = &buffer.at(i);
	    i += block.data().size();
	}
    }

    foreach (Int4 q, quartets) {
 
#define BLOCK_SIZE(A,m) (A).get(i,j)->block().size ## m()
#define BLOCK(A) (A).get(i,j)->block(q[i] - (A).base_index[i],	\
				     q[j] - (A).base_index[j])

	foreach (density_data::block &block, density_) {
	    const size_t i = block.index.first, j = block.index.second;
	    block.assign(BLOCK_SIZE(D,1), BLOCK_SIZE(D,2), BLOCK(D));
	}
	foreach (fock_data::block &block, fock_) {
	    const size_t i = block.index.first, j = block.index.second;
	    // std::cout << i << "," << j << " "
	    // 	      << block.index.first << "," << block.index.second << std::endl;
	    block.swap(BLOCK_SIZE(F,1), BLOCK_SIZE(F,2), BLOCK(F));
	}

#undef BLOCK_SIZE
#undef BLOCK

	std::fill(buffer.begin() + pos, buffer.end(), 0);

	int i = q[0], j = q[1], k = q[2], l = q[3];
	Quartet<Center> C(centers[i], centers[j], centers[k], centers[l]);
	double scale = 1.0/Quartet<Shell>::symmetry(i, j, k, l);

	fock::Transform<>::Data data(density_, fock_);
	kernel(transpose_(C), data, Eri::Parameters(parameters.cutoff));

	foreach (fock_data::block &block, fock_) {
	     fock_data::index_type &index = block.index;
	     int ij = index.first + index.second;
	     double f = ((ij == 1 || ij == 5)
			 ? parameters.cscale : parameters.xscale);
	     using namespace boost::lambda;
	     block.for_each(_1 = _1 + f*scale*_2);
	 }

    }
}

void rysq::Fock::eval(const Quartet<type> &quartet,
		      const std::vector<Int4> &quartets,
		      density_matrix_set D, fock_matrix_set F,
		      const double *eri, const Parameters &parameters) {

    foreach (const Int4 &q, quartets) {
	int i = q[0], j = q[1], k = q[2], l = q[3];

	density_set D6;
	fock_set F6;
	size_t s = 0;		

	// submatrix D_, F_;
	foreach (const size_t (&index)[2], index_list()) {
	    int i = index[0], j = index[1];
	    const double *dij = D.get(i,j)->block(q[i] - D.base_index[i],
						  q[j] - D.base_index[j]);
	    double *fij = F.get(i,j)->block(q[i] - F.base_index[i],
					    q[j] - F.base_index[j]);

	    F6[s] =  fij;
	    D6[s] = const_cast<double*>(dij);
	    ++s;
	}

	Parameters p = parameters;
	p.scale(1.0/Quartet<Shell>::symmetry(i, j, k, l));
	eval(quartet, D6, F6, eri, p);

	eri += quartet.size();
    }

}

void rysq::Fock::eval(const Quartet<type> &quartet,
		      density_set D, fock_set F,
		      const double *Eri, const Parameters &p) {

    rysq::type a = quartet[0], b = quartet[1], c = quartet[2], d = quartet[3];

#define _TYPES	(rysq::SP)(rysq::S)(rysq::P)(rysq::D)(rysq::F)

#define _FOCK(r, types) if (a == BOOST_PP_SEQ_ELEM(0, types) &&	\
			    b == BOOST_PP_SEQ_ELEM(1, types)) {	\
	fock::eval<BOOST_PP_SEQ_ELEM(0, types),			\
	    BOOST_PP_SEQ_ELEM(1, types)>(c, d, D, F, Eri, p);	\
    } else

    BOOST_PP_SEQ_FOR_EACH_PRODUCT(_FOCK, (_TYPES)(_TYPES)) {
	std::cout <<  quartet << std::endl;
	throw std::exception();
    }

#undef _FOCK
#undef _TYPES

}
