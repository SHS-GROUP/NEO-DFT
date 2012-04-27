#include <boost/thread/tss.hpp>
#include <boost/typeof/typeof.hpp>

#include "rysq/cuda.hpp"
#include "cuda/detail.hpp"
#include "transpose.hpp"
#include "foreach.hpp"

#include "boost/cuda/runtime.hpp"
#include "cuda/matrix.hpp"

using namespace rysq;
namespace host = boost::cuda;

rysq::cuda::Fock::Fock() : kernel_() {}

rysq::cuda::Fock::Fock(const rysq::Quartet<rysq::Shell> &quartet)
{
    bool t0 = (!quartet[0].L && quartet[1].L);
    bool t1 = (!quartet[2].L && quartet[3].L);
    bool t01 = (!quartet.bra().L() && quartet.ket().L());
    Transpose transpose(t0, t1, t01);

    try { this->kernel_.reset(new cuda::detail::Fock(quartet, transpose)); }
    catch (std::exception&) { this->kernel_.reset(NULL); }

    for (int i = 0; i < 4; ++i)	block_[i] = quartet[i].size;
}

rysq::cuda::Fock::~Fock() {
    boost::cuda::thread::synchronize();
}

void rysq::cuda::Fock::swap(Fock &other) {
    std::swap(transpose_, other.transpose_);
    std::swap(block_, other.block_);
    if (kernel_.get()) kernel_->reset();
    if (other.kernel_.get()) other.kernel_->reset();
    std::swap(kernel_, other.kernel_);
} 

void rysq::cuda::Fock::operator()(const Centers &centers,
				  const Quartets &quartets,
				  density_matrix_set D, fock_matrix_set F,
				  const Parameters &parameters) {
    detail::matrix_set<double> S(D, F);

    if (!stream_.get()) stream_.reset(new boost::cuda::stream());

    size_t sizes[6] = { 0 };
    {
	size_t size = 0;
	int i = 0;
	foreach (const size_t (&index)[2], index_list()){
	    sizes[i] = S.num_blocks[index[0]]*S.num_blocks[index[1]];
	    size += sizes[i];
	    ++i;
	}
	if (!mutex_.get()) mutex_.reset(new boost::cuda::vector<int>());
	if (size > mutex_.get()->size()) {
	    boost::cuda::thread::synchronize();
	    mutex_.get()->resize(std::max(size, size_t(1 << 16)), 0);
	}
    }

    int* mutex[6];
    for (int i = 0; i < 6; ++i) {
	size_t size = sizes[i];
	if (i == 0) mutex[i] = mutex_.get()->begin().get();
	if (i < 5) mutex[i+1] = mutex[i] + size;
	for (int j = 0; j < i; ++j) {
	    if (S.fock(i) == S.fock(j)) mutex[i] = mutex[j];
	}
    }

    stream_.get()->synchronize();
    quartets_.assign(quartets);

    (*kernel_)(detail::Centers(centers.begin()),
	       detail::Quartets(quartets_.begin(), quartets_.size()),
	       S, mutex, parameters, *stream_.get());

}
