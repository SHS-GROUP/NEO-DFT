#ifndef CC_TRIPLES_THREAD_HPP
#define CC_TRIPLES_THREAD_HPP

#include "cc/triples/triples.hpp"
#include "cc/triples/energy.hpp"
#include "cc/triples/clear.hpp"

#include <utility>
#include <boost/typeof/typeof.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <boost/array.hpp>

#include <boost/thread/thread.hpp>
#include <boost/thread/barrier.hpp>

#include "foreach.hpp"
#include "array/h5.hpp"
#include "array/permute.hpp"


namespace cc {
namespace triples {
namespace detail {

    template<class T3>
    struct Thread {

	struct Shared {
	    const size_t no, nv;
	    Matrix t2, Vjk, Vkj;
	    Matrix vov;
	    Matrix vv[9], vo[1];
	    Queue queue;
	    Shared(size_t no, size_t nv)
		: no(no), nv(nv), queue(nv)
	    {
		t2.resize(nv, nv, false);
		Vjk.resize(no, nv, false);
		Vkj.resize(no, nv, false);
		vov.resize(nv, no*nv, false);
		foreach(Matrix &m, vv) {
		    m.resize(nv, nv, false);
		}
		foreach(Matrix &m, vo) {
		    m.resize(nv, no, false);
		}
	    }
	};

	struct Team : boost::noncopyable {
	    Shared shared;
	    const size_t size;
	    boost::mutex mutex;
	    Task::Tuple task;
	    T3 t3;
	    Team(size_t size, size_t no, size_t nv, const T3 &t3 = T3())
		: shared(no, nv), size(size), t3(t3), barrier_(size) {}
	    void wait() { barrier_.wait(); }
	private:
	    boost::barrier barrier_;
	};

    public:
	int id, device;
	bool master;
	Correction C;
	Matrix u1;

    protected:
	size_t no_, nv_;
	Shared *shared_;
	Team *team_;

    public:
	Thread(Team *team, int id = 0, int device = -1) {
	    team_ = team;
	    shared_ = &team->shared;
	    this->id = id;
	    this->device = device;
	    this->master = (this->id == 0);

	    no_ = shared_->no;
	    nv_ = shared_->nv;

	    u1.resize(nv_,no_,false);
	    u1.clear();
	}
	virtual ~Thread() {};

	Team& team() const { return *team_; }

	void operator()(const Data &data, T3 &t3, Task::Tuple index) {
	    team_->wait();
	    clear(t3, this->id, this->team_->size);
	    team_->wait();
	    ijk(data, t3, index[0], index[1], index[2]);
	    energy(data, t3, index[0], index[1], index[2]);
	}

    protected:

	void ijk(const Data &data, T3 &t3,
		 int i, int j, int k) {

	    typedef std::pair<bool,int> terms;
	    size_t num_threads = team_->size;
	    
	    ijk1(data, t3, j, k, i, terms(0,1));
	    array::permute<1,0,2>(t3, this->id, num_threads); // cba

	    // cba
	    ijk1(data, t3, k, j, i, terms(1,1));
	    array::permute<0,2,1>(t3, this->id, num_threads); // cab
	    team_->wait();
	    array::permute<1,0,2>(t3, this->id, num_threads); // acb
	    
	    // acb
	    ijk1(data, t3, i,  k, j, terms(1,0));
	    array::permute<0,2,1>(t3, this->id, num_threads); // abc

	    // abc ijk, jik
	    ijk1(data, t3, i, j, k, terms(1,1));
	}

	void ijk1(const Data &data, T3 &t3,
		  int i, int j, int k,
		  std::pair<bool,int> terms) {
	    int no = data.no;
	    int nv = data.nv;

	    team_->wait();

	    if (terms.first && this->master) {
		data["t2(abij)"].get(shared_->t2.data().begin(),
				    index<4>(0,0,j,i), index<4>(nv,nv,j+1,i+1));
	    }

	    if (terms.second && this->master) {
		data["t2(aibj)"].get(shared_->vov.data().begin(),
				     index<4>(0,0,0,i), index<4>(nv,no,nv,i+1));
		size_t start[] = { 0,0,j,k };
		size_t stop[] = {  no, nv,j+1,k+1 };
		data["V(iajk)"].get(shared_->Vjk.data().begin(), start, stop);
		std::swap(start[2], start[3]);
		std::swap(stop[2], stop[3]);
		data["V(iajk)"].get(shared_->Vkj.data().begin(), start, stop);
	    }

	    terms.second *= 2;
	    while (terms.first || terms.second > 0) {
		if (this->master) shared_->queue.reset(nv);
		team_->wait();
		this->ijk1_impl(data,
				shared_->t2, shared_->vov, shared_->Vjk, shared_->Vkj,
				t3, i, j, k, terms, shared_->queue);
		team_->wait();
		terms.first = 0;
		terms.second--;
	    }
	}

	virtual
	void ijk1_impl(const Data &data,
		       const Matrix &tab, const Matrix &taib,
		       const Matrix &Vjk, const Matrix &Vkj,
		       T3 &t3,
		       int i, int j, int k,
		       std::pair<bool,int> terms,
		       Queue &queue) {
	    assert("method must be implemented" && 0);
	}

	void energy(const Data &data, const T3 &t3, 
		    int i, int j, int k) {
	    int no = data.no;
	    int nv = data.nv;

	    BOOST_AUTO(&vv, shared_->vv);
	    BOOST_AUTO(&t1, shared_->vo[0]);

	    if (this->master) {
		for (int q = 0; q < 9; ++q) {
		    const int indices[][2] = {{i,j}, {i,k}, {j,k},
					      {i,j}, {i,k}, {j,k},
					      {j,i}, {k,i}, {k,j}};
		    size_t r1[] = { 0, 0, indices[q][0], indices[q][1] };
		    size_t r2[] = { nv, nv, r1[2]+1, r1[3]+1 };
		    data[(q < 3) ? "t2(abij)" : "V(abij)"]
			.get(vv[q].data().begin(), r1, r2);
		}
		data["t1(ai)"]
		    .get(t1.data().begin(), index<2>(0,0), index<2>(nv, no));
	    }

	    Permutations<const Matrix&,3> t2
		(vv[0], vv[1], vv[2]);
	    Permutations<const Matrix&,6> vvoo
		(vv[3], vv[4], vv[5], vv[6], vv[7], vv[8]);

	    team_->wait();
	    energy_impl(data, t1, t2, vvoo, t3, i, j, k);
	    team_->wait();
	}

	virtual
	void energy_impl(const Data &data,
			 const Matrix &t1,
			 const Permutations<const Matrix&,3> t2,
			 const Permutations<const Matrix&,6> &V,
			 const T3 &t3,
			 int i, int j, int k) {
	    int nv = data.nv;
	    int num_threads = this->team_->size;
	    const int BLOCK = 16;
	    const int N = nv/BLOCK + (nv%BLOCK > 0);
	    typedef Energy E;
	    for (int c = 0; c < N; ++c) {
		if (c%num_threads != this->id) continue;
		// std::cout <<  this->id << std::endl;
		for (int b = 0; b < N; ++b) {
		    for (int a = 0; a < N; ++a) {
			using std::min;
			this->C += E::evaluate<BLOCK>
			    (t1, t2, V, t3, data.eh, data.ep, this->u1,
			     E::Index<'i'>(i), E::Index<'j'>(j), E::Index<'k'>(k),
			     E::Range(a*BLOCK, min(nv,(a+1)*BLOCK)),
			     E::Range(b*BLOCK, min(nv,(b+1)*BLOCK)),
			     E::Range(c*BLOCK, min(nv,(c+1)*BLOCK)));
		    }
		}
	    }
	    //std::cout << "xxx: " << this->C.etd << std::endl;
	}


    };


} // namespace detail
} // namespace triples
} // namespace cc

#endif // CC_TRIPLES_THREAD_HPP
