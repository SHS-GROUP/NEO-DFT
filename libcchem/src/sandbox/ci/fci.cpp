#include "runtime.hpp"
#include "ci.hpp"

#include <assert.h>
#include <vector>
#include <boost/multi_array.hpp>



template<typename T>
struct bitset;

template<typename T>
std::ostream& operator<<(std::ostream &os, const bitset<T> &S);


template<typename T>
struct bitset {

    explicit bitset(size_t size = 0) {
	assert(max() >= size);
	size_ = size;
	count_= 0;
	value_ = 0;
    }

    explicit bitset(const std::vector<bool> &bits) {
	assert(max() >= bits.size());
	size_ = bits.size();
	count_ = 0;
	value_ = 0;
	BOOST_AUTO(it, bits.rbegin());
	for (size_t i = 0; i < bits.size(); ++i) {
	    count_ += *it;
	    set_(i, *it);
	    ++it;
	}
    }

    size_t size() const { return size_; }
    size_t count() const { return count_; }
    T value() const { return value_; }

    bool at(size_t i) const {
	assert(i < size_);
	return (value_ & bit(i));
    }

    size_t next_zero(size_t i) const {
	while (value_ & bit(i)) ++i;
	return i;
    }

    void swap(size_t i, size_t j) {
	assert(at(i) ^ at(j));
	flip_(i);
	flip_(j);
    }

    static size_t max() { return sizeof(T)*8; }

protected:

    size_t size_, count_;
    T value_;

    void set_(size_t i, T value = 1) {
	assert(i < size_);
	value_ = (value_ | bit(i, value));
    }

    void flip_(size_t i) {
	value_ = (value_ ^ bit(i));
    }

    static T bit(size_t i, T value = 1) {
	return (T((bool)value) << i);
    }

};

template<typename T>
std::ostream& operator<<(std::ostream &os, const bitset<T> &S) {
    os << "(";
    for (size_t i = S.size(); i > 0; --i) {
	os << S.at(i-1);
	if (i > 1) os << ", ";
    }
    os << "), value = " << S.value();
    return os;
}


struct String : bitset<size_t> {

    typedef size_t type;
    typedef bitset<type> value_type;

    explicit String(size_t no, size_t ne) {
	this->size_ = no;
	this->count_ = ne;
	for (size_t i = 0; i < ne; ++i) this->set_(i);
    }

    explicit String(const bitset<type> &bits)
	: bitset<type>(bits) {}

    value_type value() const { return *this; }

    bool operator[](size_t i) const {
	return this->at(i);
    }

    size_t index() const {
	return bitset<type>::value_;
    }

    size_t bits() const {
	return bitset<type>::value_;
    }

    static std::vector<String> permutations(size_t no, size_t ne) {
	std::vector<bool> bits(no-ne, 0);
	bits.resize(no, 1);
	std::vector<String> strings;
	do {
	    String::value_type bs(bits);
	    //std::cout << bs << std::endl;
	    strings.push_back(String(bs));
	} while (std::next_permutation(bits.begin(), bits.end()));
	return strings;
    }

};



struct excitations {

    template<class H, class G, class V>
    static void evaluate12(const String &I,
			   const H &h, const G &g,
			   const V &C, V &ab) {
	apply<evaluate1>(I, h, g, C, ab);
    }

private:

    template< template<class,class,class> class F,
	      class H, class G, class V>
    static F<H,G,V> apply(const String &I,
			  const H &h, const G &g,
			  const V &C, V &ab) {
	F<H,G,V> f = { h, g, C, ab, 0 };
	for_each(I, f);
	return f;
    }

    template<class F>
    static void for_each(const String &I, F &f) {
	size_t no = I.size();
	size_t ne = I.count();
	int J[no-ne];
	for (size_t k = 0, j = 0; k < no-ne; ++k, ++j) {
	    j = I.next_zero(j);
	    J[k] = j;
	}
	for (size_t e = 0, i = 0; e < ne; ++e, ++i) {
	    while (!I[i]) ++i;
	    // std::cout << i << std::endl;
	    for (size_t k = 0; k < no-ne; ++k) {
		size_t j = J[k];
		//std::cout << i << "->" << j << std::endl;
		f(I, i, j);
	    }
	}
	//throw;
    }

    static double sgn(size_t i, size_t j) {
	return 1;
    }

    template<class H, class G, class V>
    struct evaluate1 {
	const H &h;
	const G &g;
	const V &C;
	V &ab;
	double value;
	void operator()(const String &I, size_t k, size_t l) {
	    String::value_type value = I.value();
	    value.swap(k,l);
	    String J(value);
	    double c1 = h[k][l];
	    double c2 = apply<evaluate2>(J, h, g[k][l], C, ab).value;
	    ab[I.index()] += sgn(k,l)*(c1 + 0.5*c2)*C[J.index()];
	}
    };

    template<class H, class G, class V>
    struct evaluate2 {
    	const H &h;
    	const G &g;
    	const V &C;
    	V &ab;
    	double value;
    	void operator()(const String &I, size_t i, size_t j) {
    	    value += sgn(i,j)*g[i][j];
    	}
    };

};


void fci(size_t no, size_t na, size_t nb) {

    boost::multi_array<double,4> V(boost::extents[no][no][no][no+1]);
    boost::multi_array<double,2> H(boost::extents[no][no]);

    const std::vector<String> &A = String::permutations(no,na);
    const std::vector<String> &B = String::permutations(no,nb);

    Runtime::Progress progress;
    progress.reset(A.size()*B.size());

    utility::timer t;

#pragma omp parallel
    {
	std::vector<double> C(1<<no);
	std::vector<double> ab(1<<no);

#pragma omp for schedule(dynamic,1)
	for (size_t i = 0; i < A.size(); ++i) {
	    // load C(i), ab(i)
	    for (size_t j = 0; j < B.size(); ++j) {
		const String &b = B.at(j);
		//std::cout << B << std::endl;
		excitations::evaluate12(b, H, V, C, ab);
#pragma omp critical
		progress++;
	    }
	}
    }

    std::cout << "time: " << t << std::endl;

}
