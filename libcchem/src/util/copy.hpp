#ifndef _UTIL_COPY_HPP_
#define _UTIL_COPY_HPP_

namespace util {

    template<typename D, typename S, size_t N>
    void copy(D (&dst)[N], const S *src) {
	for (int i = 0; i < N; ++i) {
	    dst[i] = src[i];
	}
    }

    template<typename T, typename F, typename I>
    void copy(const I (&dims)[2], T *to, int ldt, const F *from, int ldf) {
	for (int j = 0; j < dims[1]; ++j) {
	    for (int i = 0; i < dims[0]; ++i) {
		to[i] = from[i];
	    }
	    to += ldt;
	    from += ldf;
	}
    }

    template<typename T, typename F, typename I>
    void unpack(const I (&dims)[2], T *to, int ldim, const F *from) {
	copy(dims, to, ldim, from, dims[0]);
    }

    template<typename T, typename F, typename I>
    void pack(const I (&dims)[2], T *to, const F *from, int ldim) {
	copy(dims, to, dims[0], from, ldim);
    }

}

#endif /* _UTIL_COPY_HPP_ */
