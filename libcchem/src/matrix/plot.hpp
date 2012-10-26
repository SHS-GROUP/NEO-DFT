#ifndef _MATRIX_PLOT_HPP_
#define _MATRIX_PLOT_HPP_

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

namespace matrix {

    namespace plot {

	struct expression {
	    virtual double operator()(int i, int j) const = 0;
	    virtual int size1() const = 0;
	    virtual int size2() const = 0;
	};

	template<class AE>
	struct adapter : expression {
	    adapter(const AE &E) : E_(E) {}
	    virtual double operator()(int i, int j) const { return E_(i,j); }
	    virtual int size1() const { return E_.size1(); }
	     virtual int size2() const { return E_.size2(); }
	private: const AE &E_;
	};

	void surface_(const expression &E);

	template<class AE>
	void surface(const AE &E) { surface_(adapter<AE>(E)); }
		

	class Renderer;

	struct Surface {
	    Surface(double r0 = 0, double r1 = 0)
		: renderer_(NULL) { range_[0] = r0; range_[1] = r1; }
	    template<class AE>
	    void render(const AE &E) {
#ifdef HAVE_VTK
		if (!renderer_) initialize();
		render_(adapter<AE>(E));
#endif
	    }
	private:
	    void initialize();
	    void render_(const expression &E);
	    class Renderer *renderer_;
	    double range_[2];
 	    boost::thread thread_;
	};

    }
}

#endif /* _MATRIX_PLOT_HPP_ */
