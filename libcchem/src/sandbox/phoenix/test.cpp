#include "phoenix/thread.hpp"

int main() {
    namespace phoenix = boost::phoenix;
    BOOST_AUTO(i, phoenix::local_names::_i);
    BOOST_AUTO(j, phoenix::local_names::_j);

    namespace thread = boost::phoenix::thread;
    using phoenix::range;

    thread::group(8)
	(
	 thread::parallel_for((j = range(30), i = range(j+1)))
	 [thread::critical[std::cout << i << " " << j+1 << std::endl]]
	 );

}
