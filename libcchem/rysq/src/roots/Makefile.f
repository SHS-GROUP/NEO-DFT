roots:
	$(CXX) $(CXXFLAGS) -I../ -c roots.cpp
	$(AR) $(ARFLAGS) libroots.a roots.o
	$(CXX) $(CXXFLAGS) -I../ test.cpp libroots.a -o test

