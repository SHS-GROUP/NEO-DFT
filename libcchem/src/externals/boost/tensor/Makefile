includes = -I../ -I../../ -I/opt/boost/include
libraries = -lblas

#CXX = icc

test: test.cpp
	$(CXX) -g -O0 -Wall -o test $(includes) $(libraries) test.cpp