CXXFLAGS = -g -O0 -Wall -I ../../../..
CXXFLAGS += -I /opt/cuda/include
LDFLAGS = -L /opt/cuda/lib64 -lcudart -lcublas

test: test.cpp

clean:
	rm test

check-syntax:
	g++ -Wall -fsyntax-only -S ${CHK_SOURCES} ${CXXFLAGS}
