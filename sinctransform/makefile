
# Makefile for sinc transform programs

CXX=g++
CURRENT=/Users/hannah/Documents/Summer2017/Flatiron/sinctransform
CURRENT_FLAG=-I$(CURRENT)
FINUFFT=/Users/hannah/Documents/Summer2017/Flatiron/fi2/finufft
FINUFFT_LIB_PATH=/lib/libfinufft.a
FINUFFT_LIB=$(FINUFFT)$(FINUFFT_LIB_PATH)
FINUFFT_HEADER_FLAG=-I$(FINUFFT)/src/
FINUFFT_HEADER=$(FINUFFT)/src/finufft.h
FINUFFT_FLAGS=-lfftw3 -lm
FLAGS=-std=c++11 -g -Wall 

# Main sinc computations
sinc1d.o: sinc1d.cpp fastgl.hpp sinctransform.hpp sincutil.hpp $(FINUFFT_HEADER)
	$(CXX) $(FLAGS) -fPIC -c sinc1d.cpp $(FINUFFT_HEADER_FLAG)

sinc2d.o: sinc2d.cpp fastgl.hpp sinctransform.hpp sincutil.hpp $(FINUFFT_HEADER)
	$(CXX) $(FLAGS) -fPIC -c sinc2d.cpp $(FINUFFT_HEADER_FLAG)

sinc3d.o: sinc3d.cpp fastgl.hpp sinctransform.hpp sincutil.hpp $(FINUFFT_HEADER)
	$(CXX) $(FLAGS) -fPIC -c sinc3d.cpp $(FINUFFT_HEADER_FLAG)

# For computation of Gauss-Legendre quadrature weights
fastgl.o: fastgl.cpp
	$(CXX) $(FLAGS) -c fastgl.cpp

# Utility functions (not needed for sinc1d, sinc2d, or sinc3d)
sincutil.o: sincutil.cpp sincutil.hpp
	$(CXX) $(FLAGS) -c sincutil.cpp

# Library itself (static), with necessary sinc and Gauss-Legendre object files
libsinc.a: sinc1d.o sinc2d.o sinc3d.o fastgl.o
	ar rcs libsinc.a sinc1d.o sinc2d.o sinc3d.o fastgl.o

# For direct computation of the sinc transforms
directsinc.o: directsinc.cpp sincutil.hpp
	$(CXX) $(FLAGS) -c directsinc.cpp

# Simple examples of usage in 1d, 2d, 3d
examples/example1d: libsinc.a sincutil.o examples/example1d.cpp
	$(CXX) $(FLAGS) -o examples/example1d sincutil.o examples/example1d.cpp libsinc.a $(CURRENT_FLAG) $(FINUFFT_LIB) $(FINUFFT_FLAGS)

examples/example2d: libsinc.a sincutil.o examples/example2d.cpp
	$(CXX) $(FLAGS) -o examples/example2d sincutil.o examples/example2d.cpp libsinc.a $(CURRENT_FLAG) $(FINUFFT_LIB) $(FINUFFT_FLAGS)

examples/example3d: libsinc.a sincutil.o examples/example3d.cpp
	$(CXX) $(FLAGS) -o examples/example3d sincutil.o examples/example3d.cpp libsinc.a $(CURRENT_FLAG) $(FINUFFT_LIB) $(FINUFFT_FLAGS)

# Longer test functions using random inputs, over many requested precisions
tests/test1d: libsinc.a sincutil.o directsinc.o tests/test1d.cpp 
	$(CXX) $(FLAGS) -o tests/test1d sincutil.o directsinc.o tests/test1d.cpp libsinc.a $(CURRENT_FLAG) $(FINUFFT_LIB) $(FINUFFT_FLAGS)

tests/test2d: libsinc.a sincutil.o directsinc.o tests/test2d.cpp 
	$(CXX) $(FLAGS) -o tests/test2d sincutil.o directsinc.o tests/test2d.cpp libsinc.a $(CURRENT_FLAG) $(FINUFFT_LIB) $(FINUFFT_FLAGS)

tests/test3d: libsinc.a sincutil.o directsinc.o tests/test3d.cpp 
	$(CXX) $(FLAGS) -o tests/test3d sincutil.o directsinc.o tests/test3d.cpp libsinc.a $(CURRENT_FLAG) $(FINUFFT_LIB) $(FINUFFT_FLAGS)

clean:
	rm -f sinc*.o fastgl.o libsinc.a directsinc.o examples/example*d tests/test*d

