# Makefile for sinc transform programs.
# By Hannah Lawrence. (C) Simons Foundation 2018.
# Edits by Alex Barnett.
#
# User: please edit the FLAGS and PATHS below for your system.
#
# Tasks:
#
# make lib      - make main static library
# make examples - make and run examples
# make tests    - make and run tests (takes ~30 s to run)
# make all
# make clean

CXX=g++
CURRENT=.
CURRENT_FLAG=-I$(CURRENT)

# point to your fftw library directory (here is the usual linux one):
FFTW=/usr/lib/x86_64-linux-gnu/
FFTW_FLAG=-L$(FFTW)

# point to the top of your finufft installation:
#FINUFFT=../finufft/
FINUFFT=/Users/hannah/Documents/Flatiron18/newfinufft/finufft
#/Users/hannah/Documents/Summer2017/Flatiron/fi2/finufft

#FINUFFT_LIB_PATH=/lib/libfinufft.so
FINUFFT_LIB_PATH=/lib-static/libfinufft.a

FINUFFT_LIB=$(FINUFFT)$(FINUFFT_LIB_PATH)
FINUFFT_HEADER_FLAG=-I$(FINUFFT)/src/
FINUFFT_HEADER=$(FINUFFT)/src/finufft.h
# if multi-thread:
FINUFFT_FLAGS=-lfftw3 -lfftw3_threads -lm -lgomp
# if single-thread:
#FINUFFT_FLAGS=-lfftw3 -lm

#FLAGS=-std=c++11 -g -Wall
# -fext-numeric-literals needed for 0+1i complex literals in gcc 5.4.0:
FLAGS=-std=c++11 -g -Wall -fext-numeric-literals

EXAMPLE_DIR=examples
TEST_DIR=tests
HELPER_DIR=helper
SOURCE_DIR=src

# build list of executables (could use such automation below too..)
EXAMPLES=$(patsubst %.cpp,%,$(wildcard $(EXAMPLE_DIR)/*.cpp))
TESTS=$(patsubst %.cpp,%,$(wildcard $(TEST_DIR)/*.cpp))

# tasks:
all: lib examples tests

lib: libsinc.a

examples: $(EXAMPLES)
	$(foreach x,$(EXAMPLES),$(x);)

tests: $(TESTS)
	$(foreach x,$(TESTS),$(x);)

# Main sinc computations
$(SOURCE_DIR)/sinc1d.o: $(SOURCE_DIR)/sinc1d.cpp helper/fastgl.hpp $(SOURCE_DIR)/sinctransform.hpp helper/sincutil.hpp $(FINUFFT_HEADER)
	$(CXX) $(FLAGS) -o $(SOURCE_DIR)/sinc1d.o -c $(SOURCE_DIR)/sinc1d.cpp $(FINUFFT_HEADER_FLAG) $(CURRENT_FLAG)/$(HELPER_DIR)

$(SOURCE_DIR)/sinc2d.o: $(SOURCE_DIR)/sinc2d.cpp helper/fastgl.hpp $(SOURCE_DIR)/sinctransform.hpp helper/sincutil.hpp $(FINUFFT_HEADER)
	$(CXX) $(FLAGS) -o $(SOURCE_DIR)/sinc2d.o -c $(SOURCE_DIR)/sinc2d.cpp $(FINUFFT_HEADER_FLAG) $(CURRENT_FLAG)/$(HELPER_DIR)

$(SOURCE_DIR)/sinc3d.o: $(SOURCE_DIR)/sinc3d.cpp helper/fastgl.hpp $(SOURCE_DIR)/sinctransform.hpp helper/sincutil.hpp $(FINUFFT_HEADER)
	$(CXX) $(FLAGS) -o $(SOURCE_DIR)/sinc3d.o -c $(SOURCE_DIR)/sinc3d.cpp $(FINUFFT_HEADER_FLAG) $(CURRENT_FLAG)/$(HELPER_DIR)

# For computation of Gauss-Legendre quadrature weights
$(HELPER_DIR)/fastgl.o: $(HELPER_DIR)/fastgl.cpp
	$(CXX) $(FLAGS) -o $(HELPER_DIR)/fastgl.o -c $(HELPER_DIR)/fastgl.cpp

# Utility functions (not needed for sinc1d, sinc2d, or sinc3d)
$(HELPER_DIR)/sincutil.o: $(HELPER_DIR)/sincutil.cpp $(HELPER_DIR)/sincutil.hpp
	$(CXX) $(FLAGS) -o $(HELPER_DIR)/sincutil.o -c $(HELPER_DIR)/sincutil.cpp

# Library itself (static), with necessary sinc and Gauss-Legendre object files
libsinc.a: $(SOURCE_DIR)/sinc1d.o $(SOURCE_DIR)/sinc2d.o $(SOURCE_DIR)/sinc3d.o $(HELPER_DIR)/fastgl.o
	ar rcs libsinc.a $(SOURCE_DIR)/sinc1d.o $(SOURCE_DIR)/sinc2d.o $(SOURCE_DIR)/sinc3d.o $(HELPER_DIR)/fastgl.o

# For direct computation of the sinc transforms
$(HELPER_DIR)/directsinc.o: $(HELPER_DIR)/directsinc.cpp $(HELPER_DIR)/sincutil.hpp
	$(CXX) $(FLAGS) -o $(HELPER_DIR)/directsinc.o -c $(HELPER_DIR)/directsinc.cpp

# Simple examples of usage in 1d, 2d, 3d
$(EXAMPLE_DIR)/example1d: libsinc.a $(HELPER_DIR)/sincutil.o $(EXAMPLE_DIR)/example1d.cpp
	$(CXX) $(FLAGS) -o $(EXAMPLE_DIR)/example1d $(HELPER_DIR)/sincutil.o $(EXAMPLE_DIR)/example1d.cpp libsinc.a $(CURRENT_FLAG)/$(HELPER_DIR) $(CURRENT_FLAG)/$(SOURCE_DIR) $(FINUFFT_LIB) $(FFTW_FLAG) $(FINUFFT_FLAGS)

$(EXAMPLE_DIR)/example2d: libsinc.a $(HELPER_DIR)/sincutil.o $(EXAMPLE_DIR)/example2d.cpp
	$(CXX) $(FLAGS) -o $(EXAMPLE_DIR)/example2d $(HELPER_DIR)/sincutil.o $(EXAMPLE_DIR)/example2d.cpp libsinc.a $(CURRENT_FLAG)/$(HELPER_DIR) $(CURRENT_FLAG)/$(SOURCE_DIR) $(FINUFFT_LIB) $(FFTW_FLAG) $(FINUFFT_FLAGS)

$(EXAMPLE_DIR)/example3d: libsinc.a $(HELPER_DIR)/sincutil.o $(EXAMPLE_DIR)/example3d.cpp
	$(CXX) $(FLAGS) -o $(EXAMPLE_DIR)/example3d $(HELPER_DIR)/sincutil.o $(EXAMPLE_DIR)/example3d.cpp libsinc.a $(CURRENT_FLAG)/$(HELPER_DIR) $(CURRENT_FLAG)/$(SOURCE_DIR) $(FINUFFT_LIB) $(FFTW_FLAG) $(FINUFFT_FLAGS)

# Longer test functions using random inputs, over many requested precisions
$(TEST_DIR)/test1d: libsinc.a $(HELPER_DIR)/sincutil.o $(HELPER_DIR)/directsinc.o $(TEST_DIR)/test1d.cpp 
	$(CXX) $(FLAGS) -o $(TEST_DIR)/test1d $(HELPER_DIR)/sincutil.o $(HELPER_DIR)/directsinc.o $(TEST_DIR)/test1d.cpp libsinc.a $(CURRENT_FLAG)/$(HELPER_DIR) $(CURRENT_FLAG)/$(SOURCE_DIR) $(FINUFFT_LIB) $(FFTW_FLAG) $(FINUFFT_FLAGS)

$(TEST_DIR)/test2d: libsinc.a $(HELPER_DIR)/sincutil.o $(HELPER_DIR)/directsinc.o $(TEST_DIR)/test2d.cpp 
	$(CXX) $(FLAGS) -o $(TEST_DIR)/test2d $(HELPER_DIR)/sincutil.o $(HELPER_DIR)/directsinc.o $(TEST_DIR)/test2d.cpp libsinc.a $(CURRENT_FLAG)/$(HELPER_DIR) $(CURRENT_FLAG)/$(SOURCE_DIR) $(FINUFFT_LIB) $(FFTW_FLAG) $(FINUFFT_FLAGS)

$(TEST_DIR)/test3d: libsinc.a $(HELPER_DIR)/sincutil.o $(HELPER_DIR)/directsinc.o $(TEST_DIR)/test3d.cpp 
	$(CXX) $(FLAGS) -o $(TEST_DIR)/test3d $(HELPER_DIR)/sincutil.o $(HELPER_DIR)/directsinc.o $(TEST_DIR)/test3d.cpp libsinc.a $(CURRENT_FLAG)/$(HELPER_DIR) $(CURRENT_FLAG)/$(SOURCE_DIR) $(FINUFFT_LIB) $(FFTW_FLAG) $(FINUFFT_FLAGS)

clean:
	rm -f $(SOURCE_DIR)/sinc*.o $(HELPER_DIR)/fastgl.o libsinc.a $(HELPER_DIR)/directsinc.o $(HELPER_DIR)/sincutil.o $(EXAMPLE_DIR)/example*d $(TEST_DIR)/test*d

