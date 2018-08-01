
# Makefile for sinc transform programs

CXX=g++
CURRENT=/Users/hannah/Documents/Flatiron18/sinctransform
CURRENT_FLAG=-I$(CURRENT)
FFTW=/usr/local/lib
FFTW_FLAG=-L$(FFTW)
FINUFFT=/Users/hannah/Documents/Flatiron18/newfinufft/finufft
#/Users/hannah/Documents/Summer2017/Flatiron/fi2/finufft
FINUFFT_LIB_PATH=/lib-static/libfinufft.a
#/lib/libfinufft.a
FINUFFT_LIB=$(FINUFFT)$(FINUFFT_LIB_PATH)
FINUFFT_HEADER_FLAG=-I$(FINUFFT)/src/
FINUFFT_HEADER=$(FINUFFT)/src/finufft.h
FINUFFT_FLAGS=-lfftw3 -lm
FLAGS=-std=c++11 -g -Wall 
EXAMPLE_DIR=examples
TEST_DIR=tests
HELPER_DIR=helper
SOURCE_DIR=src

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

