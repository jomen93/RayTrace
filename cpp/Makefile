# Makefile for Ray Tracer C++ Module

# Detect OS
UNAME_S := $(shell uname -s)

# Python configuration
PYTHON := python3

# Compiler
ifeq ($(UNAME_S),Darwin)
    CXX := clang++
else
    CXX := g++
endif

# Base flags
CXXFLAGS := -O3 -std=c++14 -fPIC -Wall -Wextra -Wno-unused-parameter

# Python includes
PYTHON_INCLUDES := $(shell $(PYTHON) -m pybind11 --includes 2>/dev/null)
ifeq ($(PYTHON_INCLUDES),)
    PYTHON_INCLUDES := -I$(shell $(PYTHON) -c "import pybind11; print(pybind11.get_include())")
    PYTHON_INCLUDES += $(shell $(PYTHON)-config --includes 2>/dev/null)
endif

# Platform-specific settings
ifeq ($(UNAME_S),Darwin)
    # macOS with clang
    LDFLAGS := -shared -undefined dynamic_lookup
    
    # OpenMP on macOS - disabled for now due to compatibility issues
    # Uncomment below if you have libomp installed and want to use it
    # ifneq ($(shell brew --prefix libomp 2>/dev/null),)
    #     LIBOMP_PREFIX := $(shell brew --prefix libomp)
    #     CXXFLAGS += -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
    #     LDFLAGS += -L$(LIBOMP_PREFIX)/lib -lomp
    #     $(info Using OpenMP from $(LIBOMP_PREFIX))
    # endif
else
    # Linux with g++
    CXXFLAGS += -fopenmp
    LDFLAGS := -shared -fopenmp
endif

# Output
TARGET := raytracer_cpp.so
SOURCES := raytracer.cpp binding.cpp
OBJECTS := $(SOURCES:.cpp=.o)

.PHONY: all clean check

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^
	@echo "✓ Built: $@"

%.o: %.cpp raytracer.h
	$(CXX) $(CXXFLAGS) $(PYTHON_INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

check:
	@echo "OS: $(UNAME_S)"
	@echo "Compiler: $(CXX)"
	@echo "Python includes: $(PYTHON_INCLUDES)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"

install: $(TARGET)
	cp $(TARGET) ..
	@echo "✓ Installed to .."
