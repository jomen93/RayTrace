#!/bin/bash
# Build script for C++ Ray Tracer Module

set -e

echo "========================================"
echo "Building C++ Ray Tracer Module"
echo "========================================"

# Detect OS
UNAME_S=$(uname -s)
echo "Platform: $UNAME_S"

# Python
PYTHON=${PYTHON:-python3}
echo "Python: $PYTHON"

# Check pybind11
if ! $PYTHON -c "import pybind11" 2>/dev/null; then
    echo "Error: pybind11 not installed. Run: pip install pybind11"
    exit 1
fi
echo "✓ pybind11 found"

# Detect compiler
if [ "$UNAME_S" = "Darwin" ]; then
    CXX=clang++
    echo "✓ Compiler: clang++ (macOS)"
else
    if command -v g++ &> /dev/null; then
        CXX=g++
    elif command -v clang++ &> /dev/null; then
        CXX=clang++
    else
        echo "Error: No C++ compiler found"
        exit 1
    fi
    echo "✓ Compiler: $CXX"
fi

$CXX --version | head -1

# Base flags
CXXFLAGS="-O3 -std=c++14 -fPIC -Wall -Wextra -Wno-unused-parameter"

# Python includes
PYTHON_INCLUDES=$($PYTHON -m pybind11 --includes 2>/dev/null)
if [ -z "$PYTHON_INCLUDES" ]; then
    PYBIND_INC=$($PYTHON -c "import pybind11; print(pybind11.get_include())")
    PYTHON_INCLUDES="-I$PYBIND_INC"
    PYTHON_INCLUDES="$PYTHON_INCLUDES $($PYTHON-config --includes 2>/dev/null || echo "")"
fi
echo "✓ Python includes configured"

# Platform-specific flags
if [ "$UNAME_S" = "Darwin" ]; then
    # macOS
    LDFLAGS="-shared -undefined dynamic_lookup"
    
    # Check for OpenMP
    if brew --prefix libomp &> /dev/null; then
        LIBOMP_PREFIX=$(brew --prefix libomp)
        CXXFLAGS="$CXXFLAGS -Xpreprocessor -fopenmp -I$LIBOMP_PREFIX/include"
        LDFLAGS="$LDFLAGS -L$LIBOMP_PREFIX/lib -lomp"
        echo "✓ OpenMP: $LIBOMP_PREFIX"
    else
        echo "⚠ OpenMP not found (optional). Install: brew install libomp"
    fi
else
    # Linux
    CXXFLAGS="$CXXFLAGS -fopenmp"
    LDFLAGS="-shared -fopenmp"
    echo "✓ OpenMP enabled (Linux)"
fi

# Build
echo ""
echo "Compiling..."
cd cpp

# Clean old files
rm -f *.o *.so

echo "  → raytracer.cpp"
$CXX $CXXFLAGS $PYTHON_INCLUDES -c raytracer.cpp -o raytracer.o

echo "  → binding.cpp"  
$CXX $CXXFLAGS $PYTHON_INCLUDES -c binding.cpp -o binding.o

echo "  → Linking..."
$CXX $LDFLAGS -o raytracer_cpp.so raytracer.o binding.o

# Copy to parent
cp raytracer_cpp.so ..
echo "✓ Created: ../raytracer_cpp.so"

cd ..

echo ""
echo "========================================"
echo "Build successful!"
echo "========================================"
echo "Run: python main_cpp.py"
