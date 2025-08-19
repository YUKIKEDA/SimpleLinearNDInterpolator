# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SimpleLinearNDInterpolator is a C++ library for performing linear interpolation on scattered points in N-dimensional space using Delaunay triangulation. The library supports both scalar and vector-valued interpolation and includes nearest neighbor fallback for points outside the convex hull.

## Build and Test Commands

### Building the Project
```bash
# Configure build system
mkdir build
cd build
cmake ..

# Build the library and tests
cmake --build . --config Release

# For Windows with Visual Studio, you may need:
# cmake --build . --config Debug
```

### Running Tests
```bash
# Run all tests
ctest --output-on-failure

# Run specific test executable
./test_interpolator  # Linux/Mac
# or
.\test_interpolator.exe  # Windows
```

### Test Generation
The project includes Python scripts for generating SciPy compatibility tests:
```bash
cd scripts
python generate_scipy_test.py
```

## Architecture and Key Components

### Core Library Structure
- **Header**: `include/SimpleLinearNDInterpolator.h` - Complete API definition with extensive documentation
- **Implementation**: `src/SimpleLinearNDInterpolator.cpp` - Full implementation with detailed comments
- **Dependencies**: 
  - Qhull 2020.2 (computational geometry, located in `third_party/qhull-2020.2/`)
  - Eigen 3.4.0 (linear algebra, located in `third_party/eigen-3.4.0/`)
  - GoogleTest (testing framework, automatically downloaded)

### Key Classes and Methods
The main class is `SimpleLinearNDInterpolator` with two constructor overloads:
- Vector values: `SimpleLinearNDInterpolator(points, values)` where both are 2D vectors
- Scalar values: `SimpleLinearNDInterpolator(points, values)` where values is 1D vector

Main interpolation methods:
- `interpolate(query_point, use_nearest_neighbor_fallback=false)` - single point
- `interpolate(query_points, use_nearest_neighbor_fallback=false)` - multiple points

### Data Storage Format
The library uses "point-major" format internally:
- `points_[i][d]`: i-th point's d-th coordinate
- `values_[i][v]`: i-th point's v-th value component

### Algorithm Implementation
1. **1D Case**: Direct linear interpolation with sorted points (no triangulation)
2. **2D+ Case**: Delaunay triangulation using Qhull, then barycentric coordinate interpolation
3. **Fallback**: Nearest neighbor interpolation for points outside convex hull
4. **Linear Solver**: Uses Eigen's PartialPivLU for numerical stability

### Test Structure
Tests are located in `tests/` directory:
- `test_1d_interpolation.cpp` - 1D interpolation tests
- `test_nearest_neighbor_fallback.cpp` - nearest neighbor fallback tests  
- `generated_scipy_compat_test.cpp` - auto-generated SciPy compatibility tests

## Important Implementation Details

### Windows-Specific Build Configuration
The CMakeLists.txt includes Windows-specific settings:
- RuntimeLibrary set to "MultiThreadedDLL" for qhull compatibility
- Iterator debug level set to 0 to match qhull build
- Special debug flags to maintain debug info while avoiding conflicts

### Memory Management
- Uses `std::unique_ptr` for Qhull objects
- All containers are STL-based for automatic memory management
- No manual memory allocation/deallocation

### Numerical Considerations
- Uses epsilon tolerance (1e-12) for floating-point comparisons
- Eigen library provides numerical stability for linear equation solving
- Handles degenerate cases (duplicate points, collinear points) gracefully

### Performance Characteristics
- Construction: O(n log n) due to Delaunay triangulation
- Interpolation: O(log n) per query for simplex search
- 1D special case: O(n log n) for initial sort, then O(1) per query
- Nearest neighbor: O(n) linear search (optimization opportunity)

## Development Workflow

When adding new features or fixing bugs:
1. Modify the header file for API changes
2. Update the implementation in the .cpp file
3. Add corresponding tests in the appropriate test file
4. Run tests to ensure compatibility
5. Update documentation comments if needed

The codebase has extensive Japanese and English documentation in comments, making it well-suited for international collaboration.