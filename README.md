LAO - Linear Algebra Operations
===

Note this isn't industrial-grade, this is just a fun project.

Features:
- header-only C++20 linear algebra library.
- supports dense and sparse matrices, which can also be used to construct row and column vectors.
- arithmetic operators are implemented using expression templates.
- matrix size is determined and validated at compile-time.

The build script `scripts/build.sh` can be used to build the examples or unit tests. `scripts/run_tests.sh` runs the unit tests.