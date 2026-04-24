# spopt: Spatial Optimization for Regionalization, Facility Location, and Market Analysis

Implements spatial optimization algorithms across several problem
families: contiguity-constrained regionalization, discrete facility
location, market share analysis, and least-cost corridor and route
optimization over raster cost surfaces. Facility location problems also
accept user-supplied network travel-time matrices. Uses a 'Rust' backend
via 'extendr' for graph and routing algorithms, and the 'HiGHS' solver
via the 'highs' package for facility location mixed-integer programs.
Method-level references are provided in the documentation of the
individual functions.

## See also

Useful links:

- <https://walker-data.com/spopt-r/>

- <https://github.com/walkerke/spopt-r>

- Report bugs at <https://github.com/walkerke/spopt-r/issues>

## Author

**Maintainer**: Kyle Walker <kyle@walker-data.com>

Other contributors:

- PySAL Developers (Original Python spopt library) \[copyright holder\]
