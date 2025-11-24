# CLAUDE.md - Hessian Cosmic Structure Analysis

## Project Overview

**Hessian** is a scientific computing tool written in Fortran 90/95 for analyzing 3D cosmic structures in astrophysical simulations. It identifies and classifies large-scale structures (clusters, filaments, walls, and voids) in particle distributions using Hessian matrix eigenvalue analysis.

### Repository Information
- **Repository**: csabiu/hessian
- **Language**: Fortran 90/95
- **Primary Use**: Astrophysical/cosmological data analysis
- **Parallelization**: OpenMP

## Codebase Structure

```
hessian/
├── hessian.f90          # Main program (2596 lines)
├── kdtree2.F90          # KD-tree spatial indexing library
├── make.sh              # Build script
└── tmp.den              # Sample/test density data file
```

### File Descriptions

#### hessian.f90
The main program file containing:
- **Main program**: `density_histogram` (lines 1-393)
- **Internal subroutines** (lines 394-2596):
  - `calculate_density()`: Bins particles into a 3D grid
  - `read_data()`: Reads particle data from input files
  - `gaussian_smooth()`: Applies Gaussian smoothing to density fields
  - `gaussian_smooth2()`: Alternative smoothing implementation
  - `hessian_at_point()`: Computes Hessian matrix at a grid point
  - `heaviside()`: Heaviside step function
  - Various utility functions for structure analysis
  - RBF (Radial Basis Function) interpolation routines
  - Linear algebra routines (BLAS/LAPACK-style)

#### kdtree2.F90
Third-party library for KD-tree spatial indexing:
- **Author**: Matthew Kennel, Institute for Nonlinear Science (2004)
- **License**: Academic Free License version 1.1
- **Purpose**: Efficient nearest-neighbor searches in 3D space
- **Key modules**:
  - `kdtree2_precision_module`: Precision configuration
  - `kdtree2_priority_queue_module`: Priority queue for search results
  - Main KD-tree implementation

#### make.sh
Simple build script:
```bash
gfortran kdtree2.F90 hessian.f90 -o hessian -ffree-line-length-0 -llapack -fopenmp -O3
```

## Scientific Methodology

### Algorithm Overview

1. **Data Input**: Read particle positions (x,y,z) and masses from text file
2. **Density Field Construction**:
   - Create 3D grid based on user-specified bin size
   - Use KD-tree for efficient nearest-neighbor searches
   - Interpolate particle data onto grid using RBF interpolation
   - Convert to log-density field
3. **Multi-scale Smoothing**:
   - Apply Gaussian smoothing at multiple scales: σ = bin_size × √2^(l-1)
   - Default: 8 smoothing scales (`nsmooth = 8`)
4. **Hessian Analysis**:
   - Compute Hessian matrix of smoothed density field at each grid point
   - Calculate eigenvalues (λ₁ ≤ λ₂ ≤ λ₃)
   - Scale by σ² (smoothing scale)
5. **Structure Classification**:
   - **Clusters**: All eigenvalues negative (matter converges in all directions)
   - **Filaments**: Two negative eigenvalues (matter flows along one axis)
   - **Walls**: One negative eigenvalue (matter confined to a plane)
   - **Voids**: All positive eigenvalues (matter diverges)
6. **Signature Calculation**:
   - `scluster`: Cluster signature strength
   - `sfil`: Filament signature strength
   - `swall`: Wall signature strength
7. **Threshold Determination**: Find optimal thresholds by maximizing mass in each structure type
8. **Output**: Write classification and signature fields to text files

### Key Formulas

**Filament signature** (hessian.f90:199):
```fortran
sfil = (λ₂²/|λ₁|) × (1 - |λ₃/λ₁|) × H(1 - |λ₃/λ₁|) × H(-λ₂) × H(-λ₁)
```

**Wall signature** (hessian.f90:200):
```fortran
swall = |λ₁| × (1 - |λ₃/λ₁|) × (1 - |λ₂/λ₁|) × H(...) × H(...)
```

**Cluster signature** (hessian.f90:198):
```fortran
scluster = (λ₃²/|λ₁|) × H(-λ₃) × H(-λ₂) × H(-λ₁)
```

Where H(x) is the Heaviside function.

## Building and Running

### Dependencies
- **Compiler**: gfortran (GNU Fortran compiler)
- **Libraries**:
  - LAPACK (Linear Algebra PACKage) - for eigenvalue computation
  - OpenMP - for parallel processing
- **System**: Linux (tested on Linux 4.4.0)

### Build Instructions

```bash
# Make build script executable (if needed)
chmod +x make.sh

# Build the program
./make.sh

# This produces the 'hessian' executable
```

**Compiler flags explained**:
- `-ffree-line-length-0`: Allow unlimited line length (Fortran free-form)
- `-llapack`: Link against LAPACK library
- `-fopenmp`: Enable OpenMP parallelization
- `-O3`: Maximum optimization level

### Running the Program

**Usage**:
```bash
./hessian <filename> <bin_size> <calc_density>
```

**Arguments**:
1. `filename`: Input file path containing particle data
2. `bin_size`: Grid resolution (physical units per bin)
3. `calc_density`: Boolean flag (`.true.` or `.false.`)
   - `.true.`: Calculate density from scratch using binning
   - `.false.`: Use RBF interpolation with KD-tree (default, recommended)

**Example**:
```bash
./hessian tmp.den 1.1 .false.
```

### Input File Format

Text file with one particle per line:
```
x1 y1 z1 mass1
x2 y2 z2 mass2
...
```

Each line contains 4 space-separated values: x, y, z coordinates and mass.

### Output Files

The program generates several output files:

| File | Contents | Format |
|------|----------|--------|
| `density.txt` | Log-density field | `i j k density(i,j,k)` |
| `smoothed_density_N.txt` | Smoothed density at scale N | `i j k smoothed_density(i,j,k)` |
| `eigenvalues.txt` | Hessian eigenvalues (optional) | `i j k λ₁ λ₂ λ₃` |
| `sfil.txt` | Filament signature strength | `i j k sfil(i,j,k)` |
| `swall.txt` | Wall signature strength | `i j k swall(i,j,k)` |
| `scluster.txt` | Cluster signature strength | `i j k scluster(i,j,k)` |
| `structure.txt` | Final structure classification | `i j k type` |

**Structure type codes** in `structure.txt`:
- `0`: Void
- `1`: Wall
- `2`: Filament
- `3`: Cluster

## Code Organization and Key Sections

### Main Program Flow (hessian.f90:1-393)

1. **Initialization** (lines 27-52):
   - Print OpenMP thread count
   - Parse command-line arguments
   - Set flags: `write_eig`, `write_density`, `calc_density`

2. **Data Loading** (lines 55-78):
   - Count particles in input file
   - Allocate arrays
   - Call `read_data()`

3. **Grid Setup** (lines 82-113):
   - Determine bounding box
   - Calculate grid dimensions (nx, ny, nz)
   - Allocate density, eigenvalue, and signature arrays

4. **Density Field Construction** (lines 117-164):
   - **If calc_density = .true.**: Simple binning via `calculate_density()`
   - **If calc_density = .false.** (default):
     - Build KD-tree from particle data
     - Parallel loop over grid points (OpenMP)
     - Find 10 nearest neighbors (`nn2=10`)
     - Compute RBF weights
     - Interpolate log(mass) to grid

5. **Multi-scale Analysis** (lines 169-236):
   - Loop over smoothing scales (1 to `nsmooth`)
   - Apply Gaussian smoothing
   - **Parallel Hessian computation** (lines 184-219):
     - Compute Hessian at each interior grid point
     - Solve for eigenvalues using LAPACK's `dsyev`
     - Calculate structure signatures (sfil, swall, scluster)
   - Optionally write smoothed density to file

6. **Structure Classification** (lines 240-309):
   - Find maximum signatures across all scales
   - Determine thresholds by maximizing mass in each structure type
   - Classify each grid point

7. **Output** (lines 311-393):
   - Write eigenvalues (if `write_eig = .true.`)
   - Write signature fields (sfil.txt, swall.txt, scluster.txt)
   - Write final structure classification (structure.txt)
   - Write density field (if `write_density = .true.`)

### Important Subroutines

#### read_data() (line 417)
- Reads particle positions and masses from input file
- Populates `x(3, count)` and `mass(count)` arrays

#### gaussian_smooth() (line 436)
- Applies 3D Gaussian smoothing to density field
- Uses direct convolution (not FFT)
- Handles periodic boundary conditions
- **Parameters**: input field, output field, grid dimensions, sigma

#### hessian_at_point() (line 538)
- Computes 3×3 Hessian matrix from 3×3×3 local density field
- Uses finite difference approximation
- **Input**: 3×3×3 density cube, bin_size
- **Output**: 3×3 Hessian matrix

#### RBF Functions (lines 2026-2303)
- `phi1` to `phi5`: Different radial basis functions
- Default uses `phi5` (thin-plate spline variant)
- `rbf_weight()`: Compute RBF interpolation weights
- `rbf_interp_nd()`: Perform RBF interpolation

### Parallelization

The code uses OpenMP for parallel processing in two main sections:

1. **Density interpolation** (lines 134-156):
   ```fortran
   !$OMP PARALLEL DO PRIVATE(i,j,k,gridxyz,wgt,resultsb) SHARED(density,tree)
   ```
   - Parallelizes over grid points
   - Each thread has private KD-tree search results

2. **Hessian computation** (lines 184-219):
   ```fortran
   !$OMP PARALLEL DO PRIVATE(i,j,k,hessian,eigenvalues,work,info,field) SHARED(...)
   ```
   - Parallelizes eigenvalue computation
   - Each thread has private work arrays for LAPACK

**Performance notes**:
- Thread count is auto-detected by OpenMP
- Reported at program start
- Both sections are embarrassingly parallel (no dependencies)

## Development Guidelines for AI Assistants

### Code Style and Conventions

1. **Fortran Style**:
   - Free-form Fortran 90/95 syntax
   - 4-space indentation
   - Array indexing starts at 1 (Fortran convention)
   - Implicit typing disabled (`implicit none`)
   - Double precision throughout (`dp` kind parameter)

2. **Naming Conventions**:
   - Lowercase for all variables and subroutines
   - Underscores for multi-word names (snake_case)
   - Physics-related names: `density`, `eigenvalue1`, `hessian`, `sfil`, `swall`
   - Grid indices: `i, j, k` for spatial dimensions; `l` for smoothing scale

3. **Array Conventions**:
   - 3D spatial arrays: `array(nx, ny, nz)` or `array(nx, ny, nz, nsmooth)`
   - Particle data: `x(3, count)` - first index is dimension
   - Fortran column-major order (differs from C/Python)

4. **Commenting**:
   - Fortran comments start with `!`
   - Minimal inline comments (code should be self-documenting)
   - Important physics formulas are commented out for reference

### Common Gotchas

1. **Array Indexing**:
   - Fortran arrays are 1-indexed, not 0-indexed
   - Array slices are inclusive: `array(i-1:i+1)` includes all three elements

2. **OpenMP Thread Safety**:
   - Ensure all loop variables and temporaries are private
   - Shared arrays must not have race conditions
   - KD-tree searches are read-only, so safe to share tree pointer

3. **Precision**:
   - Use `_dp` suffix for literal constants: `1.0_dp`, not `1.0`
   - Mixing precisions can cause subtle bugs

4. **LAPACK**:
   - `dsyev` destroys input matrix (hessian)
   - Work array size must be correct (here: 3×3 = 9)
   - `info` parameter indicates success/failure

5. **Boundary Conditions**:
   - Hessian computed only for interior points: `i=2,nx-1`
   - Boundary points are not analyzed (missing neighbors)

6. **File I/O**:
   - Unit numbers (e.g., `10`, `20`) must be unique
   - Always close files after writing
   - Use `status='replace'` to overwrite existing output

### Testing and Debugging

1. **Test Data**:
   - `tmp.den` is included as test/example data
   - Contains particle distribution from cosmological simulation

2. **Validation**:
   - Check that mass is conserved: `sum(10**density)` should normalize to 1
   - Eigenvalues should be sorted: λ₁ ≤ λ₂ ≤ λ₃
   - Structure classifications should be mutually exclusive

3. **Performance**:
   - Most time spent in KD-tree searches and eigenvalue computation
   - Increasing `nn2` (nearest neighbors) improves accuracy but slows down
   - Grid resolution (bin_size) has cubic impact on memory and runtime

4. **Common Issues**:
   - **Segmentation fault**: Usually array allocation failure (too large grid)
   - **LAPACK error**: Check `info` parameter from `dsyev`
   - **Poor classification**: Try adjusting smoothing scales or bin_size

### Making Changes

#### Adding New Features

1. **New Structure Types**:
   - Add signature calculation in Hessian analysis loop (around line 198)
   - Update threshold determination logic (around line 275)
   - Add new output file writing section

2. **Different Smoothing**:
   - Modify `gaussian_smooth()` or add new smoothing subroutine
   - Note: FFT-based smoothing is commented out but not implemented

3. **Alternative RBF**:
   - Currently uses `phi5` (thin-plate spline)
   - Can switch to `phi1`-`phi4` by changing function pointer

#### Bug Fixes

1. **Recent Fixes** (from git log):
   - `f7cc080`: "bugs, updates and good stuff"
   - `cf6a6db`: "BIG bug fixed, field was inverted"
   - `9217c09`: "omp speedup"

2. **Known Issues**:
   - No comprehensive error handling for file I/O
   - No validation of command-line arguments
   - Memory not explicitly deallocated (relies on program exit)

#### Extending Functionality

**Safe to modify**:
- Smoothing parameters (`nsmooth`, sigma calculation)
- Output formatting
- Structure signature formulas
- Threshold determination method

**Modify with caution**:
- KD-tree usage (well-tested library)
- LAPACK calls (easy to get wrong)
- OpenMP directives (race conditions)
- RBF interpolation (numerically sensitive)

**Do not modify**:
- kdtree2.F90 (third-party library)
- Linear algebra subroutines (BLAS/LAPACK implementations)

## Git Workflow

### Branch Strategy

- **Main development**: Direct commits to main branch (small project)
- **Claude branches**: AI-assisted work uses branches like `claude/claude-md-*`
- **Commit style**: Brief, informal messages describing changes

### Current State

- Working directory: Clean (no uncommitted changes)
- Current branch: `claude/claude-md-mice5iy4chzt5o8k-01F27aoviwfhy5vLU5KcfqPG`
- Recent commits:
  - `f7cc080`: bugs, updates and good stuff
  - `cf6a6db`: BIG bug fixed, field was inverted
  - `9217c09`: omp speedup
  - `d3c03a1`: init

### Making Commits

1. **Before committing**:
   - Test with `./make.sh && ./hessian tmp.den 1.1 .false.`
   - Verify output files are reasonable

2. **Commit message style**:
   - Brief description of what changed
   - No need for formal conventional commits
   - Examples: "fix eigenvalue sorting bug", "add new output file"

3. **Push to remote**:
   ```bash
   git push -u origin <branch-name>
   ```
   - Use full branch name including `claude/` prefix
   - Retry with exponential backoff if network issues

## Performance Characteristics

### Computational Complexity

- **KD-tree construction**: O(N log N) where N = particle count
- **Density interpolation**: O(nx × ny × nz × nn2 × log N)
- **Gaussian smoothing**: O(nx × ny × nz × kernel_size³)
- **Hessian eigenvalues**: O(nx × ny × nz × nsmooth)

### Memory Requirements

Main arrays (all double precision, 8 bytes):
- Particle data: `x(3,N) + mass(N) + wgt(N)` = 40N bytes
- Density fields: `density + smoothed_density` = 2 × 8 × nx × ny × nz bytes
- Eigenvalues: 3 × 8 × nx × ny × nz bytes
- Signatures: 3 × nsmooth × 8 × nx × ny × nz bytes
- **Total**: ~40N + (5 + 3×nsmooth) × 8 × nx × ny × nz bytes

**Example**:
- N = 10⁶ particles, 100³ grid, nsmooth = 8
- Memory ≈ 40 MB + 232 MB = ~270 MB

### Typical Runtime

On modern workstation (single-threaded equivalent):
- 10⁶ particles, 100³ grid: ~5-10 minutes
- 10⁷ particles, 200³ grid: ~1-2 hours

OpenMP scaling: Near-linear up to ~8-16 threads (depends on grid size)

## Scientific Context

### Cosmic Web

The universe's large-scale structure resembles a cosmic web with four main components:

1. **Clusters**: Dense nodes where galaxy clusters form (~3% of volume, ~30% of mass)
2. **Filaments**: Elongated structures connecting clusters (~10% of volume, ~40% of mass)
3. **Walls**: Sheet-like structures (~15% of volume, ~20% of mass)
4. **Voids**: Underdense regions (~70% of volume, ~10% of mass)

### Hessian Analysis Method

Based on eigenvalues of the Hessian matrix of the density field:
- Developed by Hahn et al. (2007), Forero-Romero et al. (2009)
- Also called "NEXUS" or "T-web" method variants
- Multi-scale approach handles hierarchical structure

### Applications

- Understanding galaxy formation in different environments
- Testing cosmological simulations
- Measuring large-scale structure statistics
- Tracing dark matter distribution

## References and Resources

### Key Papers

1. Hahn, O., et al. (2007): "Properties of Dark Matter Haloes in Clusters, Filaments, Sheets and Voids"
2. Forero-Romero, J., et al. (2009): "A Dynamical Classification of the Cosmic Web"
3. Cautun, M., et al. (2013): "NEXUS: Tracing the Cosmic Web Connection"

### External Dependencies

- **LAPACK**: [www.netlib.org/lapack](http://www.netlib.org/lapack)
- **KD-tree library**: Matthew Kennel's implementation (included)
- **gfortran**: [gcc.gnu.org/fortran](https://gcc.gnu.org/fortran/)

### Related Tools

- NEXUS/NEXUS+: More sophisticated implementations
- DisPerSE: Morse theory-based structure finder
- SpineWeb: Alternative Hessian-based tool

## Quick Reference

### File Reading Quick Ref
```fortran
! Check lines in file
wc -l <filename>

! View subroutines
grep "^subroutine" hessian.f90

! Find OpenMP sections
grep "OMP PARALLEL" hessian.f90
```

### Key Variables
- `nx, ny, nz`: Grid dimensions
- `bin_size`: Grid resolution (physical units)
- `count`: Number of particles
- `x(3,count)`: Particle positions
- `mass(count)`: Particle masses (log10 after reading)
- `density(nx,ny,nz)`: Log-density field
- `smoothed_density(nx,ny,nz)`: Smoothed log-density
- `eigenvalue1/2/3(nx,ny,nz)`: Hessian eigenvalues
- `sfil/swall/scluster(nx,ny,nz,nsmooth)`: Structure signatures
- `structure(nx,ny,nz)`: Final classification (0=void, 1=wall, 2=fil, 3=cluster)
- `sigma`: Smoothing scale (physical units)
- `nsmooth`: Number of smoothing scales (default 8)
- `nn2`: Number of nearest neighbors for interpolation (default 10)

### Important Line Numbers
- Main program start: line 2
- Command-line parsing: lines 38-52
- KD-tree creation: line 127
- Density interpolation: lines 134-156
- Multi-scale loop: lines 169-236
- Hessian computation: lines 184-219
- Structure classification: lines 295-309
- Output writing: lines 311-393
- Subroutines start: line 394

---

**Last Updated**: 2025-11-24
**For**: AI assistants working with csabiu/hessian repository
