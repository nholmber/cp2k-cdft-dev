# How to compile the CP2K code

##  1. Acquire the code:

See  https://www.cp2k.org/download

For users, the preferred method is to download a release.
For developers, the preferred method is to download it from Git.

## 2. Install Prerequisites

Sub-points here discuss prerequisites needed to build CP2K. Most of these can be conveniently installed via the [toolchain script](./tools/toolchain). Copies of the recommended versions of 3rd party software can be downloaded from https://www.cp2k.org/static/downloads/.

### 2a. GNU make (required, build system)

GNU make should be on your system (gmake or make on linux) and used for the build, go to https://www.gnu.org/software/make/make.html download from https://ftp.gnu.org/pub/gnu/make/ also Python (2.x) is required for building.

### 2b. Python (required, build system)
Python 2.x is needed to run the dependency generator. On most system Python is already installed. For more information visit: https://www.python.org/

### 2c. Fortran and C Compiler (required, build system)
A Fortran 2003 compiler and matching C compiler should be installed on your system. We have good experience with gcc/gfortran (gcc >=4.6 works, later version recommended). Be aware that some compilers have bugs that might cause them to fail (internal compiler errors, segfaults) or, worse, yield a mis-compiled CP2K. Report bugs to compiler vendors; they (and we) have an interest in fixing them. Always run a `make -j test` (See point 5.) after compilation to identify these problems.

### 2d. BLAS and LAPACK (required, base functionality)
BLAS and LAPACK should be installed.  Using vendor-provided libraries can make a very significant difference (up to 100%, e.g., ACML, MKL, ESSL), not all optimized libraries are bug free. Use the latest versions available, use the interfaces matching your compiler, and download all patches!

  * The canonical BLAS and LAPACK can be obtained from the Netlib repository:
    * http://www.netlib.org/blas/
    * http://www.netlib.org/lapack/ and see also
    * http://www.netlib.org/lapack-dev/
  * Open fast alternatives, include:
    * http://www.openblas.net/
    * http://math-atlas.sourceforge.net/
    * https://www.tacc.utexas.edu/research-development/tacc-software/gotoblas2

If compiling with OpenMP support then it is recommended to use a non-threaded version of BLAS. In particular if compiling with MKL and using OpenMP you must define `-D__MKL` to ensure the code is thread-safe. MKL with multiple OpenMP threads in CP2K requires that CP2K was compiled with the Intel compiler (and `-D__INTEL_COMPILER` is defined if explicit pre-processing is performed using cpp instead of the compiler).

On the Mac, BLAS and LAPACK may be provided by Apple's Accelerate framework. If using this framework, `-D__ACCELERATE` must be defined to account for some interface incompatibilities between Accelerate and reference BLAS/LAPACK.

When building on/for Windows using the Minimalist GNU for Windows (MinGW) environment, you must set `-D__MINGW`,  `-D__NO_STATM_ACCESS` and `-D__NO_IPI_DRIVER` to avoid undefined references during linking, respectively errors while printing the statistics.

## 2e. MPI and SCALAPACK (optional, required for MPI parallel builds)
MPI (version 2) and SCALAPACK are needed for parallel code. (Use the latest versions available and download all patches!).

:warning: Note that your MPI installation must match the used Fortran compiler. If your computing platform does not provide MPI, there are several freely available alternatives:

  * MPICH2 MPI:  http://www-unix.mcs.anl.gov/mpi/mpich/
  * OpenMPI MPI: http://www.open-mpi.org/
  * ScaLAPACK:
    * http://www.netlib.org/scalapack/
    * http://www.netlib.org/lapack-dev/
    * ScaLAPACK can be part of ACML or cluster MKL. These libraries are recommended if available.
    * Recently a [ScaLAPACK installer](http://www.netlib.org/scalapack/scalapack_installer.tgz) has been added that simplifies the installation.

CP2K assumes that the MPI library implements MPI version 3. If you have an older version of MPI (e.g. MPI 2.0) available you must define `-D__MPI_VERSION=2` in the arch file.

## 2f. FFTW (optional, improved performance of FFTs)
FFTW can be used to improve FFT speed on a wide range of architectures. It is strongly recommended to install and use FFTW3. The current version of CP2K works with FFTW 3.X (use `-D__FFTW3`). It can be downloaded from http://www.fftw.org/

:warning: Note that FFTW must know the Fortran compiler you will use in order to install properly (e.g., `export F77=gfortran` before configure if you intend to use gfortran).

:warning: Note that on machines and compilers which support SSE you can configure FFTW3 with `--enable-sse2`. Compilers/systems that do not align memory (NAG f95, Intel IA32/gfortran) should either not use `--enable-sse2` or otherwise set the define `-D__FFTW3_UNALIGNED` in the arch file. When building an OpenMP parallel version of CP2K (ssmp or psmp), the FFTW3 threading library libfftw3_threads (or libfftw3_omp) is required.

## 2g. LIBINT (optional, enables methods including HF exchange)
Hartree-Fock exchange (optional, use `-D__LIBINT`) requires the libint package to be installed.
  * Download from http://sourceforge.net/projects/libint/files/v1-releases/libint-1.1.4.tar.gz/download.
  * Additional information can be found in [README_LIBINT](./tools/hfx_tools/libint_tools/README_LIBINT).
  * Tested against libinit-1.1.4 and currently hardcoded to the default angular momentum LIBINT_MAX_AM 5.
  * Use `-D__LIBINT_MAX_AM` and `-D__LIBDERIV_MAX_AM1` to match the values in `include/libint/libint.h`.
  * `-D__MAX_CONTR=4` (default=2) can be used to compile efficient contraction kernels up to l=4, but the build time will increase accordingly.
  * :warning: Do **NOT** use libinit-1.1.3, which was buggy.

### 2h. libsmm (optional, improved performance for matrix multiplication)
  * A library for small matrix multiplies can be built from the included source (see tools/build_libsmm/README).  Usually only the double precision real and perhaps complex is needed.  Link to the generated libraries. For a couple of architectures prebuilt libsmm are available at https://www.cp2k.org/static/downloads/libsmm/.
  * Add `-D__HAS_smm_dnn` to the defines to make the code use the double precision real library.  Similarly use `-D__HAS_smm_snn` for single precision real and `-D__HAS_smm_znn` / `-D__HAS_smm_cnn` for double / single precision complex.
  * Add `-D__HAS_smm_vec` to enable the new vectorized interfaces of libsmm.

### 2i. libxsmm (optional, improved performance for matrix multiplication)
  * A library for matrix operations and deep learning primitives: https://github.com/hfp/libxsmm/
  * Add `-D__LIBXSMM` to enable it (with suitable include and library paths)

### 2j. CUDA (optional, improved performance on GPU systems)
  * `-D__ACC` needed to enable accelerator support.
  * Use the `-D__DBCSR_ACC` to enable accelerator support for matrix multiplications.
  * Add `-lcudart -lrt` to LIBS.
  * Use `-D__PW_CUDA` for CUDA support for PW (gather/scatter/fft) calculations.
  * CUFFT 7.0 has a known bug and is therefore disabled by default. NVidia's webpage list a patch (an upgraded version cufft i.e. >= 7.0.35) - use this together with `-D__HAS_PATCHED_CUFFT_70`.
  * Use `-D__CUDA_PROFILING` to turn on Nvidia Tools Extensions.
  * Link to a blas/scalapack library that accelerates large DGEMMs (e.g. libsci_acc)

### 2k. libxc (optional, wider choice of xc functionals)
  * The version 4.0.3 (or later) of libxc can be downloaded from http://www.tddft.org/programs/octopus/wiki/index.php/Libxc.
  * During the installation, the directories `$(LIBXC_DIR)/lib` and `$(LIBXC_DIR)/include` are created.
  * Add `-D__LIBXC` to DFLAGS, `-I$(LIBXC_DIR)/include` to FCFLAGS and `-L$(LIBXC_DIR)/lib -lxcf03 -lxc` to LIBS.
  * :warning: Note that the deprecated flags `-D__LIBXC2` and `-D__LIBXC3` are ignored.

### 2l. ELPA (optional, improved performance for diagonalization)
Library ELPA for the solution of the eigenvalue problem
  * ELPA replaces the ScaLapack SYEVD to improve the performance of the diagonalization
  * A version of ELPA can to be downloaded from http://elpa.rzg.mpg.de/software.
  * During the installation the libelpa.a (or libelpa_mt.a if omp active) is created.
  * Add `-D__ELPA=YYYYMM` to DFLAGS, where `YYYYMM` denotes the release date of the library.
  * Currently supported versions are: `201112`, `201308`, `201311`, `201406`, `201502`, `201505`, `201511`, `201605`, and `201611`.
  * Add `-I$(ELPA_INCLUDE_DIR)/modules` to FCFLAGS
  * Add `-I$(ELPA_INCLUDE_DIR)/elpa` to FCFLAGS
  * Add `-L$(ELPA_DIR)` to `LDFLAGS`
  * Add `-lelpa` to LIBS
  * For specific architectures it can be better to install specifically
    optimized kernels (see BG) and/or employ a higher optimization level to compile it.

### 2m. PEXSI (optional, low scaling SCF method)
The Pole EXpansion and Selected Inversion (PEXSI) method requires the PEXSI library and two dependencies (ParMETIS or PT-Scotch and SuperLU_DIST).
  * Download PEXSI (www.pexsi.org) and install it and its dependencies by following its README.md.
  * PEXSI versions 0.10.x have been tested with CP2K. Older versions are not supported.
  * PEXSI needs to be built with `make finstall`.

In the arch file of CP2K:
  * Add `-lpexsi_${SUFFIX} -llapack -lblas -lsuperlu_dist_3.3 -lparmetis -lmetis`, and their paths (with `-L$(LIB_DIR)`) to LIBS.
  * It is important that a copy of LAPACK and BLAS is placed before and after these libraries  (replace `-llapack` and `-lblas` with the optimized versions as needed).
  * In order to link in PT-Scotch instead of ParMETIS replace `-lparmetis -lmetis` with: `-lptscotchparmetis -lptscotch -lptscotcherr -lscotchmetis -lscotch -lscotcherr`
  * Add `-I$(PEXSI_DIR)/fortran/` to FCFLAGS.
  * Add `-D__LIBPEXSI` to DFLAGS.

Below are some additional hints that may help in the compilation process:
  * For building PT-Scotch, the flag `-DSCOTCH_METIS_PREFIX` in `Makefile.inc` must not be set and the flag `-DSCOTCH_PTHREAD` must be removed.
  * For building SuperLU_DIST with PT-Scotch, you must set the following in `make.inc`:

```
METISLIB = -lscotchmetis -lscotch -lscotcherr
PARMETISLIB = -lptscotchparmetis -lptscotch -lptscotcherr
```

### 2n. QUIP (optional, wider range of interaction potentials)
QUIP - QUantum mechanics and Interatomic Potentials Support for QUIP can be enabled via the flag `-D__QUIP`.

For more information see http://www.libatoms.org/ .

### 2o. PLUMED (optional, enables various enhanced sampling methods)
CP2K can be compiled with PLUMED 2.x (`-D__PLUMED2`).

See https://cp2k.org/howto:install_with_plumed for full instructions.

### 2p. spglib (optional, crystal symmetries tools)
A library for finding and handling crystal symmetries
  * The spglib can be downloaded from https://github.com/atztogo/spglib
  * For building CP2K with the spglib add `-D__SPGLIB` to DFLAGS

### 2q. JSON-Fortran (optional, required for SIRIUS)
JSON-Fortran is a Fortran 2008 JSON API.
  * The code is available at https://github.com/jacobwilliams/json-fortran
  * For building CP2K with JSON-Fortran add `-D__JSON` to DFLAGS.

### 2r. SIRIUS (optional, plane wave calculations)
SIRIUS is a domain specific library for electronic structure calculations.
  * The code is available at https://github.com/electronic-structure/SIRIUS
  * For building CP2K with SIRIUS add `-D__SIRIUS` to DFLAGS.
  * Furthermore, SIRIUS depends on JSON-Fortran.
  * See https://electronic-structure.github.io/SIRIUS/ for more information.

## 3. Compile

### 3a. ARCH files
The location of compiler and libraries needs to be specified. Examples for a number of common architectures examples can be found in [arch folder](./arch/). The names of these files match `architecture.version` e.g., [Linux-x86-64-gfortran.sopt](./arch/Linux-x86-64-gfortran.sopt). Alternatively https://dashboard.cp2k.org/ provides sample arch files as part of the testing reports (click on the status field, search for 'ARCH-file').
  * With -DNDEBUG assertions may be stripped ("compiled out").
  * NDEBUG is the ANSI-conforming symbol name (not __NDEBUG).
  * Regular release builds may carry assertions for safety.

Conventionally, there are six versions:

| Acronym |        Meaning          |         Recommended for            |
|---------|-------------------------|------------------------------------|
| sdbg    | serial                  | single core testing and debugging  |
| sopt    | serial                  | general single core usage          |
| ssmp    | parallel (only OpenMP)  | optimized, single node, multi core |
| pdbg    | parallel (only MPI)     | multinode testing and debugging    |
| popt    | parallel (only MPI)     | general usage, no threads          |
| psmp    | parallel (MPI + OpenMP) | general usage, threading might improve scalability and memory usage |

You'll need to modify one of these files to match your system's settings.

You can now build CP2K using these settings (where -j N allows for a parallel build using N processes):
```
> cd cp2k/makefiles
> make -j N ARCH=architecture VERSION=version
```
e.g.
```
> make -j N ARCH=Linux-x86-64-gfortran VERSION=sopt
```
as a short-cut, you can build several version of the code at once
```
> make -j N ARCH=Linux-x86-64-gfortran VERSION="sopt popt ssmp psmp"
```
An executable should appear in the `./exe/` folder.

All compiled files, libraries, executables, .. of all architectures and versions can be removed with
```
> make distclean
```
To remove only objects and mod files (i.e., keep exe) for a given ARCH/VERSION use, e.g.,
```
> make ARCH=Linux-x86-64-gfortran VERSION=sopt clean
```
to remove everything for a given ARCH/VERSION use, e.g.,
```
> make ARCH=Linux-x86-64-gfortran VERSION=sopt realclean
```

### 3b. Compilation Flags

The following flags should be present (or not) in the arch file, partially depending on installed libraries (see 2.)
  * `-D__parallel -D__SCALAPACK` parallel runs
  * `-D__LIBINT` use libint (needed for HF exchange)
  * `-D__LIBXC` use libxc
  * `-D__ELPA` use ELPA in place of SYEVD  to solve the eigenvalue problem
  * `-D__FFTW3` FFTW version 3 is recommended
  * `-D__PW_CUDA` CUDA FFT and associated gather/scatter on the GPU
  * `-D__MKL` link the MKL library for linear algebra and/or FFT

  * with `-D__GRID_CORE=X` (with X=1..6) specific optimized core routines can be selected.  Reasonable defaults are [provided](./src/grid/collocate_fast.f90) but trial-and-error might yield (a small ~10%) speedup.
  * with `-D__HAS_LIBGRID` (and `-L/path/to/libgrid.a` in LIBS) tuned versions of integrate and collocate routines can be [generated](./tools/autotune_grid/README).
  * `-D__PILAENV_BLOCKSIZE`: can be used to specify the blocksize (e.g. `-D__PILAENV_BLOCKSIZE=1024`), which is a hack to overwrite (if the linker allows this) the PILAENV function provided by Scalapack. This can lead to much improved PDGEMM performance. The optimal value depends on hardware (GPU?) and precise problem. Alternatively, Cray provides an environment variable to this effect (e.g. `export LIBSCI_ACC_PILAENV=4000`)
  * `-D__STATM_RESIDENT` or `-D__STATM_TOTAL` toggles memory usage reporting between resident memory and total memory
  * `-D__CRAY_PM_ACCEL_ENERGY` or `-D__CRAY_PM_ENERGY` switch on energy profiling on Cray systems
  * `-D__NO_ABORT` to avoid calling abort, but STOP instead (useful for coverage testing, and to avoid core dumps on some systems)

Features useful to deal with legacy systems
  * `-D__NO_MPI_THREAD_SUPPORT_CHECK`  - Workaround for MPI libraries that do not declare they are thread safe (funneled) but you want to use them with OpenMP code anyways.
  * `-D__HAS_NO_MPI_MOD` - workaround if mpi has been built for a different (version of the) Fortran compiler, rendering the MPI module unreadable (reverts to f77 style mpif.h includes)
  * `-D__NO_IPI_DRIVER` disables the socket interface in case of troubles compiling on systems that do not support POSIX sockets.
  * `-D__HAS_IEEE_EXCEPTIONS` disables trapping temporarily for libraries like scalapack.
  * The Makefile automatically compiles in the path to the data directory via the flag `-D__DATA_DIR`. If you want to compile in a different path, set the variable `DATA_DIR` in your arch-file.
  * `-D__HAS_NO_OMP_3` CP2K assumes that compilers support OpenMP 3.0. If this is not the case specify this flag to compile. Runtime performance will be poorer on low numbers of processors
  * `-D__HAS_NO_CUDA_STREAM_PRIORITIES` - Needed for CUDA sdk version < 5.5
  * `-D__NO_STATM_ACCESS` - Do not try to read from /proc/self/statm to get memory usage information. This is otherwise attempted on several. Linux-based architectures or using with the NAG, gfortran, compilers.
  * `-D__F2008` Allow for conformity check with the Fortran 2008 standard when using the GFortran compiler flag `-std=f2008`

## 4. If it doesn't work?
If things fail, take a break... go back to 2a (or skip to step 6).

## 5. Regtesting

If compilation works fine, it is recommended to test the generated binary, to exclude errors in libraries, or miscompilations, etc.
```
make -j ARCH=... VERSION=... test
```

should work if you can locally execute CP2K without the need for e.g. batch submission.

In the other case, you might need to configure the underlying testing script as described more systematically at https://www.cp2k.org/dev:regtesting

## 6. Talk to us
In any case please tell us your comments, praise, criticism, thanks,... see https://www.cp2k.org/

## 7. Manual
A reference manual of CP2K can be found on the web: https://manual.cp2k.org/ or can be generated using the cp2k executable, see https://manual.cp2k.org/trunk/generate_manual_howto.html

## 8. Happy computing!

 The CP2K team.
