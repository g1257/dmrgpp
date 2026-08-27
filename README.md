# Quick Start

## Description

DMRG++ is a free and open source implementation of the
density matrix renormalization group (DMRG) algorithm.
You are welcomed to use it and publish data
 obtained with DMRG++. If you do,
*please cite this work* (see [`CITATION.cff`](CITATION.cff)).

## License and Disclaimers

The full software license for DMRG++ v7.00
can be found in file LICENSE in the root directory of the code,
along with disclaimers.

## Building and Running DMRG++

### Required Software

* GNU C++ or LLVM CLANG++ (C++17 is used)

* The BLAS and LAPACK library or equivalent

* HDF5 devel with C++ support

* boost-devel (boost-spirit) for Ainur
Only headers files are used; boost runtime is not used.

* cmake and dependencies

* GSL (GNU Scientific Library)

* Kokkos and Kokkos-Kernels [automatic download possible]

* Catch2 for testing [automatic download possible]

### Downloading DMRG++
Create a directory somewhere and cd to it.

<pre>
git clone https://github.com/dmrgpp-project/dmrgpp.git
cd dmrgpp/
</pre>
Pull requests should be opened against the `master` branch.

### Configure

Use the following command to configure DMRG++
```
cmake -B builddir [<options...>]
```
`-B builddir` creates a build directory named `builddir` (you can choose a
different name if you prefer).

The `[<options...>]` part is where you specify other configuration options.

#### Common CMake Options
These options are generally useful for any CMake project:

* `-DCMAKE_CXX_COMPILER=<compiler>`: Specifies the full path to the C++
  compiler.

  Example: `-DCMAKE_CXX_COMPILER=clang++`

* `-DCMAKE_CXX_STANDARD=<standard>`: Sets the C++ standard. The default is `17`.

  Example: `-DCMAKE_CXX_STANDARD=20`

* `-DCMAKE_BUILD_TYPE=<type>`: Controls optimization level and debugging
  information. Common options are `Debug`, `Release`, `RelWithDebInfo`
  (default), and `MinSizeRel`.

* `-G<generator>`: Specifies the generator to create a native build system
  (e.g., `"Unix Makefiles"` for standard UNIX makefiles, `"Ninja"` for
  `build.ninja`  or `"Visual Studio 17 2022"` for Visual Studio.

  Example: `-GNinja`

#### External dependencies (Third-Party Libraries)
To control where the project should find external dependencies installed on
your system, you may pass `-D<PackageName>_ROOT=/path/to/package/install` to
find a specific `<PackageName>` installed at `/path/to/package/install`.

External dependencies are: `BLAS`, `Boost`, `GSL`, `HDF5`, `Kokkos`, `KokkosKernels`, `LAPACK` and `MPI`.

Example: `-DBoost_ROOT=/path/to/boost`

> **Note on Catch2 and Kokkos:**
> We use the Catch2 library as our unit test framework and Kokkos as required dependency.
> By default, our CMake is set up to use `FetchContent` to download these libraries if they aren't found on
> your system.
>
> To force the use of a locally installed library and disable the automatic
> download fallback, use the following options (where Catch2 serves as an example):
>
> ```bash
> -DCMAKE_REQUIRE_FIND_PACKAGE_Catch2=ON -DCatch2_ROOT=/path/to/catch2
> ```
>
> By default, `CMake` tries to check externally dependencies that have been installed via FetchContent for updates.
> This can be avoided by configuring with
> ```bash
> -DFETCHCONTENT_UPDATES_DISCONNECTED=ON
> ```
> which might be especially useful when there is limited access to the internet.

### Build and Test
After the configuration step succeeded, build using
```
cmake --build builddir
```
and test with
```
ctest --test-dir builddir --output-on-failure
```

### Install
After building, install with
```
cmake --install builddir --prefix installdir
```

### Running DMRG++

Assuming you are in your build directory
copy input2.ain to it and then you may run with
```
./dmrg/dmrg -f input.ain
```

You will now have two files a data2.hdf5 and an ASCII file runForinput2.cout.
The name data2 is obtained from the corresponding label in the input file,
in this case input2.ain. Normally the code writes stdout to
runForinput2.cout for an input called input2.inp, and stderr to the
terminal. If you would like to override the default inferred
name runForinput2.cout you may use
<code>./dmrg/dmrg -f ../TestSuite/inputs/input2.ain -l myoutputfile</code>
If you would like stdout be written to the terminal say -l -

### Profiling

DMRG++ uses Kokkos for its most compute intensive components. Thus, it supports profiling through [Kokkos Tools](https://github.com/kokkos/kokkos-tools).
To use, e.g., the [Space Time Stack](https://github.com/kokkos/kokkos-tools/wiki/Space-Time-Stack) tool set the `KOKKOS_TOOLS_LIBS` environment variable
to the corresponding library path
```
export KOKKOS_TOOLS_LIBS=${YOUR_KOKKOS_TOOLS_INSTALL_DIR}/lib/libkp_space_time_stack.so
```
before running `dmrg`.
