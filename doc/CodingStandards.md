# Coding Standards for DMRG++

## Introduction
This document describes coding standards that are used in the **DMRG++**
project. Adhering to these standards ensures code consistency, readability, and
maintainability across the entire codebase.

---

## Source Code Formatting

All source code files must be automatically formatted using the specified
tools. We enforce these formatting rules using a continuous integration (CI)
build that checks for their successful application via `pre-commit`.

### C++ Source Code
All C++ source code must be formatted with
[clang-format](https://clang.llvm.org/docs/ClangFormat.html).
* Version of Truth: We use version **20.1.8** as the source of truth for
  the ClangFormat configuration (defined in the project's `.clang-format`
  file).

### CMake Files
All CMake files (`CMakeLists.txt`, `.cmake` files, etc.) must be formatted
with
[cmake-format](https://cmake-format.readthedocs.io/en/latest/cmake-format.html).

---

## Enforcement: Using `pre-commit`
We use [pre-commit](https://pre-commit.com/) to automatically check and apply
formatting rules before you push your code. This is the recommended and easiest
way to ensure your contributions meet the standards.

### How to Use pre-commit

1. **Installation:** Install `pre-commit` (e.g., via `pip install pre-commit`).

2. **Setup:** Run `pre-commit install` in the root of the DMRG++ repository.
   This installs the git hooks.

Once installed, every time you run `git commit`, `pre-commit` will
automatically execute `clang-format` and `cmake-format` (along with any other
checks defined in the `.pre-commit-config.yaml` file) on the files you've
staged.

* If a file is reformatted, the commit will fail, and you must `git add` the
re-formatted file(s) and commit again.
* If all checks pass, the commit will proceed normally.

Using `pre-commit` ensures that your changes will pass the formatting checks in
the CI build!
