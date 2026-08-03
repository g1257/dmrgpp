
# Internal Packages and Includes

This memo specifies how internal dependencies are managed.
It does not address third party libraries (TPLs)

## General Rules
Do not use the special directory "..", and follow
https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#rs-incform

## Packages
DMRG++ is composed of the following packages: PsimagLite, LanczosPlusPlus,
dmrg, and cincuenta. TODO: Maybe names should all be capitalized.

![Packages](codemap.png)

## How includes should be done
### Package includes one of its files
If package (say LanczosPlusPlus) wants to include a file
from within itself it should include it directly or may prepend it
with a directory internal to the package. These examples are correct.

```cpp
#include "MyFile.hpp"
#include "Engine/MyOtherFile.hpp"
```

It must not use the special directory "..", and must not
prepend it with its own name.
These examples for package LanczosPlusPlus are all incorrect:

```cpp
#include "../MyFile.hpp"
#include "LanczosPlusPlus/MyFile.hpp"
#include "LanczosPlusPlus/Engine/MyOtherFile.hpp"
```

### Package includes a file from another Package
If a package wants to include a file from another package it must
prepend the include with package's name and inner location **in angle brackets**.
These examples for package dmrg are correct:

```cpp
#include <LanczosPlusPlus/MyFile.hpp>
#include <LanczosPlusPlus/Engine/MyOtherFile.hpp>
```

These examples for package dmrg are all incorrect:

```cpp
#include "MyFile.hpp"
#include "Engine/MyOtherFile.hpp"
#include "../LanczosPlusPlus/MyFile.hpp"
#include "LanczosPlusPlus/YetAnotherFile.hpp"
```

These includes are allowed:

```cpp
#include <PsimagLite/Matrix.h>
#include <PsimagLite/Io/IoNg.h>
```
This is bad
```cpp
#include <Matrix.h>
#include <IoNg.h>
```

