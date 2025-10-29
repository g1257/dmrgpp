AC_INIT(PsiDrivers, 1.0)

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AC_PROG_CXX

AX_LAPACK(, echo "ERROR: LAPACK not Found"; exit)

dnl Process Makefile.in to create Makefile
AC_OUTPUT(Makefile)

