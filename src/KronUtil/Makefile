include ../Config.make 

all: libkronutil.a test1 test2

CPPFLAGS += -I../../../PsimagLite/src
LDFLAGS += -L. -L../../../PsimagLite/lib

SRC=    kronUtilComplex.cpp \
	coord2csr.cpp  \
	csc_nnz.cpp \
	csc_kron_mult.cpp  \
	csc_matmul_post.cpp  \
	csc_matmul_pre.cpp  \
	csr2den.cpp \
	csr_den_kron_mult.cpp \
	csr_nnz.cpp \
	csr_transpose.cpp \
	csr_kron_mult.cpp  \
	csr_kron_submatrix.cpp  \
	csr_matmul_post.cpp  \
	csr_matmul_pre.cpp  \
	csr_submatrix.cpp  \
	den2csr.cpp \
	den_csr_kron_mult.cpp \
	den_copymat.cpp \
	den_zeros.cpp \
	den_transpose.cpp \
	den_submatrix.cpp \
        den_gen_matrix.cpp \
	den_kron_submatrix.cpp \
	den_kron_form.cpp  \
	den_kron_mult.cpp  \
	den_matmul_post.cpp  \
	den_matmul_pre.cpp  \
	den_nnz.cpp  \
	estimate_kron_cost.cpp

OBJ= \
        KronUtilComplex.o \
	coord2csr.o  \
	csc_nnz.o \
	csc_kron_mult.o  \
	csc_matmul_post.o  \
	csc_matmul_pre.o  \
	csr2den.o \
	csr_den_kron_mult.o \
	csr_nnz.o \
	csr_transpose.o \
	csr_kron_mult.o  \
	csr_kron_submatrix.o  \
	csr_matmul_post.o  \
	csr_matmul_pre.o  \
	csr_submatrix.o  \
	den2csr.o \
	den_csr_kron_mult.o \
	den_copymat.o \
	den_zeros.o \
	den_transpose.o \
	den_submatrix.o \
        den_gen_matrix.o \
	den_kron_submatrix.o \
	den_kron_form.o  \
	den_kron_mult.o  \
	den_matmul_post.o  \
	den_matmul_pre.o  \
	den_nnz.o  \
	estimate_kron_cost.o

libkronutil.a: $(OBJ)
	ar rc libkronutil.a $(OBJ)
	ranlib libkronutil.a

test1: libkronutil.a test1.o $(OBJ)
	$(CXX) $(CFLAGS) -o test1 test1.o libkronutil.a $(LDFLAGS)

test2: libkronutil.a test2.o
	$(CXX) $(CFLAGS) -o test2 test2.o libkronutil.a $(LDFLAGS)

clean:
	touch csr_kron_mult.o libkronutil.a test1 test2
	rm *.o libkronutil.a test1 test2
