include ../Config.make 

SRC= \
	coord2csr.c  \
	csc_nnz.c \
	csc_kron_mult.c  \
	csc_matmul_post.c  \
	csc_matmul_pre.c  \
	csr2den.c \
	csr_den_kron_mult.c \
	csr_nnz.c \
	csr_transpose.c \
	csr_kron_mult.c  \
	csr_kron_submatrix.c  \
	csr_matmul_post.c  \
	csr_matmul_pre.c  \
	csr_submatrix.c  \
	den2csr.c \
	den_csr_kron_mult.c \
	den_copymat.c \
	den_zeros.c \
	den_transpose.c \
	den_submatrix.c \
        den_gen_matrix.c \
	den_kron_submatrix.c \
	den_kron_form.c  \
	den_kron_mult.c  \
	den_matmul_post.c  \
	den_matmul_pre.c  \
	den_nnz.c  \
	estimate_kron_cost.c

OBJ= \
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
	ranlibkronutil libkronutil.a

test1: libkronutil.a test1.o
	$(CC) $(CFLAGS) -o test1 test1.o libkronutil.a $(LIBS)

test2: libkronutil.a test2.o
	$(CC) $(CFLAGS) -o test2 test2.o libkronutil.a $(LIBS)

clean:
	touch csr_kron_mult.o libkronutil.a test1 test2
	rm *.o libkronutil.a test1 test2
