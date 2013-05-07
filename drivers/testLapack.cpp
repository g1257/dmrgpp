#include <iostream>
#include <vector>
#include <cstdlib>

extern "C" void dsyev_(char *,char *,int *,double *,int *, double *,double *,int *,int *);

int main()
{
	int n = 100;
	typename Vector<double>::Type m(n*n);
	// fill "matrix"
	for (size_t i=0;i<size_t(n*n);i++) m[i] =  5*drand48();
	// symmetrize:
	for (size_t i=0;i<size_t(n);i++)
		for (size_t j=i+1;j<size_t(n);j++)
			m[i+j*n] = m[j+i*n];

	typename Vector<double>::Type eigs(n);
	char jobz='V';
	char uplo='U';
	int lda=n;
	typename Vector<double>::Type work(3);
	int info = 0;
	int lwork= -1;

	// query:
	dsyev_(&jobz,&uplo,&n,&(m[0]),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
	if (info!=0) {
		std::cerr<<"diag: dsyev_: failed with info="<<info<<"\n";
		return 1;
	}
	lwork = int(work[0])+1;
	work.resize(lwork+1);
	
	// real work:
	dsyev_(&jobz,&uplo,&n,&(m[0]),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
	if (info!=0) {
		std::cerr<<"diag: dsyev_: failed with info="<<info<<"\n";
		return 1;
	}
}
