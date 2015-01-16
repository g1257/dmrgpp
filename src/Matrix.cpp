#include "Matrix.h"
#include "String.h"

namespace PsimagLite {

void checkBlasStatus(int info,PsimagLite::String msg)
{
	if (info == 0) return;

	PsimagLite::String str = msg;
	str += " failed with info = " + ttos(info) + "\n";
	throw RuntimeError(str);
}

void diag(Matrix<double> &m,Vector<double> ::Type& eigs,char option)
{
#ifdef NO_LAPACK
	throw RuntimeError("diag: dsyev_: NO LAPACK!\n");
#else
	char jobz=option;
	char uplo='U';
	int n=m.n_row();
	int lda=m.n_col();
	Vector<double>::Type work(3);
	int info,lwork= -1;

	if (lda<=0) throw RuntimeError("lda<=0\n");

	eigs.resize(n);

	// query:
	dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}
	lwork = int(work[0])+1;
	work.resize(lwork+1);
	// real work:
	dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}
#endif
}

void diag(Matrix<std::complex<double> > &m,Vector<double> ::Type&eigs,char option)
{
#ifdef NO_LAPACK
	throw RuntimeError("diag: zheev: NO LAPACK!\n");
#else
	char jobz=option;
	char uplo='U';
	int n=m.n_row();
	int lda=m.n_col();
	Vector<std::complex<double> >::Type work(3);
	Vector<double>::Type rwork(3*n);
	int info,lwork= -1;

	eigs.resize(n);

	// query:
	zheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(eigs[0]),&(work[0]),&lwork,&(rwork[0]),&info);
	lwork = int(real(work[0]))+1;
	work.resize(lwork+1);
	// real work:
	zheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(eigs[0]),&(work[0]),&lwork,&(rwork[0]),&info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		throw RuntimeError("diag: zheev: failed with info!=0.\n");
	}
#endif
}

void diag(Matrix<float> &m,Vector<float> ::Type& eigs,char option)
{
#ifdef NO_LAPACK
	throw RuntimeError("diag: dsyev_: NO LAPACK!\n");
#else
	char jobz=option;
	char uplo='U';
	int n=m.n_row();
	int lda=m.n_col();
	Vector<float>::Type work(3);
	int info,lwork= -1;

	if (lda<=0) throw RuntimeError("lda<=0\n");

	eigs.resize(n);

	// query:
	ssyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}
	lwork = int(work[0])+1;
	work.resize(lwork+1);
	// real work:
	ssyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}
#endif
}

void diag(Matrix<std::complex<float> > &m,Vector<float> ::Type& eigs,char option)
{
#ifdef NO_LAPACK
	throw RuntimeError("diag: cheev: NO LAPACK!\n");
#else
	char jobz=option;
	char uplo='U';
	int n=m.n_row();
	int lda=m.n_col();
	Vector<std::complex<float> >::Type work(3);
	Vector<float>::Type rwork(3*n);
	int info,lwork= -1;

	eigs.resize(n);

	// query:
	cheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(eigs[0]),&(work[0]),&lwork,&(rwork[0]),&info);
	lwork = int(real(work[0]))+1;
	work.resize(lwork+1);
	// real work:
	cheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(eigs[0]),&(work[0]),&lwork,&(rwork[0]),&info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		throw RuntimeError("diag: cheev: failed with info!=0.\n");
	}
#endif
}

void geev(char jobvl,
          char jobvr,
          Matrix<std::complex<double> >& a,
          Vector<std::complex<double> >::Type& w,
          Matrix<std::complex<double> >& vl,
          Matrix<std::complex<double> >& vr)
{
#ifdef NO_LAPACK
	throw RuntimeError("diag: geev: NO LAPACK!\n");
#else
	int n = a.n_row();
	int lda = a.n_col();
	int ldvl = vl.n_row();
	int ldvr = vr.n_row();
	int info = 0;
	Vector<std::complex<double> >::Type work(10,0);
	Vector<double>::Type rwork(2*n+1,0);
	int lwork = -1;
	zgeev_(&jobvl,
	      &jobvr,
	      &n,
	      &(a(0,0)),
	      &lda,
	      &(w[0]),
	        &(vl(0,0)),
	        &ldvl,
	        &(vr(0,0)),
	        &ldvr,
	        &(work[0]),
	        &lwork,
	        &(rwork[0]),
	        &info);

	checkBlasStatus(info,"zgeev_");

	lwork = std::real(work[0]);
	work.resize(lwork,0.0);

	zgeev_(&jobvl,
	      &jobvr,
	      &n,
	      &(a(0,0)),
	      &lda,
	      &(w[0]),
	        &(vl(0,0)),
	        &ldvl,
	        &(vr(0,0)),
	        &ldvr,
	        &(work[0]),
	        &lwork,
	        &(rwork[0]),
	        &info);

	checkBlasStatus(info,"zgeev_");
#endif
}

} // namespace PsimagLite

