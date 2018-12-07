#include "Matrix.h"

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
	int n=m.rows();
	int lda=m.cols();
	Vector<double>::Type work(3);
	int info;
	int lwork= -1;

	if (lda<=0) throw RuntimeError("lda<=0\n");

	eigs.resize(n);

	// query:
	dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}

	const int NB = 256;
	lwork = std::max(1 + static_cast<int>(work[0]), (NB + 2)*n);
	work.resize(lwork);
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
	int n=m.rows();
	int lda=m.cols();
	Vector<std::complex<double> >::Type work(3);
	Vector<double>::Type rwork(3*n);
	int info;
	int lwork= -1;

	eigs.resize(n);

	// query:
	zheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(eigs[0]),&(work[0]),&lwork,&(rwork[0]),&info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		throw RuntimeError("diag: zheev_: failed with info!=0.\n");
	}

	const int NB = 256;
	lwork = std::max(1 + static_cast<int>(std::real(work[0])), (NB + 2)*n);
	work.resize(lwork);
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
	int n=m.rows();
	int lda=m.cols();
	Vector<float>::Type work(3);
	int info;
	int lwork= -1;

	if (lda<=0) throw RuntimeError("lda<=0\n");

	eigs.resize(n);

	// query:
	ssyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		throw RuntimeError("diag: dsyev_: failed with info!=0.\n");
	}

	const int NB = 256;
	lwork = std::max(1 + static_cast<int>(work[0]), (NB + 2)*n);
	work.resize(lwork);

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
	int n=m.rows();
	int lda=m.cols();
	Vector<std::complex<float> >::Type work(3);
	Vector<float>::Type rwork(3*n);
	int info,lwork= -1;

	eigs.resize(n);

	// query:
	cheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(eigs[0]),&(work[0]),&lwork,&(rwork[0]),&info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		throw RuntimeError("diag: cheev_: failed with info!=0.\n");
	}

	const int NB = 256;
	lwork = std::max(1 + static_cast<int>(std::real(work[0])), (NB + 2)*n);
	work.resize(lwork);

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
	int n = a.rows();
	int lda = a.cols();
	int ldvl = vl.rows();
	int ldvr = vr.rows();
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

	const int NB = 256;
	lwork = std::max(1 + static_cast<int>(std::real(work[0])), (NB + 2)*n);
	work.resize(lwork);

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

void inverse(Matrix<std::complex<double> > &m)
{
#ifdef NO_LAPACK
	throw RuntimeError("inverse: NO LAPACK!\n");
#else
	int n = m.rows();
	int info = 0;
	Vector<int>::Type ipiv(n,0);
	psimag::LAPACK::zgetrf_(&n,&n,&(m(0,0)),&n,&(ipiv[0]),&info);
	int lwork = -1;
	Vector<std::complex<double> >::Type work(2);
	psimag::LAPACK::zgetri_(&n,&(m(0,0)),&n,&(ipiv[0]),&(work[0]),&lwork,&info);
	lwork = static_cast<int>(PsimagLite::real(work[0]));
	work.resize(lwork+2);
	psimag::LAPACK::zgetri_(&n,&(m(0,0)),&n,&(ipiv[0]),&(work[0]),&lwork,&info);
	String s = "zgetri_ failed\n";
	if (info!=0) throw RuntimeError(s.c_str());
#endif
}

void inverse(Matrix<double> &m)
{
#ifdef NO_LAPACK
	throw RuntimeError("inverse: NO LAPACK!\n");
#else
	int n = m.rows();
	int info = 0;
	Vector<int>::Type ipiv(n,0);
	psimag::LAPACK::dgetrf_(&n,&n,&(m(0,0)),&n,&(ipiv[0]),&info);
	int lwork = -1;
	Vector<double>::Type work(2);
	psimag::LAPACK::dgetri_(&n,&(m(0,0)),&n,&(ipiv[0]),&(work[0]),&lwork,&info);
	lwork = static_cast<int>(work[0]);
	work.resize(lwork+2);
	psimag::LAPACK::dgetri_(&n,&(m(0,0)),&n,&(ipiv[0]),&(work[0]),&lwork,&info);
	String s = "dgetri_ failed\n";
	if (info!=0) throw RuntimeError(s.c_str());
#endif
}

void svd(char jobz,
         Matrix<double>& a,
         Vector<double>::Type& s,
         Matrix<double>& vt,
         SvdRealFuncPtr<double>::F backend,
         int counter)
{
#ifdef NO_LAPACK
	throw RuntimeError("svd: dgesdd_: NO LAPACK!\n");
#else
	if (jobz != 'A' && jobz != 'S') {
		String msg("svd: jobz must be either A or S");
		String jobzString = " ";
		jobzString[0] = jobz;
		throw RuntimeError(msg + ", not " + jobzString + "\n");
	}

	int m = a.rows();
	int n = a.cols();
	int lda = m;
	int min = (m<n) ? m : n;

	s.resize(min);
	int ldu = m;
	int ucol = (jobz == 'A') ? m : min;
	Matrix<double> u(ldu,ucol);
	int ldvt = (jobz == 'A') ? n : min;
	vt.resize(ldvt,n);

	Vector<double>::Type work(100,0);
	int info = 0;
	Vector<int>::Type iwork(8*min,0);

	// query optimal work
	int lwork = -1;
	backend(&jobz,
	        &m,
	        &n,
	        &(a(0,0)),
	        &lda,
	        &(s[0]),
	        &(u(0,0)),
	        &ldu,
	        &(vt(0,0)),
	        &ldvt,
	        &(work[0]),
	        &lwork,
	        &(iwork[0]),
	        &info);
	if (info!=0) {
		String str(__FILE__);
		str += " svd(...) failed at workspace size calculation ";
		str += "with info=" + ttos(info) + "\n";
		throw RuntimeError(str.c_str());
	}

	lwork = int(work[0]) + (m+n)*256;
	work.resize(lwork+10);
	// real work:
	backend(&jobz,
	        &m,
	        &n,
	        &(a(0,0)),
	        &lda,
	        &(s[0]),
	        &(u(0,0)),
	        &ldu,
	        &(vt(0,0)),
	        &ldvt,
	        &(work[0]),
	        &lwork,
	        &(iwork[0]),
	        &info);
	if (info != 0) {
		String str(__FILE__);
		str += " " + ttos(__LINE__);
		str += " svd(...) failed with info=" + ttos(info);
		str += " matrix is " + ttos(a.rows()) + " " + ttos(a.cols()) + "\n";
		if (info < 0 || counter > 0)
			throw RuntimeError(str);

		std::cerr<<str;
		std::cerr<<"Will try with fallback...\n";
		svd(jobz, a, s, vt, psimag::LAPACK::dgesdd_, 1);
	}

	a = u;
#endif
}

void svd(char jobz,
         Matrix<float>& a,
         Vector<float>::Type& s,
         Matrix<float>& vt,
         SvdRealFuncPtr<float>::F backend,
         int counter)
{
#ifdef NO_LAPACK
	throw RuntimeError("svd: dgesdd_: NO LAPACK!\n");
#else
	if (jobz != 'A' && jobz != 'S') {
		String msg("svd: jobz must be either A or S");
		String jobzString = " ";
		jobzString[0] = jobz;
		throw RuntimeError(msg + ", not " + jobzString + "\n");
	}

	int m = a.rows();
	int n = a.cols();
	int lda = m;
	int min = (m<n) ? m : n;

	s.resize(min);
	int ldu = m;
	int ucol = (jobz == 'A') ? m : min;
	Matrix<float> u(ldu,ucol);
	int ldvt = (jobz == 'A') ? n : min;
	vt.resize(ldvt,n);

	Vector<float>::Type work(100,0);
	int info = 0;
	Vector<int>::Type iwork(8*min,0);

	// query optimal work
	int lwork = -1;
	backend(&jobz,
	        &m,
	        &n,
	        &(a(0,0)),
	        &lda,
	        &(s[0]),
	        &(u(0,0)),
	        &ldu,
	        &(vt(0,0)),
	        &ldvt,
	        &(work[0]),
	        &lwork,
	        &(iwork[0]),
	        &info);
	if (info!=0) {
		String str(__FILE__);
		str += " " + ttos(__LINE__);
		str += " svd(...) failed at workspace size calculation ";
		str += "with info=" + ttos(info) + "\n";
		throw RuntimeError(str.c_str());
	}

	lwork = int(work[0]) + (m+n)*256;
	work.resize(lwork+10);
	// real work:
	backend(&jobz,
	        &m,
	        &n,
	        &(a(0,0)),
	        &lda,
	        &(s[0]),
	        &(u(0,0)),
	        &ldu,
	        &(vt(0,0)),
	        &ldvt,
	        &(work[0]),
	        &lwork,
	        &(iwork[0]),
	        &info);
	if (info != 0) {
		String str(__FILE__);
		str += " " + ttos(__LINE__);
		str += " svd(...) failed with info=" + ttos(info);
		str += " matrix is " + ttos(a.rows()) + " " + ttos(a.cols()) + "\n";
		if (info < 0 || counter > 0)
			throw RuntimeError(str);

		std::cerr<<str;
		std::cerr<<"Will try with fallback...\n";
		svd(jobz, a, s, vt, psimag::LAPACK::sgesdd_, 1);
	}

	a = u;
#endif
}

void svd(char jobz,
         Matrix<std::complex<double> >& a,
         Vector<double>::Type& s,
         Matrix<std::complex<double> >& vt,
         SvdCmplxFuncPtr<double>::F backend,
         int counter)
{
#ifdef NO_LAPACK
	throw RuntimeError("svd: dgesdd_: NO LAPACK!\n");
#else
	if (jobz != 'A' && jobz != 'S') {
		String msg("svd: jobz must be either A or S");
		String jobzString = " ";
		jobzString[0] = jobz;
		throw RuntimeError(msg + ", not " + jobzString + "\n");
	}

	int m = a.rows();
	int n = a.cols();
	int lda = m;
	int min = (m<n) ? m : n;

	s.resize(min);
	int ldu = m;
	int ucol = (jobz == 'A') ? m : min;
	Matrix<std::complex<double> > u(ldu,ucol);
	int ldvt = (jobz == 'A') ? n : min;
	vt.resize(ldvt,n);
	int lrwork = 2.0*min*std::max(5*min+7,2*std::max(m,n)+2*min+1);
	Vector<double>::Type rwork(lrwork, 0.0);

	Vector<std::complex<double> >::Type work(100,0);
	int info = 0;
	Vector<int>::Type iwork(8*min,0);

	// query optimal work
	int lwork = -1;
	backend(&jobz,
	        &m,
	        &n,
	        &(a(0,0)),
	        &lda,
	        &(s[0]),
	        &(u(0,0)),
	        &ldu,
	        &(vt(0,0)),
	        &ldvt,
	        &(work[0]),
	        &lwork,
	        &(rwork[0]),
	        &(iwork[0]),
	        &info);
	if (info != 0) {
		String str(__FILE__);
		str += " " + ttos(__LINE__);
		str += " svd(...) failed at workspace size calculation ";
		str += "with info=" + ttos(info) + "\n";
		throw RuntimeError(str.c_str());
	}

	lwork = int(std::real(work[0])) + (m+n)*256;
	work.resize(lwork+10);
	// real work:
	backend(&jobz,
	        &m,
	        &n,
	        &(a(0,0)),
	        &lda,
	        &(s[0]),
	        &(u(0,0)),
	        &ldu,
	        &(vt(0,0)),
	        &ldvt,
	        &(work[0]),
	        &lwork,
	        &(rwork[0]),
	        &(iwork[0]),
	        &info);
	if (info != 0) {
		String str(__FILE__);
		str += " " + ttos(__LINE__);
		str += " svd(...) failed with info=" + ttos(info);
		str += " matrix is " + ttos(a.rows()) + " " + ttos(a.cols()) + "\n";
		if (info < 0 || counter > 0)
			throw RuntimeError(str);

		std::cerr<<str;
		std::cerr<<"Will try with fallback...\n";
		svd(jobz, a, s, vt, psimag::LAPACK::zgesdd_, 1);
	}

	a = u;
#endif
}

void svd(char jobz,
         Matrix<std::complex<float> >& a,
         Vector<float>::Type& s,
         Matrix<std::complex<float> >& vt,
         SvdCmplxFuncPtr<float>::F backend,
         int counter)
{
#ifdef NO_LAPACK
	throw RuntimeError("svd: dgesdd_: NO LAPACK!\n");
#else
	if (jobz != 'A' && jobz != 'S') {
		String msg("svd: jobz must be either A or S");
		String jobzString = " ";
		jobzString[0] = jobz;
		throw RuntimeError(msg + ", not " + jobzString + "\n");
	}

	int m = a.rows();
	int n = a.cols();
	int lda = m;
	int min = (m<n) ? m : n;

	s.resize(min);
	int ldu = m;
	int ucol = (jobz == 'A') ? m : min;
	Matrix<std::complex<float> > u(ldu,ucol);
	int ldvt = (jobz == 'A') ? n : min;
	vt.resize(ldvt,n);
	int lrwork = 2.0*min*std::max(5*min+7,2*std::max(m,n)+2*min+1);
	Vector<float>::Type rwork(lrwork, 0.0);

	Vector<std::complex<float> >::Type work(100,0);
	int info = 0;
	Vector<int>::Type iwork(8*min,0);

	// query optimal work
	int lwork = -1;
	backend(&jobz,
	        &m,
	        &n,
	        &(a(0,0)),
	        &lda,
	        &(s[0]),
	        &(u(0,0)),
	        &ldu,
	        &(vt(0,0)),
	        &ldvt,
	        &(work[0]),
	        &lwork,
	        &(rwork[0]),
	        &(iwork[0]),
	        &info);
	if (info!=0) {
		String str(__FILE__);
		str += " " + ttos(__LINE__);
		str += " svd(...) failed at workspace size calculation ";
		str += "with info=" + ttos(info) + "\n";
		throw RuntimeError(str.c_str());
	}

	lwork = int(std::real(work[0])) + (m+n)*256;
	work.resize(lwork+10);
	// real work:
	backend(&jobz,
	        &m,
	        &n,
	        &(a(0,0)),
	        &lda,
	        &(s[0]),
	        &(u(0,0)),
	        &ldu,
	        &(vt(0,0)),
	        &ldvt,
	        &(work[0]),
	        &lwork,
	        &(rwork[0]),
	        &(iwork[0]),
	        &info);
	if (info != 0) {
		String str(__FILE__);
		str += " " + ttos(__LINE__);
		str += " svd(...) failed with info=" + ttos(info);
		str += " matrix is " + ttos(a.rows()) + " " + ttos(a.cols()) + "\n";
		if (info < 0 || counter > 0)
			throw RuntimeError(str);

		std::cerr<<str;
		std::cerr<<"Will try with fallback...\n";
		svd(jobz, a, s, vt, psimag::LAPACK::cgesdd_, 1);
	}

	a = u;
#endif
}

} // namespace PsimagLite

