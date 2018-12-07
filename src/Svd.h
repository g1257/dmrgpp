#ifndef SVD_H
#define SVD_H
#include "Matrix.h"

namespace PsimagLite {

template<typename ComplexOrRealType>
class Svd {

public:

	typedef typename Real<ComplexOrRealType>::Type RealType;

	Svd(String name = "gesdd") : name_(name)
	{}

	bool canTryAgain() const
	{
		return (name_ == "gesdd");
	}

	String name() const { return name_; }

	void operator()(char jobz,
	                Matrix<ComplexOrRealType>& a,
	                typename Vector<RealType>::Type& s,
	                Matrix<ComplexOrRealType>& vt)
	{
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
		Matrix<ComplexOrRealType> u(ldu,ucol);
		int ldvt = (jobz == 'A') ? n : min;
		vt.resize(ldvt,n);
		int lrwork = 2.0*min*std::max(5*min+7,2*std::max(m,n)+2*min+1);
		typename Vector<typename Real<ComplexOrRealType>::Type>::Type rwork(lrwork, 0.0);

		typename Vector<ComplexOrRealType>::Type work(100,0);
		int info = 0;
		Vector<int>::Type iwork(8*min,0);

		// query optimal work
		int lwork = -1;
		mycall(&jobz,
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
			str += " svd(...) failed at workspace size calculation ";
			str += "with info=" + ttos(info) + "\n";
			throw RuntimeError(str.c_str());
		}

		RealType lworkReal = PsimagLite::real(work[0]);
		lwork = static_cast<int>(lworkReal) + (m+n)*256;
		work.resize(lwork+10);

		// real work:
		mycall(&jobz,
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
			if (info < 0 || !canTryAgain())
				throw RuntimeError(str);

			std::cerr<<str;
			std::cerr<<"Will try with fallback...\n";
			name_ = "gesvd";
			operator()(jobz, a, s, vt);
		}

		a = u;
	}

private:

	void mycall(char* jobz,
	            int* m,
	            int* n,
	            ComplexOrRealType* a, // T*,
	            int* lda,
	            RealType* s,
	            ComplexOrRealType* u, //T*,
	            int* ldu,
	            ComplexOrRealType* vt, // T*,
	            int* ldvt,
	            ComplexOrRealType* work, // T*,
	            int* lwork,
	            RealType* rwork, // nothing
	            int* iwork,
	            int *info)
	{
		if (name_ == "gesdd") {
			psimag::LAPACK::GESDD(jobz,
			                      m,
			                      n,
			                      a,
			                      lda,
			                      s,
			                      u,
			                      ldu,
			                      vt,
			                      ldvt,
			                      work,
			                      lwork,
			                      rwork,
			                      iwork,
			                      info);
		} else if (name_ == "gesvd") {
			psimag::LAPACK::GESVD(jobz,
			                      jobz,
			                      m,
			                      n,
			                      a,
			                      lda,
			                      s,
			                      u,
			                      ldu,
			                      vt,
			                      ldvt,
			                      work,
			                      lwork,
			                      rwork,
			                      info);
		} else {
			throw PsimagLite::RuntimeError("Unknown backend " + name_ + "\n");
		}
	}

	Svd(const Svd&);

	Svd& operator=(const Svd&);

	String name_;
};

}
#endif // SVD_H
