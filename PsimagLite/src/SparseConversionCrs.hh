#ifndef SPARSECONVERSIONCRS_HH
#define SPARSECONVERSIONCRS_HH
#include "CrsMatrix.h"
#include "Matrix.h"

namespace PsimagLite {

template <typename T> class SparseConversions {
public:

	static CrsMatrix<value_type> toCRS(const T& a) { return a.toCRS(); }
};

template <typename T> class SparseConversions<CrsMatrix<T>> {
public:

	static CrsMatrix<T> toCRS(const CrsMatrix<T>& crs) { return crs; }
};
}
#endif // SPARSECONVERSIONCRS_HH
