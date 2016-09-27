#ifndef USE_C_VECTOR_SIZE_H
#define USE_C_VECTOR_SIZE_H

#ifdef USE_COMPRESSED_VECTORS
typedef PsimagLite::Vector<SizeType>::Type CvectorSizeType;
#else
typedef PsimagLite::Vector<SizeType>::Type CvectorSizeType;
#endif

#endif // USE_C_VECTOR_SIZE_H

