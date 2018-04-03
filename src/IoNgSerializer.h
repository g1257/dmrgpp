#ifndef IONGSERIALIZER_H
#define IONGSERIALIZER_H
#include "H5Cpp.h"
#include "Vector.h"

namespace PsimagLite {

class IoNgSerializer {

public:

	IoNgSerializer(H5::H5File* hdf5file, H5::Group* hdf5group)
	    : hdf5file_(hdf5file), hdf5group_(hdf5group)
	{}

private:

	H5::H5File* hdf5file_;
	H5::Group* hdf5group_;
};
}
#endif // IONGSERIALIZER_H
