#ifndef IONGSERIALIZER_H
#define IONGSERIALIZER_H
#include "H5Cpp.h"
#include "Vector.h"

namespace PsimagLite {

template<typename T>
struct ToH5 {
	static const H5::PredType type;
	static const H5T_class_t super;
};

class IoNgSerializer {

public:

	IoNgSerializer(H5::H5File* hdf5file, H5::Group* hdf5group)
	    : hdf5file_(hdf5file), hdf5group_(hdf5group)
	{}

	void createGroup(String group)
	{
		hdf5file_->createGroup(group);
	}

	void writeToTag(String name, SizeType what)
	{
		hsize_t dims[1];
		dims[0] = 1;
		H5::DataSpace *dataspace = new H5::DataSpace(1, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->createDataSet(name,
		                                                                ToH5<SizeType>::type,
		                                                                *dataspace,
		                                                                dsCreatPlist));
		void* ptr = static_cast<SizeType*>(&what);
		dataset->write(ptr, ToH5<SizeType>::type);
		delete dataset;
		delete dataspace;
	}

	void writeToTag(String name, String what)
	{
		hsize_t dims[1];
		dims[0] = what.length();
		H5::DataSpace *dataspace = new H5::DataSpace(1, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->createDataSet(name,
		                                                                ToH5<char>::type,
		                                                                *dataspace,
		                                                                dsCreatPlist));
		void* ptr = static_cast<void*>(&what[0]);
		dataset->write(ptr, ToH5<char>::type);
		delete dataset;
		delete dataspace;
	}

private:

	H5::H5File* hdf5file_;
	H5::Group* hdf5group_;
};
}
#endif // IONGSERIALIZER_H
