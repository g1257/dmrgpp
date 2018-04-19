#ifndef IONGSERIALIZER_H
#define IONGSERIALIZER_H
#include "H5Cpp.h"
#include "../Vector.h"
#include "../TypeToString.h"
#include "TypeToH5.h"

namespace PsimagLite {

class IoNgSerializer {

public:

	enum WriteMode {NO, ALLOW_OVERWRITE};

	IoNgSerializer(H5::H5File* hdf5file, H5::Group* hdf5group)
	    : hdf5file_(hdf5file), hdf5group_(hdf5group)
	{}

	void createGroup(String group)
	{
		hdf5file_->createGroup("Def/" + group);
	}

	void write(String name2, SizeType what, WriteMode allowOverwrite = NO)
	{
		String name = "Def/" + name2;
		void* ptr = static_cast<SizeType*>(&what);

		if (allowOverwrite) {
			overwrite<SizeType>(name, ptr);
		} else {
			hsize_t dims[1];
			dims[0] = 1;
			writeNew<SizeType>(name, ptr, dims, 1);
		}
	}

	void write(String name2, String what, WriteMode allowOverwrite = NO)
	{
		overwriteNotSupported(allowOverwrite);
		String name = "Def/" + name2;
		hsize_t dims[1];
		dims[0] = what.length();
		H5::DataSpace *dataspace = new H5::DataSpace(1, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->createDataSet(name,
		                                                                TypeToH5<char>::type,
		                                                                *dataspace,
		                                                                dsCreatPlist));
		void* ptr = static_cast<void*>(&what[0]);
		dataset->write(ptr, TypeToH5<char>::type);
		delete dataset;
		delete dataspace;
	}

	template<typename T1, typename T2>
	void write(String name2,
	           const std::pair<T1, T2>& what,
	           WriteMode allowOverwrite = NO)
	{
		overwriteNotSupported(allowOverwrite);
		createGroup(name2);
		write(name2 + "/0", what.first);
		write(name2 + "/1", what.second);
	}

	void write(String name2,
	           const std::vector<bool>&,
	           WriteMode allowOverwrite = NO)
	{
		overwriteNotSupported(allowOverwrite);
		std::cerr<<"Vector of booleans with name "<<name2<<" cannot be printed ";
		std::cerr<<"FIXME TODO WARNING\n";
	}

	template<typename T>
	void write(String name2,
	           const std::vector<T>& what,
	           WriteMode allowOverwrite = NO,
	           typename EnableIf<Loki::TypeTraits<T>::isArith, int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		String name = "Def/" + name2;
		hsize_t dims[1];
		dims[0] = what.size();
		H5::DataSpace *dataspace = new H5::DataSpace(1, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		internalWrite<T>(&(what[0]), name, *dataspace, dsCreatPlist);
		delete dataspace;
	}

	template<typename T>
	void write(String name2,
	           const std::vector<std::complex<T> >& what,
	           WriteMode allowOverwrite = NO,
	           typename EnableIf<Loki::TypeTraits<T>::isArith, int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		String name = "Def/" + name2;
		hsize_t dims[1];
		dims[0] = 2*what.size();
		H5::DataSpace *dataspace = new H5::DataSpace(1, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		internalWrite<T>(&(what[0]), name, *dataspace, dsCreatPlist);
		delete dataspace;
	}

	template<typename T>
	void write(String name2,
	           const std::vector<std::vector<T> >& what,
	           WriteMode allowOverwrite = NO,
	           typename EnableIf<Loki::TypeTraits<typename Real<T>::Type>::isArith,
	           int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		SizeType n = what.size();
		createGroup(name2);
		write(name2 + "/Size", n);
		for (SizeType i = 0; i < n; ++i)
			write(name2 + "/" + typeToString(i), what[i]);
	}

	template<typename T1, typename T2>
	void write(String name2,
	           const std::vector<std::pair<T1, T2> >& what,
	           WriteMode allowOverwrite = NO,
	           typename EnableIf<
	           Loki::TypeTraits<T1>::isArith && Loki::TypeTraits<T2>::isArith,
	           int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		SizeType n = what.size();
		createGroup(name2);
		write(name2 + "/Size", n);
		for (SizeType i = 0; i < n; ++i)
			write(name2 + "/" + typeToString(i), what[i]);
	}

	template<typename T>
	void write(String name2,
	           const std::vector<T>& what,
	           WriteMode allowOverwrite = NO,
	           typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith,
	           int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		SizeType n = what.size();
		createGroup(name2);
		write(name2 + "/Size", n);
		for (SizeType i = 0; i < n; ++i)
			what[i].write(name2 + "/" + typeToString(i), *this);
	}

	template<typename T>
	void write(String name2,
	           const std::vector<T*>& what,
	           WriteMode allowOverwrite = NO,
	           typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith,
	           int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		SizeType n = what.size();
		createGroup(name2);
		write(name2 + "/Size", n);
		for (SizeType i = 0; i < n; ++i)
			what[i]->write(name2 + "/" + typeToString(i), *this);
	}

private:

	void overwriteNotSupported(WriteMode mode) const
	{
		if (mode == NO) return;
		throw RuntimeError("Overwrite not supported for this type\n");
	}

	template<typename SomeType>
	void overwrite(String name, void* ptr)
	{
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet(name));
		dataset->write(ptr, TypeToH5<SomeType>::type);
		delete dataset;
	}

	template<typename SomeType>
	void writeNew(String name, void* ptr, hsize_t dims[], SizeType ndims)
	{
		H5::DataSpace *dataspace = new H5::DataSpace(ndims, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		internalWrite<SomeType>(ptr, name, *dataspace, dsCreatPlist);
		delete dataspace;
	}

	template<typename T>
	void internalWrite(const void *ptr,
	                   String name,
	                   H5::DataSpace& dataspace,
	                   H5::DSetCreatPropList& dsCreatPlist)
	{

		H5::DataSet* dataset = new H5::DataSet(hdf5file_->createDataSet(name,
		                                                                TypeToH5<T>::type,
		                                                                dataspace,
		                                                                dsCreatPlist));
		dataset->write(ptr, TypeToH5<T>::type);
		delete dataset;
	}

	H5::H5File* hdf5file_;
	H5::Group* hdf5group_;
};
}
#endif // IONGSERIALIZER_H
