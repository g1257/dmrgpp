#ifndef IONGSERIALIZER_H
#define IONGSERIALIZER_H
#include "H5Cpp.h"
#include "../Vector.h"
#include "../TypeToString.h"
#include "TypeToH5.h"
#include <cassert>
#include <stack>

namespace PsimagLite {

class IoNgSerializer {

public:

	enum WriteMode {NO_OVERWRITE, ALLOW_OVERWRITE};

	IoNgSerializer(H5::H5File* hdf5file)
	    : hdf5file_(hdf5file)
	{}

	void createGroup(String group)
	{
		hdf5file_->createGroup("Def/" + group);
	}

	/* write functions START */

	void write(String name2, SizeType what, WriteMode allowOverwrite = NO_OVERWRITE)
	{
		String name = "Def/" + name2;
		void* ptr = static_cast<SizeType*>(&what);

		if (allowOverwrite) {
			overwrite<SizeType>(name, ptr);
		} else {
			hsize_t dims[1];
			dims[0] = 1;
			internalWrite<SizeType>(name, ptr, dims, 1);
		}
	}

	void write(String name2, String what, WriteMode allowOverwrite = NO_OVERWRITE)
	{
		overwriteNotSupported(allowOverwrite);
		String name = "Def/" + name2;
		hsize_t dims[1];
		dims[0] = what.length();
		void* ptr = static_cast<void*>(&what[0]);
		internalWrite<char>(name, ptr, dims, 1);
	}

	template<typename T1, typename T2>
	void write(String name2,
	           const std::pair<T1, T2>& what,
	           WriteMode allowOverwrite = NO_OVERWRITE)
	{
		overwriteNotSupported(allowOverwrite);
		createGroup(name2);
		write(name2 + "/0", what.first);
		write(name2 + "/1", what.second);
	}

	void write(String name2,
	           const std::vector<bool>&,
	           WriteMode allowOverwrite = NO_OVERWRITE)
	{
		overwriteNotSupported(allowOverwrite);
		std::cerr<<"Vector of booleans with name "<<name2<<" cannot be printed ";
		std::cerr<<"FIXME TODO WARNING\n";
	}

	template<typename T>
	void write(String name2,
	           const std::vector<T>& what,
	           WriteMode allowOverwrite = NO_OVERWRITE,
	           typename EnableIf<Loki::TypeTraits<T>::isArith, int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		if (what.size() == 0) return;
		String name = "Def/" + name2;
		hsize_t dims[1];
		dims[0] = what.size();
		assert(0 < what.size());
		const void* ptr = static_cast<const void*>(&what[0]);
		internalWrite<T>(name, ptr, dims, 1);
	}

	template<typename T>
	void write(String name2,
	           const std::vector<std::complex<T> >& what,
	           WriteMode allowOverwrite = NO_OVERWRITE,
	           typename EnableIf<Loki::TypeTraits<T>::isArith, int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		String name = "Def/" + name2;
		hsize_t dims[1];
		dims[0] = 2*what.size();
		assert(0 < what.size());
		const void* ptr = static_cast<const void*>(&what[0]);
		internalWrite<T>(name, ptr, dims, 1);
	}

	template<typename T>
	void write(String name2,
	           const std::vector<std::vector<T> >& what,
	           WriteMode allowOverwrite = NO_OVERWRITE,
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
	           WriteMode allowOverwrite = NO_OVERWRITE,
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
	           WriteMode allowOverwrite = NO_OVERWRITE,
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
	           WriteMode allowOverwrite = NO_OVERWRITE,
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

	/* write functions END */

	// read functions START

	template<typename SomeType>
	void read(SomeType& value,
	          String name,
	          typename EnableIf<Loki::TypeTraits<SomeType>::isArith, int>::Type = 0)
	{
		void* ptr = static_cast<void *>(&value);
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet(name));
		dataset->read(ptr, TypeToH5<SomeType>::type);
		delete dataset;
	}

	void read(String& value,
	          String name)
	{
		void* ptr = static_cast<void *>(&(value[0])); // FIXME CHECK SIZE
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet("Def/" + name));
		dataset->read(ptr, TypeToH5<char>::type);
		delete dataset;
	}

	void read(std::vector<bool>&,
	          String name)
	{
		std::cerr<<"Vector of booleans with name "<<name<<" cannot be read ";
		std::cerr<<"FIXME TODO WARNING\n";
	}

	template<typename SomeType>
	void read(SomeType& value,
	          String name,
	          typename EnableIf<IsEnum<SomeType>::True, int>::Type = 0)
	{
		throw RuntimeError("Cannot read label " + name + " (enums not supported yet)\n");
	}

	template<typename T>
	void read(std::vector<T>& what,
	          String name,
	          typename EnableIf<Loki::TypeTraits<T>::isArith, int>::Type = 0)
	{
		void* ptr = static_cast<void *>(&(what[0]));  // FIXME CHECK SIZE
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet(name));
		dataset->read(ptr, TypeToH5<T>::type);
		delete dataset;
	}

	template<typename T>
	void read(std::vector<std::complex<T> >& what,
	          String name,
	          typename EnableIf<Loki::TypeTraits<T>::isArith, int>::Type = 0)
	{
		void* ptr = static_cast<void *>(&(what[0]));  // FIXME CHECK SIZE
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet(name));
		dataset->read(ptr, TypeToH5<T>::type);
		delete dataset;
	}

	template<typename T>
	void read(std::vector<T>& what,
	          String name,
	          typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith,
	          int>::Type = 0)
	{
		throw RuntimeError("Cannot read " + name + " (vector<compound> not yet supported)\n");
	}

	template<typename T>
	void read(std::stack<T>& what,
	          String name)
	{
		throw RuntimeError("Cannot read label " + name + " (stacks not supported yet)\n");
	}

	// read functions END

private:

	void overwriteNotSupported(WriteMode mode)
	{
		if (mode == NO_OVERWRITE) return;
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
	void internalWrite(String name, const void* ptr, hsize_t dims[], SizeType ndims)
	{
		H5::DataSpace *dataspace = new H5::DataSpace(ndims, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->createDataSet(name,
		                                                                TypeToH5<SomeType>::type,
		                                                                *dataspace,
		                                                                dsCreatPlist));
		dataset->write(ptr, TypeToH5<SomeType>::type);
		delete dataset;
		delete dataspace;
	}

	H5::H5File* hdf5file_;
};
}
#endif // IONGSERIALIZER_H
