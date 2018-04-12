#ifndef IONGSERIALIZER_H
#define IONGSERIALIZER_H
#include "H5Cpp.h"
#include "Vector.h"
#include "TypeToString.h"

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
		hdf5file_->createGroup("Def/" + group);
	}

	void overwrite(String name2, SizeType what)
	{
		String name = "Def/" + name2;
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet(name));
		void* ptr = static_cast<SizeType*>(&what);
		dataset->write(ptr, ToH5<SizeType>::type);
		delete dataset;
	}

	void writeToTag(String name2, SizeType what)
	{
		String name = "Def/" + name2;
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

	void writeToTag(String name2, String what)
	{
		String name = "Def/" + name2;
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

	template<typename T1, typename T2>
	void writeToTag(String name2,
	                const std::pair<T1, T2>& what)
	{
		createGroup(name2);
		writeToTag(name2 + "/0", what.first);
		writeToTag(name2 + "/1", what.second);
	}

	void writeToTag(String name2,
	                const std::vector<bool>&)
	{
		std::cerr<<"Vector of booleans with name "<<name2<<" cannot be printed ";
		std::cerr<<"FIXME TODO WARNING\n";
	}

	template<typename T>
	void writeToTag(String name2,
	                const std::vector<T>& what,
	                typename EnableIf<Loki::TypeTraits<T>::isArith, int>::Type = 0)
	{
		String name = "Def/" + name2;
		hsize_t dims[1];
		dims[0] = what.size();
		H5::DataSpace *dataspace = new H5::DataSpace(1, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		internalWrite<T>(&(what[0]), name, *dataspace, dsCreatPlist);
		delete dataspace;
	}

	template<typename T>
	void writeToTag(String name2,
	                const std::vector<std::complex<T> >& what,
	                typename EnableIf<Loki::TypeTraits<T>::isArith, int>::Type = 0)
	{
		String name = "Def/" + name2;
		hsize_t dims[1];
		dims[0] = 2*what.size();
		H5::DataSpace *dataspace = new H5::DataSpace(1, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		internalWrite<T>(&(what[0]), name, *dataspace, dsCreatPlist);
		delete dataspace;
	}

	template<typename T>
	void writeToTag(String name2,
	                const std::vector<std::vector<T> >& what,
	                typename EnableIf<Loki::TypeTraits<typename Real<T>::Type>::isArith,
	                int>::Type = 0)
	{
		SizeType n = what.size();
		createGroup(name2);
		writeToTag(name2 + "/SIZE", n);
		for (SizeType i = 0; i < n; ++i)
			writeToTag(name2 + "/" + typeToString(i), what[i]);
	}

	template<typename T1, typename T2>
	void writeToTag(String name2,
	                const std::vector<std::pair<T1, T2> >& what,
	                typename EnableIf<
	                Loki::TypeTraits<T1>::isArith && Loki::TypeTraits<T2>::isArith,
	                int>::Type = 0)
	{
		SizeType n = what.size();
		createGroup(name2);
		writeToTag(name2 + "/SIZE", n);
		for (SizeType i = 0; i < n; ++i)
			writeToTag(name2 + "/" + typeToString(i), what[i]);
	}

	template<typename T>
	void writeToTag(String name2,
	                const std::vector<T>& what,
	                typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith,
	                int>::Type = 0)
	{
		SizeType n = what.size();
		createGroup(name2);
		writeToTag(name2 + "/SIZE", n);
		for (SizeType i = 0; i < n; ++i)
			what[i].serialize(name2 + "/" + typeToString(i), *this);
	}

	template<typename T>
	void writeToTag(String name2,
	                const std::vector<T*>& what,
	                typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith,
	                int>::Type = 0)
	{
		SizeType n = what.size();
		createGroup(name2);
		writeToTag(name2 + "/SIZE", n);
		for (SizeType i = 0; i < n; ++i)
			what[i]->serialize(name2 + "/" + typeToString(i), *this);
	}

private:

	template<typename T>
	void internalWrite(const void *ptr,
	                   String name,
	                   H5::DataSpace& dataspace,
	                   H5::DSetCreatPropList& dsCreatPlist)
	{

		H5::DataSet* dataset = new H5::DataSet(hdf5file_->createDataSet(name,
		                                                                ToH5<T>::type,
		                                                                dataspace,
		                                                                dsCreatPlist));
		dataset->write(ptr, ToH5<T>::type);
		delete dataset;
	}

	H5::H5File* hdf5file_;
	H5::Group* hdf5group_;
};
}
#endif // IONGSERIALIZER_H
