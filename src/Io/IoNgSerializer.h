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

	template<typename T>
	void write(String name2,
	           const T& what,
	           WriteMode allowOverwrite = NO_OVERWRITE,
	           typename EnableIf<Loki::TypeTraits<T>::isArith || IsEnum<T>::True, int>::Type = 0)
	{
		String name = "Def/" + name2;
		const void* ptr = static_cast<const T*>(&what);

		if (allowOverwrite) {
			overwrite<T>(name, ptr);
		} else {
			hsize_t dims[1];
			dims[0] = 1;
			internalWrite<T>(name, ptr, dims, 1);
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

	void write(String name2, bool b, WriteMode allowOverwrite = NO_OVERWRITE)
	{
		overwriteNotSupported(allowOverwrite);
		String name = "Def/" + name2;
		hsize_t dims[1];
		dims[0] = 1;
		unsigned char tmp[1];
		tmp[0] = (b) ? '1' : '0';
		const void* ptr = static_cast<const void*>(tmp);
		internalWrite<unsigned char>(name, ptr, dims, 1);
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

	// Note: THIS WILL EMPTY THE STACK OBJECT!
	template<typename T>
	void write(String name,
	           std::stack<T>& what,
	           WriteMode allowOverwrite = NO_OVERWRITE,
	           typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith,
	           int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		createGroup(name);
		write(name + "/Size", what.size());
		SizeType i = 0;
		while (what.size() > 0) {
			const T& t = what.top();
			t.write(name + "/" + ttos(i++), *this);
			what.pop();
		}
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
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet("Def/" + name));
		dataset->read(ptr, typeToH5<SomeType>());
		delete dataset;
	}

	void read(String& what,
	          String name)
	{
		readInternal(what, name);
	}

	void read(bool& value, String name)
	{
		unsigned char tmp[1];
		tmp[0] = ' ';
		void* ptr = static_cast<void *>(tmp);
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet("Def/" + name));
		dataset->read(ptr, typeToH5<unsigned char>());
		delete dataset;
		value = (tmp[0] == '1');
	}

	template<typename T1, typename T2>
	void read(std::pair<T1, T2>& what,
	          String name)
	{
		read(what.first, name + "/0");
		read(what.second, name + "/1");
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
		SizeType x = 0;
		read(x, name);
		value = static_cast<SomeType>(x);
	}

	template<typename T>
	void read(std::vector<T>& what,
	          String name,
	          typename EnableIf<Loki::TypeTraits<T>::isArith, int>::Type = 0)
	{
		readInternal(what, name);
	}

	template<typename T>
	void read(std::vector<std::complex<T> >& what,
	          String name,
	          typename EnableIf<Loki::TypeTraits<T>::isArith, int>::Type = 0)
	{
		readInternal(what, name);
	}

	template<typename T1, typename T2>
	void read(std::vector<std::pair<T1, T2> >& what,
	          String name,
	          typename EnableIf<
	          Loki::TypeTraits<T1>::isArith && Loki::TypeTraits<T2>::isArith,
	          int>::Type = 0)
	{
		SizeType size = 0;
		read(size, name + "/Size");
		what.resize(size);
		for (SizeType i = 0; i < size; ++i)
			read(what[i], name + "/" + typeToString(i));
	}

	template<typename T>
	void read(std::vector<T>& what,
	          String name,
	          typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith,
	          int>::Type = 0)
	{
		SizeType size = 0;
		read(size, name + "/Size");
		what.resize(size);
		for (SizeType i = 0; i < size; ++i)
			what[i].read(name + "/" + ttos(i), *this);
	}

	template<typename T>
	void read(std::vector<std::vector<T> >& what,
	          String name)
	{
		SizeType size = 0;
		read(size, name + "/Size");
		what.resize(size);
		for (SizeType i = 0; i < size; ++i)
			read(what[i], name + "/" + ttos(i));
	}

	template<typename T>
	void read(std::stack<T>& what,
	          String name)
	{
		throw RuntimeError("Cannot read label " + name + " (stacks not supported yet)\n");
	}

	// read functions END

private:

	IoNgSerializer(const IoNgSerializer&);

	IoNgSerializer& operator=(const IoNgSerializer&);

	template<typename SomeVectorType>
	void readInternal(SomeVectorType& what, String name)
	{
		typedef typename Real<typename SomeVectorType::value_type>::Type UnderlyingType;

		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet("Def/" + name));
		const H5::DataSpace& dspace = dataset->getSpace();
		const int ndims = dspace.getSimpleExtentNdims();
		if (ndims != 1)
			throw RuntimeError("IoNgSerializer: problem reading vector<arith> (ndims)\n");

		hsize_t* dims = new hsize_t[ndims];
		dspace.getSimpleExtentDims(dims);

		if (dims[0] == 0)
			throw RuntimeError("IoNgSerializer: problem reading vector<arith> (dims)\n");

		what.resize(dims[0], 0);
		void* ptr = static_cast<void *>(&(what[0]));
		dataset->read(ptr, typeToH5<UnderlyingType>());
		delete[] dims;
		delete dataset;
	}

	void overwriteNotSupported(WriteMode mode)
	{
		if (mode == NO_OVERWRITE) return;
		throw RuntimeError("Overwrite not supported for this type\n");
	}

	template<typename SomeType>
	void overwrite(String name, const void* ptr)
	{
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet(name));
		dataset->write(ptr, typeToH5<SomeType>());
		delete dataset;
	}

	template<typename SomeType>
	void internalWrite(String name, const void* ptr, hsize_t dims[], SizeType ndims)
	{
		H5::DataSpace *dataspace = new H5::DataSpace(ndims, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		H5::DataSet* dataset = new H5::DataSet(hdf5file_->createDataSet(name,
		                                                                typeToH5<SomeType>(),
		                                                                *dataspace,
		                                                                dsCreatPlist));
		dataset->write(ptr, typeToH5<SomeType>());
		delete dataset;
		delete dataspace;
	}

	H5::H5File* hdf5file_;
};
}
#endif // IONGSERIALIZER_H
