#ifndef IONGSERIALIZER_H
#define IONGSERIALIZER_H
#include "H5Cpp.h"
#include "../Vector.h"
#include "../TypeToString.h"
#include "TypeToH5.h"
#include <cassert>
#include <stack>
#include "../Complex.h"

namespace PsimagLite {

class IoNgSerializer {

	typedef std::vector<unsigned char> VectorOfBoolInternalType;

	static const unsigned char CANARY_VALUE = 170;

public:

	/*
		H5F_ACC_TRUNC - Truncate file, if it already exists,
		erasing all data previously stored in the file.
		H5F_ACC_EXCL - Fail if file already exists. H5F_ACC_TRUNC
		and H5F_ACC_EXCL are mutually exclusive
		H5F_ACC_RDONLY - Open file as read-only, if it already exists, and fail, otherwise
		H5F_ACC_RDWR - Open file for read/write, if it already exists, and fail, otherwise
	*/

	enum WriteMode {NO_OVERWRITE, ALLOW_OVERWRITE};

	IoNgSerializer(String filename, unsigned int mode)
	    : hdf5file_(0), filename_(filename), mode_(mode)
	{
#ifdef NDEBUG
		H5::Exception::dontPrint();
#endif

		try {
			hdf5file_ = new H5::H5File(filename, mode);
		} catch(H5::Exception& e) {
			delete hdf5file_;
			hdf5file_ = 0;
			throw e;
		}

		if (mode == H5F_ACC_TRUNC) {
			try {
				createGroup("");
			} catch(H5::Exception& e) {
				filename_ = "";
				delete hdf5file_;
				hdf5file_ = 0;
				throw e;
			}
		}

		if (mode == H5F_ACC_RDONLY) {
			try {
				readCanary();
			} catch(H5::Exception& e) {
				filename_ = "";
				delete hdf5file_;
				hdf5file_ = 0;
				throw e;
			}
		}
	}

	~IoNgSerializer()
	{
		if (hdf5file_ && mode_ != H5F_ACC_RDONLY)
			writeCanary();

		filename_ = "";
		delete hdf5file_;
		hdf5file_ = 0;
	}

	void open(String filename,
	          unsigned int mode)
	{
		if (hdf5file_)
			throw RuntimeError("IoNgSerializer::open(): object already open\n");

		filename_ = filename;
		hdf5file_ = new H5::H5File(filename, mode);
		if (hdf5file_ && mode_ != H5F_ACC_RDONLY)
			writeCanary();
	}

	void close()
	{
		if (hdf5file_ && mode_ != H5F_ACC_RDONLY)
			writeCanary();

		hdf5file_->close();
		delete hdf5file_;
		hdf5file_ = 0;
		filename_ = "";
	}

	void flush()
	{
		if (hdf5file_ && mode_ != H5F_ACC_RDONLY)
			writeCanary();

		hdf5file_->flush(H5F_SCOPE_GLOBAL);
	}

	const String& filename() const
	{
		return filename_;
	}

	void createGroup(String group)
	{
		hdf5file_->createGroup("Def/" + group);
	}

	bool doesGroupExist(String groupName)
	{
		groupName = "Def/" + groupName;

		try         {
			H5::Group group = hdf5file_->openGroup(groupName.c_str());
			group.close();
		} catch (...) {
			return false;
		}

		return true;
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
	           WriteMode allowOverwrite = NO_OVERWRITE,
	           typename EnableIf<Loki::TypeTraits<T1>::isArith &&
	           Loki::TypeTraits<T2>::isArith, int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		createGroup(name2);
		write(name2 + "/0", what.first);
		write(name2 + "/1", what.second);
	}

	template<typename T1, typename T2>
	void write(String name2,
	           const std::pair<T1, T2>& what,
	           WriteMode allowOverwrite = NO_OVERWRITE,
	           typename EnableIf<Loki::TypeTraits<T1>::isArith &&
	           !Loki::TypeTraits<T2>::isArith, int>::Type = 0)
	{
		overwriteNotSupported(allowOverwrite);
		createGroup(name2);
		write(name2 + "/0", what.first);
		what.second.write(name2 + "/1", *this);
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

	void write(String name,
	           const std::vector<bool>& what,
	           WriteMode allowOverwrite = NO_OVERWRITE)
	{
		overwriteNotSupported(allowOverwrite);
		VectorOfBoolInternalType converted = convertFromBoolean(what);
		write(name, converted, allowOverwrite);
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
		if (what.size() == 0) return;

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
	           typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith
	           && !IsPairLike<T>::True, int>::Type = 0)
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
	           const std::vector<T>& what,
	           WriteMode allowOverwrite = NO_OVERWRITE,
	           typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith
	           && IsPairLike<T>::True, int>::Type = 0)
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
	          String name,
	          typename EnableIf<Loki::TypeTraits<T1>::isArith &&
	          Loki::TypeTraits<T2>::isArith, int>::Type = 0)
	{
		read(what.first, name + "/0");
		read(what.second, name + "/1");
	}

	template<typename T1, typename T2>
	void read(std::pair<T1, T2>& what,
	          String name,
	          typename EnableIf<Loki::TypeTraits<T1>::isArith &&
	          !Loki::TypeTraits<T2>::isArith, int>::Type = 0)
	{
		read(what.first, name + "/0");
		what.second.read(name + "/1", *this);
	}

	void read(std::vector<bool>& what,
	          String name)
	{
		VectorOfBoolInternalType original;
		read(original, name);
		convertToBoolean(what, original);
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
	          typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith
	          && !IsPairLike<T>::True, int>::Type = 0)
	{
		SizeType size = 0;
		read(size, name + "/Size");
		what.resize(size);
		for (SizeType i = 0; i < size; ++i)
			what[i].read(name + "/" + ttos(i), *this);
	}

	template<typename T>
	void read(std::vector<T>& what,
	          String name,
	          typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith
	          && IsPairLike<T>::True, int>::Type = 0)
	{
		SizeType size = 0;
		read(size, name + "/Size");
		what.resize(size);
		for (SizeType i = 0; i < size; ++i)
			read(what[i], name + "/" + ttos(i));
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
		SizeType x = 0;
		read(x, name + "/Size");
		for (SizeType i = 0; i < x; ++i) {
			T t;
			t.read(name + "/" + ttos(x - i - 1), *this);
			what.push(t);
		}
	}

	static void unitTest(std::vector<bool>& x)
	{
		VectorOfBoolInternalType y = convertFromBoolean(x);
		std::cout<<"convertFromBoolean\n";
		std::cout<<"Input ";
		for (SizeType i = 0; i < x.size(); ++i)
			std::cout<<x[i]<<" ";
		std::cout<<"\nOutput: ";
		for (SizeType i = 0; i < y.size(); ++i)
			std::cout<<static_cast<unsigned short int>(y[i])<<" ";
		std::cout<<"\n\nconvertToBoolean\n";
		x.clear();
		convertToBoolean(x, y);
		std::cout<<"Output ";
		for (SizeType i = 0; i < x.size(); ++i)
			std::cout<<x[i]<<" ";
		std::cout<<"\n";
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

		SizeType complexSize = getComplexSize<typename SomeVectorType::value_type>(dims[0]);

		what.resize(complexSize, 0);
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

	void writeCanary()
	{
		hsize_t dims[1];
		dims[0] = 1;
		static const String name = "/Def/Canary";
		H5::DataSpace *dataspace = new H5::DataSpace(1, dims); // create new dspace
		H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
		H5::DataSet* dataset = 0;

		try {
			dataset = new H5::DataSet(hdf5file_->openDataSet(name));
		} catch (H5::Exception&) {
			dataset = new H5::DataSet(hdf5file_->createDataSet(name,
			                                                   typeToH5<unsigned char>(),
			                                                   *dataspace,
			                                                   dsCreatPlist));
		}

		unsigned char c = CANARY_VALUE;
		dataset->write(&c, typeToH5<unsigned char>());
		delete dataset;
		delete dataspace;
	}

	void readCanary()
	{
		static const String name = "/Def/Canary";

		H5::DataSet* dataset = new H5::DataSet(hdf5file_->openDataSet(name));
		const H5::DataSpace& dspace = dataset->getSpace();
		const int ndims = dspace.getSimpleExtentNdims();
		if (ndims != 1)
			throw RuntimeError("IoNgSerializer: problem reading vector<arith> (ndims)\n");

		hsize_t* dims = new hsize_t[ndims];
		dspace.getSimpleExtentDims(dims);

		if (dims[0] == 0)
			throw RuntimeError("IoNgSerializer: problem reading vector<arith> (dims)\n");

		unsigned char c = 0;
		void* ptr = static_cast<void *>(&c);
		dataset->read(ptr, typeToH5<unsigned char>());
		delete[] dims;
		delete dataset;

		if (c != CANARY_VALUE)
			throw RuntimeError("File " + filename_ + " is not valid (dead canary)\n");
	}

	static VectorOfBoolInternalType convertFromBoolean(const std::vector<bool>& src)
	{
		typedef VectorOfBoolInternalType::value_type ValueType;
		SizeType total = src.size();

		if (total == 0)
			return VectorOfBoolInternalType(booleanEncodedSize_, 0);

		SizeType bytesNeeded = total/8;
		bytesNeeded += 5;

		VectorOfBoolInternalType c(bytesNeeded, 0);
		encodeBooleanSize(c, total);
		SizeType blockSize = sizeof(ValueType);

		ValueType mask = 1;
		SizeType j = booleanEncodedStart_;
		SizeType bytes = 0;
		for (SizeType i = 0; i < total; ++i) {
			assert(j < c.size());
			assert(mask > 0);
			if (src[i]) c[j] |= mask;
			mask <<= 1;
			if (i > 0 && ((i + 1) % 8 == 0)) ++bytes;
			if (bytes == blockSize) {
				bytes = 0;
				++j;
				mask = 1;
			}
		}

		return c;
	}

	static void convertToBoolean(std::vector<bool>& dest,
	                             const VectorOfBoolInternalType& x)
	{
		typedef VectorOfBoolInternalType::value_type ValueType;
		SizeType numberOfBits = sizeof(ValueType)*8*x.size();
		SizeType blockSize = sizeof(ValueType);

		SizeType encodedSize = decodeBooleanSize(x);
		assert(encodedSize <= numberOfBits);

		numberOfBits = encodedSize;

		dest.resize(numberOfBits);

		ValueType mask = 1;
		SizeType j = booleanEncodedStart_;
		SizeType bytes = 0;
		for (SizeType i = 0; i < numberOfBits; ++i) {
			assert(j < x.size());
			assert(mask > 0);
			dest[i] = (x[j] & mask);
			mask <<= 1;
			if (i > 0 && ((i + 1) % 8 == 0)) ++bytes;
			if (bytes == blockSize) {
				bytes = 0;
				++j;
				mask = 1;
			}
		}
	}

	static void encodeBooleanSize(VectorOfBoolInternalType& x, SizeType total)
	{
		static short int byteSize = 256;
		assert(x.size() >= booleanEncodedSize_);
		SizeType tmp = total;
		std::fill(x.begin(), x.begin() + booleanEncodedSize_, 0);
		for (SizeType i = 0; i < booleanEncodedSize_; ++i) {
			x[i] = (tmp % byteSize);
			tmp >>= 8;
			if (tmp == 0) break;
		}
	}

	static SizeType decodeBooleanSize(const VectorOfBoolInternalType& x)
	{
		static short int byteSize = 256;
		assert(x.size() >= booleanEncodedSize_);
		SizeType tmp = 0;
		SizeType level = 1;
		for (SizeType i = 0; i < booleanEncodedSize_; ++i) {
			tmp += x[i]*level;
			level *= byteSize;
		}

		return tmp;
	}

	template<typename T>
	static typename EnableIf<!IsComplexNumber<T>::True, hsize_t>::Type
	getComplexSize(hsize_t x)
	{
		return x;
	}

	template<typename T>
	static typename EnableIf<IsComplexNumber<T>::True, hsize_t>::Type
	getComplexSize(hsize_t x)
	{
		if (x&1)
			throw RuntimeError("FATAL: Complex vector with uneven fp size\n");

		return x/2;
	}

	H5::H5File* hdf5file_;
	String filename_;
	unsigned int mode_;
	static const SizeType booleanEncodedSize_ = 4;
	static const SizeType booleanEncodedStart_ = 4;
};
}
#endif // IONGSERIALIZER_H
