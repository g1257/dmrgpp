#ifndef AINURREADABLE_H
#define AINURREADABLE_H
#include "AinurStore.h"

namespace PsimagLite {

class AinurReadable {

public:

	typedef Store StoreType;
	typedef Vector<StoreType>::Type VectorStoreType;
	typedef StoreType::AinurLexicalType AinurLexicalType;
	typedef AinurLexicalType::VectorStringType VectorStringType;
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;

	AinurReadable(const VectorStringType& names,
	              const VectorStoreType& storage)
	    : names_(names), storage_(storage)
	{}

	int storageIndexByName(String name) const
	{
		VectorStringType::const_iterator it = std::find(names_.begin(),
		                                                names_.end(),
		                                                name);
		if (it == names_.end())
			return -1;
		return it - names_.begin();
	}

	String& prefix() { return prefix_; }

	const String& prefix() const { return prefix_; }

	void printUnused(std::ostream& os) const
	{
		SizeType n = storage_.size();
		for (SizeType i = 0; i < n; ++i) {
			if (storage_[i].used() > 0 || storage_[i].valueSize() == 0)
				continue;
			assert(i < names_.size());
			os<<"WARNING: Unused label "<<names_[i]<<"\n";
		}
	}

	void readValue(long int& t, String s) const
	{
		int t2 = 0;
		readValue(t2, s);
		t = t2;
	}

	void readValue(SizeType& t, String s) const
	{
		int t2 = 0;
		readValue(t2, s);
		t = t2;
	}

	void readValue(int& t, String s) const
	{
		s = prefix_ + s;
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		const Store& store = storage_[x];
		if (store.type() != Store::SCALAR && store.subType() != Store::INTEGER)
			err("In input, " + s + " must be an integer\n");
		store.increaseUsage();
		t = atoi(store.value(0, names_[x]).c_str());
	}

	void readValue(RealType& t, String s) const
	{
		s = prefix_ + s;
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		const Store& store = storage_[x];
		if (store.type() != Store::SCALAR && store.subType() != Store::REAL)
			err("In input, " + s + " must be a real\n");
		store.increaseUsage();
		t = atof(store.value(0, names_[x]).c_str());
	}

	void readValue(String& t, String s) const
	{
		s = prefix_ + s;
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		const Store& store = storage_[x];
		if (store.type() != Store::SCALAR && store.subType() != Store::STRING)
			err("In input, " + s + " must be a string\n");
		store.increaseUsage();
		t = store.value(0, names_[x]);
	}

	// read vectors
	template<typename VectorLikeType>
	typename EnableIf<IsVectorLike<VectorLikeType>::True,void>::Type
	readValue(VectorLikeType& v, String s) const
	{
		s = prefix_ + s;
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		const Store& store = storage_[x];
		if (store.type() != Store::VECTOR)
			err("In input, " + s + " must be a vector\n");
		store.increaseUsage();
		SizeType n = store.valueSize();
		if (n == 0)
			err("In input, vector " + s + " has 0 entries\n");

		store.increaseUsage();

		String tmp = (n == 2) ? store.value(1, names_[x]) : "";
		AinurLexicalType::removeTrailingBlanks(tmp);
		SizeType start = tmp.find("...");
		SizeType times = (start != String::npos && tmp.length() > 3)  ?
		            atoi(tmp.substr(start + 3, tmp.length() - 3).c_str()) : 0;

		if (n == 2 && start != String::npos) {
			assert(static_cast<SizeType>(x) < names_.size());
			if (v.size() < 3 && times == 0)
				err("Ellipsis cannot be used for this vector, " + names_[x] + "\n");
			if (times > 0) v.resize(times);
			n = v.size();
			for (SizeType i = 0; i < n; ++i)
				getEntryFromString(v[i], store.value(0, names_[x]));
			return;
		}

		if (v.size() != n) {
			v.clear();
			v.resize(n);
		}

		for (SizeType i = 0; i < n; ++i)
			getEntryFromString(v[i], store.value(i, names_[x]));
	}

	// read matrices
	template<typename FloatingType>
	typename EnableIf<Loki::TypeTraits<FloatingType>::isFloat,void>::Type
	readValue(Matrix<FloatingType>& m, String s) const
	{
		s = prefix_ + s;
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		const Store& store = storage_[x];
		if (store.type() != Store::MATRIX)
			err("In input, " + s + " must be a matrix\n");
		store.increaseUsage();
		SizeType n = store.valueSize();
		if (n < 2)
			err("In input, matrix " + s + " internal storage error I\n");

		store.increaseUsage();

		SizeType rows = atoi(store.value(0, names_[x]).c_str());
		SizeType cols = atoi(store.value(1, names_[x]).c_str());

		m.clear();
		if (rows == 0 && cols == 0) return;
		if (rows*cols == 0)
			err("Matrix with one of {rows, cols} 0 must have both 0\n");

		m.resize(rows, cols);

		if (rows*cols +2 != n)
			err("In input, matrix " + s + " internal storage error II\n");

		for (SizeType i = 0; i < rows; ++i)
			for (SizeType j = 0; j < cols; ++j)
				getEntryFromString(m(i,j), store.value(i + j*rows + 2, names_[x]));
	}

private:

	void getEntryFromString(SizeType& entry, String s) const
	{
		entry = atoi(s.c_str());
	}

	void getEntryFromString(RealType& entry, String s) const
	{
		entry = atof(s.c_str());
	}

	void getEntryFromString(ComplexType& entry, String s) const
	{
		err("getEntryFromString not implemented for complex\n");
	}

	const VectorStringType& names_;
	const VectorStoreType& storage_;
	String prefix_;
}; // class AinurReadable

} // namespace PsimagLite
#endif // AINURREADABLE_H
