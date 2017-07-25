#ifndef AINURSTATEMENT_H
#define AINURSTATEMENT_H
#include "AinurLexical.h"
#include "AinurStore.h"
#include <complex>

namespace PsimagLite {

class AinurStatements {

public:

	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef AinurLexical AinurLexicalType;
	typedef AinurLexicalType::VectorStringType VectorStringType;
	typedef Vector<Store>::Type VectorStoreType;

	void push(const String& s2)
	{
		String s = s2;
		AinurLexicalType::removeTrailingWhitespace(s);
		if (s == "") return;
		SizeType storageIndex = 0;

		VectorStringType leftAndRight;
		split(leftAndRight, s, "=");
		AinurLexicalType::removeTrailingWhitespace(leftAndRight);
		if (leftAndRight.size() != 1 && leftAndRight.size() != 2)
			err("Syntax error: " + s + "\n");

		procLeftEquality(storageIndex, leftAndRight[0], s);
		if (leftAndRight.size() == 1) return;

		if (storageIndex >= storage_.size())
			err("StorageIndex too big\n");

		unescape(leftAndRight[1]);
		storage_[storageIndex].setRhs(leftAndRight[1]);
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
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		Store store = storage_[x];
		if (store.type() != Store::SCALAR && store.subType() != Store::INTEGER)
			err("In input, " + s + " must be an integer\n");
		store.increaseUsage();
		t = atoi(store.value(0).c_str());
	}

	void readValue(RealType& t, String s) const
	{
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		Store store = storage_[x];
		if (store.type() != Store::SCALAR && store.subType() != Store::REAL)
			err("In input, " + s + " must be a real\n");
		store.increaseUsage();
		t = atof(store.value(0).c_str());
	}

	void readValue(String& t, String s) const
	{
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		Store store = storage_[x];
		if (store.type() != Store::SCALAR && store.subType() != Store::STRING)
			err("In input, " + s + " must be a string\n");
		store.increaseUsage();
		t = store.value(0);
	}

	// read vectors
	template<typename VectorLikeType>
	typename EnableIf<IsVectorLike<VectorLikeType>::True,void>::Type
	readValue(VectorLikeType& v, String s) const
	{
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		Store store = storage_[x];
		if (store.type() != Store::VECTOR)
			err("In input, " + s + " must be a vector\n");
		store.increaseUsage();
		SizeType n = store.valueSize();
		if (n == 0)
			err("In input, vector " + s + " has 0 entries\n");

		v.clear();
		v.resize(n);
		for (SizeType i = 0; i < n; ++i)
			getEntryFromString(v[i], store.value(i));
		store.increaseUsage();
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
		err("Reading complex entry for vector unimplemented\n");
	}

	void procLeftEquality(SizeType& y,
	                      String s,
	                      String context)
	{
		VectorStringType dotified;
		VectorStringType lhs;
		split(lhs,s," ");
		SizeType l = lhs.size();
		// require vector.vector.integer FiniteLoops
		// matrix.integer FiniteLoops
		// identifier
		if (l < 1 || l > 3)
			err("Too much or too little on left? -- " + context + " --\n");
		int x = -1;

		if (l == 1) { // identifier
			split(dotified,lhs[0],".");
			if (dotified.size() != 1)
				err("Dotified failed " + context + "\n");
			x = storageIndexByName(dotified[0]);
		} else if (l == 2) { // matrix.integer FiniteLoops
			split(dotified,lhs[1],".");
			if (dotified.size() != 1)
				err("Dotified failed " + context + "\n");
			x = assignStorageByName(dotified[0]);
			storage_.push_back(Store(lhs[0]));
		} else if (l == 3) {
			// require vector.vector.integer FiniteLoops
			split(dotified,lhs[2],".");
			if (dotified.size() != 1)
				err("Dotified failed " + context + "\n");
			x = assignStorageByName(dotified[0]);
			storage_.push_back(Store(lhs[1], lhs[0]));
		}

		if (x < 0)
			err("Undeclared variable " + dotified[0] + "\n");

		y = x;
	}

	Store::Attribute getAttribute(String s, String context) const
	{
		if (s == "let" || s == "function") return Store::NONE;
		if (s == "require") return Store::REQUIRED;
		if (s == "const") return Store::CONST;

		err("Expected let require or const " + context + "\n");
		return Store::NONE;
	}

	int assignStorageByName(String name)
	{
		int x = storageIndexByName(name);
		if (x >= 0)
			err("Already in scope " + name + "\n");
		names_.push_back(name);
		return names_.size() - 1;
	}

	int storageIndexByName(String name) const
	{
		VectorStringType::const_iterator it = std::find(names_.begin(), names_.end(), name);
		if (it == names_.end())
			return -1;
		return it - names_.begin();
	}

	void unescape(String& s) const
	{

	}

	VectorStringType names_;
	VectorStoreType storage_;
}; // class AinurStatements
} // namespace PsimagLite
#endif // AINURSTATEMENT_H
