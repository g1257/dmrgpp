#ifndef AINURSTATEMENT_H
#define AINURSTATEMENT_H
#include "AinurLexical.h"
#include "AinurStore.h"

namespace PsimagLite {

class AinurStatement {

public:

	typedef AinurLexical AinurLexicalType;
	typedef AinurLexicalType::VectorStringType VectorStringType;
	typedef Vector<Store>::Type VectorStoreType;

	AinurStatement(const String& s2)
	{
		String s = s2;
		AinurLexicalType::removeTrailingWhitespace(s);
		if (s == "") return;
		SizeType storageIndex = 0;
		VectorStringType dotified;

		VectorStringType leftAndRight;
		split(leftAndRight, s, "=");
		AinurLexicalType::removeTrailingWhitespace(leftAndRight);
		if (leftAndRight.size() != 1 && leftAndRight.size() != 2)
			err("Syntax error: " + s + "\n");

		procLeftEquality(storageIndex, dotified, leftAndRight[0], s);
		if (leftAndRight.size() == 1) return;

		if (storageIndex >= storage_.size())
			err("StorageIndex too big\n");

		unescape(leftAndRight[1]);
		storage_[storageIndex].procDotified(dotified, leftAndRight[1]);
	}

private:

	void procLeftEquality(SizeType& y,
	                      VectorStringType& dotified,
	                      String s,
	                      String context)
	{
		VectorStringType lhs;
		split(lhs,s," ");
		SizeType l = lhs.size();
		if (l == 0 || l > 2)
			err("Nothing or too much on left? " + context + "\n");
		int x = -1;

		if (l == 1) {
			split(dotified,lhs[0],".");
			if (dotified.size() == 0)
				err("Name too short " + context + "\n");
			x = storageIndexByName(dotified[0]);
		} else if (l == 2) {
			Store::Attribute attr = getAttribute(lhs[0], context);
			split(dotified,lhs[1],".");
			if (dotified.size() == 0)
				err("Name too short " + context + "\n");
			x = assignStorageByName(dotified[0]);
			storage_.push_back(Store(Store::UNKNOWN, Store::UNDEFINED, attr));
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

}; // class AinurStatement
} // namespace PsimagLite
#endif // AINURSTATEMENT_H
