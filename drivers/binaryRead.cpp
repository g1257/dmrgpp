#include "IoBinary.h"
#include "Vector.h"

typedef std::pair<SizeType,SizeType>PairType;

int main(int argc,char *argv[])
{

	if (argc<2) return 1;
	String myfile(argv[1]);
	PsimagLite::IoBinary::In fin(myfile);

	while(true) {
		String label = fin.readNextLabel();
		if (label=="NOTFOUND") break;
		std::cout<<"label="<<label<<"\n";
		char check = 0;
		SizeType total = 0;
		PairType type;
		fin.readCheckTotalAndType(check,total,type);
		int check1 = check;
		std::cout<<"check="<<check1<<" total="<<total<<" type="<<fin.nameOfType(type)<<"\n";
	}
}

