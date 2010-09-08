#include "Vector.h"
#include "AkimaSpline.h"

int main(int argc,char *argv[])
{
	typedef double FieldType;
	typedef Vector<FieldType> VectorType;
	VectorType x,s;
	readTwoColumnData(argv[1],x,s);
	size_t total = atoi(argv[2]);

	AkimaSpline<VectorType> akimaSpline(x,s);
	
	FieldType dx = (akimaSpline.end()-akimaSpline.start())/total
	for (size_t i=0;i<total;i++) {
		FieldType x = i*dx + akimaSpline.start();
		std::cout<<x<<" "<<akimaSpline(x)<<"\n";
	}
}



