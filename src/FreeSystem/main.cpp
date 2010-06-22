#include "FreeSystemCore.h"

using namespace Dmrg;

typedef double FieldType;
typedef FreeSystemCore<FieldType> FreeSystemCoreType;
typedef FreeSystemCoreType::HilbertVectorType HilbertVectorType;
typedef FreeSystemCoreType::FreeOperatorType FreeOperatorType;

int main()
{
	size_t n = 10;
	size_t dof = 2; // spin up and down
	
	psimag::Matrix<FieldType> t(n,n);
	for (size_t i=0;i<n;i++) {
		for (size_t j=0;j<n;j++) {
			if (i-j==1 || j-i==1) t(i,j) = 1.0;
		}
	}
	
	FreeSystemCoreType fsc(t,dof,true);
	std::vector<size_t> ne(dof,5); // 5 up and 5 down
	HilbertVectorType gs = fsc.newGroundState(ne);
	std::cout<<gs;
	
	size_t site = 0;
	size_t flavor = 0;
	FreeOperatorType myOp = fsc.newSimpleOperator("creation",site,flavor);
	
	HilbertVectorType phi = fsc.newState();
	
	myOp.apply(phi,gs);
	std::cout<<"-----------------\n";
	std::cout<<phi;
	
	FreeOperatorType myOp2 = fsc.newSimpleOperator("creation",site,flavor);
	HilbertVectorType phi2 = fsc.newState();
	myOp2.apply(phi2,phi);
	std::cout<<"-----------------\n";
	std::cout<<phi2;
}	

