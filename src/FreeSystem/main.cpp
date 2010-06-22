#include "FreeSystemCore.h"

int main()
{
	size_t n = 2;
	typedef double FieldType;
	psimag::Matrix<FieldType> t(n,n);
	throw std::runtime_error("Needs to fill matrix\n");
	size_t ne = 5;
	Dmrg::FreeSystemCore<FieldType> fsc(t,ne);
}	

