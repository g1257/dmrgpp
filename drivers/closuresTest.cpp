#include "Vector.h"

int main(int argc, char** argv)
{
	int n = 4;
	std::vector<double> v1(n,1.0);
	std::vector<double> v2(n,1.1);
	std::vector<double> v3;

	// multiplication by scalar
	v3 <= v2*1.2;
	std::cout<<v3;

//	// plus
	v3 <= v1 + v2;
	std::cout<<v3;
	v3 <= v1 + 1.3*v2;
	std::cout<<v3;
	v3 <= 1.3*v2 + v1;
	std::cout<<v3;

//	// minus
//	v3 <= v1 - v2;
//	std::cout<<v3;
//	v3 <= v1 - 0.5*v2;
//	std::cout<<v3;
//	v3 <= 0.5*v2 - v1;
//	std::cout<<v3;
}

