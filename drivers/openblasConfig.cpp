#include <iostream>
#include <openblas/cblas.h>
//#include <openblas/openblas_config.h>

int main()
{
	std::cout<<"openblas_get_parallel= "<<openblas_get_parallel()<<"\n";
	std::cout<<"openblas_get_config= "<<openblas_get_config()<<"\n";
}

