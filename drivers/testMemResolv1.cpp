#include <iostream>
#include <cstdlib>
#include "TestMemResolv1.h"
#include <unistd.h>

typedef TestMemResolv1<int> TestMemResolv1Type;

void recover()
{
	std::cout<<"From disk\n";
	PsimagLite::String filename = "test10.mem";
	PsimagLite::MemResolv mresolv(filename,"TestMemResolv1<int>");
	std::cout<<mresolv;

	TestMemResolv1Type* ptr = reinterpret_cast<TestMemResolv1Type*>(mresolv.get());
	std::cout<<"From disk: "<<ptr->get(1)<<"\n";
}

void save(TestMemResolv1Type* ptr)
{
	PsimagLite::MemResolv mresolv(ptr);
	void* ptrCopy = mresolv.dup();
	std::cout<<"copy= "<<ptrCopy<<"\n";

	TestMemResolv1Type* TestMemResolv1 = (TestMemResolv1Type*) ptrCopy;
	std::cout<<TestMemResolv1->get(1)<<"\n";

	mresolv.save("test10.mem","TestMemResolv1<int>");
}

int main(int argc, char *argv[])
{
	if (argc < 2) return 1;
	int n = atoi(argv[1]);
	if (n < 2) return 1;

	TestMemResolv1Type* ptr = new TestMemResolv1Type(n);
	TestMemResolv1Type& mytest = *ptr;
	mytest.setTo(10);
	std::cout<<mytest.get(1)<<"\n";

	save(ptr);
	delete ptr;
	sleep(1);
	recover();
}

