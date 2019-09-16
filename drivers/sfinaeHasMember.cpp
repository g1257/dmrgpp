#include <type_traits>
#include <iostream>

// SFINAE test
// source: https://stackoverflow.com/questions/257288/is-it-possible-to-write-a-template-to-check-for-a-functions-existence
template <typename T>
class has_helloworld
{
    typedef char one;
    struct two { char x[2]; };

    template <typename C> static one test( typeof(&C::helloworld) ) ;
    template <typename C> static two test(...);

public:
    enum { value = sizeof(test<T>(0)) == sizeof(char) };
};

class A {
public:

	int helloworld() const { return 42; }
};

class B {};

template<typename T>
typename std::enable_if<has_helloworld<T>::value, int>::type
f(const T& t)
{
	return t.helloworld();
}

int main()
{
	A a;
	std::cout<<f(a)<<"\n";
	//B b;
	//f(b);
}

