#include "Adder.hpp"

using namespace addns;

int main(){
	std::printf("Hi\n");
	std::vector<double> Input = {1,2,3,4,5};
	Adder A(Input);
	A.PlusOne();
	A.PlusTwo();
	std::printf("%s\n",A.Print().c_str());
	std::printf("%s\n",A.sayHello().c_str());
	twoInts B = A.returntwoInts();
	std::printf("Two ints: %d, %d\n", B.a, B.b);

	char *s = "Hello, World!";
	std::string newString = A.convertToStdString(s);
	std::printf("%s\n",newString.c_str());

}