#include "Adder.hpp"

int main(){
	std::printf("Hi\n");
	std::vector<double> Input = {1,2,3,4,5};
	Adder A(Input);
	A.PlusOne();
	std::printf("%s\n",A.Print().c_str());
}