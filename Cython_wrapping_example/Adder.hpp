#include <vector>
#include <string>

namespace addns{
struct twoInts{
	int a;
	int b;	
};

class Adder{
public:
	Adder();
	Adder(std::vector<double> Input);
	std::vector<double> ReturnVector();
	void PlusOne();
	void PlusTwo(); //See if can call PlusOne from Python without declaring
	void PlusVector(std::vector<double> vector_to_add);
	void stringIn(std::string Input);
	std::string stringOut();
	std::string Print();
	std::string sayHello();
	twoInts returntwoInts();
	std::string convertToStdString(char *pyStr);
private:
	std::vector<double> internal;
	std::string privatestring; //Don't declare this attribute to Cython -- see if it can be accessed
	twoInts structint;
	std::string inString;
};
}