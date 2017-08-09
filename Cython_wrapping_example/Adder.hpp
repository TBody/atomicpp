#include <vector>
#include <string>

class Adder{
public:
	Adder();
	Adder(std::vector<double> Input);
	std::vector<double> ReturnVector();
	void PlusOne();
	void PlusTwo(); //See if can call PlusOne from Python without declaring
	void PlusVector(std::vector<double> vector_to_add);
	std::string Print();
	std::string sayHello();
private:
	std::vector<double> internal;
	std::string privatestring; //Don't declare this attribute to Cython -- see if it can be accessed
};