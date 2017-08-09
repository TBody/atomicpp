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
private:
	std::vector<double> internal;
};