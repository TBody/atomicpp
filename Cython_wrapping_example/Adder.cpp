#include "Adder.hpp"
#include <vector>
#include <string>
#include <cassert>

using namespace addns;

Adder::Adder():internal(1,0.0){};
Adder::Adder(std::vector<double> Input):internal(Input.size(),0.0){
	for(int i=0; i< internal.size(); ++i){
		internal[i] = Input[i];
	}
	privatestring = "Hello";
	structint.a = 1;
	structint.b = 3;
};
std::vector<double> Adder::ReturnVector(){
	return internal;
};
void Adder::PlusOne(){
	for(int i=0; i< internal.size(); ++i){
		internal[i] += 1;
	}
};
void Adder::PlusTwo(){
	PlusOne();
	PlusOne();
};
void Adder::PlusVector(std::vector<double> vector_to_add){
	assert( internal.size()== vector_to_add.size());
	for(int i=0; i< internal.size(); ++i){
		internal[i] += vector_to_add[i];
	}
};
std::string Adder::Print(){
	std::string printstring = "[";
	for(int i=0; i< internal.size(); ++i){
		printstring += (std::to_string(internal[i]) + ", ");
	}
	printstring += "]";
	return printstring;
};
void Adder::stringIn(std::string Input){
	inString = Input;
};
std::string Adder::stringOut(){
	return inString;
};
std::string Adder::sayHello(){
	return privatestring;
};
twoInts Adder::returntwoInts(){
	return structint;
}
std::string Adder::convertToStdString(char *pyStr){
	std::string cppStr(pyStr);
	return cppStr;
}