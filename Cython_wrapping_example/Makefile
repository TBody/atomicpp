compiler     = g++
atomicpp     = atomicpp
warn_flags   = -Og -Wall
run_flags    = -O3
help_flags   = --enable-checking -v -da -Q

objects      = testAdder.o Adder.o

run: $(objects)
	$(compiler) -o testAdder $(objects)
	./testAdder

%.o: %.cpp
	$(compiler) $(flags) -c $< -o $@

testAdder.o: testAdder.cpp Adder.hpp
Adder.o: Adder.hpp Adder.cpp