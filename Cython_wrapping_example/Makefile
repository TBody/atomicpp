compiler     = g++
atomicpp     = atomicpp
warn_flags   = -Og -Wall
run_flags    = -O3
help_flags   = --enable-checking -v -da -Q

objects      = testAdder.o Adder.o

run_python:
	python setup.py build_ext --inplace --verbose
	python python_test.py

clean:
	rm -f add.cpp
	rm -rf build
	rm -f *.so
	rm -f *.o

run: $(objects)
	$(compiler) -o testAdder $(objects)
	./testAdder

%.o: %.cpp
	$(compiler) $(flags) -c $< -o $@

testAdder.o: testAdder.cpp Adder.hpp
Adder.o: Adder.hpp Adder.cpp