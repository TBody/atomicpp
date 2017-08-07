# Autogeneration of dependancies: see http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#tldr

compiler     = g++
atomicpp     = atomicpp
warn_flags   = -Og -Wall
run_flags    = -O3
help_flags   =  --enable-checking -v -da -Q
# flags        = $(warn_flags) -std=c++11
flags        = -Og -std=c++11

objects      = Prad.o $(atomicpp)/sharedFunctions.o $(atomicpp)/ImpuritySpecies.o $(atomicpp)/RateCoefficient.o $(atomicpp)/SD1DData.o $(atomicpp)/RateEquations.o 

run: $(objects)
	$(compiler) -o Prad $(objects)
	./Prad

%.o: %.cpp
	$(compiler) $(flags) -c $< -o $@

Prad.o: Prad.cpp $(atomicpp)/sharedFunctions.hpp $(atomicpp)/ImpuritySpecies.hpp
$(atomicpp)sharedFunctions.o: $(atomicpp)/sharedFunctions.cpp $(atomicpp)/sharedFunctions.hpp
$(atomicpp)ImpuritySpecies.o: $(atomicpp)/ImpuritySpecies.cpp $(atomicpp)/ImpuritySpecies.hpp
$(atomicpp)/RateCoefficient.o: $(atomicpp)/RateCoefficient.cpp $(atomicpp)/RateCoefficient.hpp
$(atomicpp)/SD1DData.o: $(atomicpp)/SD1DData.cpp $(atomicpp)/SD1DData.hpp
$(atomicpp)/RateEquations.o: $(atomicpp)/RateEquations.cpp $(atomicpp)/RateEquations.hpp
