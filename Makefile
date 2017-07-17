# Autogeneration of dependancies: see http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#tldr

compiler     = g++
warn_flags   = -Og -Wall
run_flags    = -O3
flags        = $(warn_flags) -std=c++11

# objects      = Prad.o sharedFunctions.o ImpuritySpecies.o RateCoefficient.o SD1DData.o
objects      = Prad.o sharedFunctions.o ImpuritySpecies.o

run: $(objects)
	$(compiler) -o Prad $(objects)
	./Prad

%.o: %.cpp
	$(compiler) $(flags) -c $< -o $@

Prad.o: Prad.cpp sharedFunctions.hpp ImpuritySpecies.hpp
sharedFunctions.o: sharedFunctions.cpp sharedFunctions.hpp
ImpuritySpecies.o: ImpuritySpecies.cpp ImpuritySpecies.hpp
# RateCoefficient.o: RateCoefficient.cpp RateCoefficient.hpp
# SD1DData.o: SD1DData.cpp SD1DData.hpp