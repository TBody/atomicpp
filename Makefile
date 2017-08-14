compiler     = g++
atomicpp     = atomicpp
warn_flags   = -g -O0 -Wall -fno-inline
run_flags    = -O3
# flags        = $(warn_flags) -std=c++11
flags        = -std=c++11
# flags        = -Og -std=c++11

objects      = Prad.o $(atomicpp)/sharedFunctions.o $(atomicpp)/ImpuritySpecies.o $(atomicpp)/RateCoefficient.o $(atomicpp)/BilinearSpline.o $(atomicpp)/SD1DData.o $(atomicpp)/RateEquations.o 

cpp_run: $(objects)
	$(compiler) $(flags) -o Prad $(objects)
	./Prad

cpp: $(objects)
	$(compiler) $(flags) -o Prad $(objects)

%.o: %.cpp
	$(compiler) $(flags) -c $< -o $@

py_run: $(objects) $(atomicpp)/setup.py $(atomicpp)/atomicpy.pyx
	cd atomicpp; python setup.py build_ext --inplace --verbose
	python Prad.py

py: $(objects) $(atomicpp)/setup.py $(atomicpp)/atomicpy.pyx
	cd atomicpp; python setup.py build_ext --inplace --verbose

clean:
	rm -f $(atomicpp)/atomicpy.cpp
	rm -rf $(atomicpp)/build
	rm -rf $(atomicpp)/__pycache__
	rm -f $(atomicpp)/*.so
	rm -f $(atomicpp)/*.o
	rm -f **/*.o
	rm -rf **/*.dSYM
	rm -f *.html
	rm -f Prad

Prad.o: Prad.cpp $(atomicpp)/sharedFunctions.hpp $(atomicpp)/ImpuritySpecies.hpp
$(atomicpp)sharedFunctions.o: $(atomicpp)/sharedFunctions.cpp $(atomicpp)/sharedFunctions.hpp
$(atomicpp)ImpuritySpecies.o: $(atomicpp)/ImpuritySpecies.cpp $(atomicpp)/ImpuritySpecies.hpp
$(atomicpp)/RateCoefficient.o: $(atomicpp)/RateCoefficient.cpp $(atomicpp)/RateCoefficient.hpp
$(atomicpp)/BilinearSpline.o: $(atomicpp)/BilinearSpline.cpp $(atomicpp)/BilinearSpline.hpp
$(atomicpp)/SD1DData.o: $(atomicpp)/SD1DData.cpp $(atomicpp)/SD1DData.hpp
$(atomicpp)/RateEquations.o: $(atomicpp)/RateEquations.cpp $(atomicpp)/RateEquations.hpp

