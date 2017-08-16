compiler     = g++
# atomicpp     = atomicpp
warn_flags   = -g -O0 -Wall -fno-inline
run_flags    = -O3
# flags        = $(warn_flags) -std=c++11
flags        = -std=c++11 -O2
# flags        = -Og -std=c++11

objects      = alglibinternal.o fasttransforms.o solvers.o alglibmisc.o integration.o specialfunctions.o ap.o interpolation.o spline.o dataanalysis.o linalg.o statistics.o diffequations.o optimization.o

cpp_run: $(objects)
	$(compiler) $(flags) -o Prad $(objects)
	./Prad

cpp: $(objects)
	$(compiler) $(flags) -o Prad $(objects)

%.o: %.cpp
	$(compiler) $(flags) -c $< -o $@

clean:
	rm -f **/*.o

spline.o: spline.cpp
alglibinternal.o: alglibinternal.cpp alglibinternal.h
fasttransforms.o: fasttransforms.cpp fasttransforms.h
solvers.o: solvers.cpp solvers.h
alglibmisc.o: alglibmisc.cpp alglibmisc.h
integration.o: integration.cpp integration.h
specialfunctions.o: specialfunctions.cpp specialfunctions.h
ap.o: ap.cpp ap.h
interpolation.o: interpolation.cpp interpolation.h
dataanalysis.o: dataanalysis.cpp dataanalysis.h
linalg.o: linalg.cpp linalg.h
statistics.o: statistics.cpp statistics.h
diffequations.o: diffequations.cpp diffequations.h
optimization.o: optimization.cpp optimization.h

