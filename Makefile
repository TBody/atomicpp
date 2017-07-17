COMPILER = g++
WARN_FLAGS = -std=c++11 -O0 -g -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror -Wno-unused
RUN_TIME = -std=c++11 -O3

RUN: Prad.cpp
	$(COMPILER) $(RUN_TIME) Prad.cpp -c
	$(COMPILER) $(RUN_TIME) Prad.o -o Prad
	./Prad

Testing: Prad.cpp
	$(COMPILER) $(WARN_FLAGS) Prad.cpp -c
	$(COMPILER) $(WARN_FLAGS) Prad.o -o Prad
	./Prad

clean:
	rm -f Laplacian *.o
