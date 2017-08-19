#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include <iostream>
#include <vector>

void handler(int sig) {
  void *array[20];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 20);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

void baz() {
  std::vector<double> a (10, 0.0);

  std::cout<< a.at(11) << std::endl;

 int *foo = (int*)-1; // make a bad pointer
  printf("%d\n", *foo);       // causes segfault
}

void bar() { baz(); }
void foo() { bar(); }


int main(int argc, char **argv) {
  // signal(SIGTERM, handler); //termination request, sent to the program
  // signal(SIGSEGV, handler); //invalid memory access (segmentation fault)
  // signal(SIGINT , handler); //external interrupt, usually initiated by the user
  // signal(SIGILL , handler); //invalid program image, such as invalid instruction
  // signal(SIGABRT, handler); //abnormal termination condition, as is e.g. initiated by std::abort()
  // signal(SIGFPE , handler); //erroneous arithmetic operation such as divide by zero
  // signal(SIGSEGV, handler);   // install our handler
  // signal(vector::_M_range_check, handler);
  foo(); // this will call foo, bar, and baz.  baz segfaults.
}