CXX = $(shell if [[ `command -v g++` == "" ]]; then echo clang++; else echo g++; fi)

main:
	#Run 'make tests' for testing. Run 'make install' for install.

tests: test
	./test

test: tests/tests.cpp inc/libQBD.hpp inc/base.hpp inc/stationary.hpp
	$(CXX) tests/tests.cpp -lboost_unit_test_framework -std=c++11 -Og -g -ffast-math -march=native -o test
