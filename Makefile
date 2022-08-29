CXX = $(shell if [ -z "`which g++`" ] ; then echo "clang++"; else echo "g++"; fi)
#CXX = clang++

CXX_FLAGS_COMMON = -Wall -Wextra -Wpedantic -Wctor-dtor-privacy -Wnon-virtual-dtor \
-Wold-style-cast -Woverloaded-virtual -Wsign-promo -Wfloat-equal  -Wcast-qual \
-Wconversion -Wzero-as-null-pointer-constant -Wextra-semi -Wsign-conversion
CXX_FLAGS_GCC = -Wduplicated-branches -Wduplicated-cond -Wshadow=compatible-local -Wlogical-op

CXX_FLAGS = $(shell if [ "${CXX}" == "g++" ] ; then echo ${CXX_FLAGS_COMMON} ${CXX_FLAGS_GCC}; else echo ${CXX_FLAGS_COMMON}; fi)

main:
	#Run 'make tests' for testing.

tests: test
	./test

# -fsanitize=address,undefined 
test: tests/tests.cpp inc/libQBD.hpp inc/base.hpp inc/stationary.hpp inc/transient.hpp
	$(CXX) $(CXX_FLAGS) tests/tests.cpp -lboost_unit_test_framework -std=c++11 -Og -fsanitize=address,undefined -ffp-contract=fast -march=native -o test

clean:
	rm test

# --bug-hunting
check:
	cppcheck --enable=all --bug-hunting --language=c++ --max-ctu-depth=100 --std=c++11 .