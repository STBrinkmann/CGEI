# Stand-alone tests of CGEI's C++ engines under sanitizers (no R required).
#
#   make release   optimised build (-O2), plain run
#   make asan      AddressSanitizer + UndefinedBehaviorSanitizer (gcc, libgomp)
#   make tsan      ThreadSanitizer (clang with the LLVM OpenMP runtime and its
#                  Archer tool, which tells TSan about OpenMP synchronisation)
#   make all       all of the above
#
# Ubuntu packages: g++, clang, libomp-dev, libclang-rt-dev

CXX     ?= g++
CLANGXX ?= clang++
ARCHER  ?= $(shell find /usr/lib/llvm-*/lib -name libarcher.so 2>/dev/null | head -n 1)
FLAGS    = -std=c++17 -Wall -Wextra -g
SRC      = engine_tests.cpp
HDRS     = ../../src/viewshed_engine.h ../../src/los_geometry.h ../../src/boxfilter.h \
           ../../src/natural_breaks.h

.PHONY: all release asan tsan clean
all: release asan tsan

release: $(SRC) $(HDRS)
	$(CXX) $(FLAGS) -O2 -fopenmp $(SRC) -o engine_tests_release
	./engine_tests_release

asan: $(SRC) $(HDRS)
	$(CXX) $(FLAGS) -O1 -fno-omit-frame-pointer -fsanitize=address,undefined \
	  -fno-sanitize-recover=undefined -fopenmp $(SRC) -o engine_tests_asan
	ASAN_OPTIONS=detect_leaks=1:abort_on_error=1 UBSAN_OPTIONS=print_stacktrace=1 ./engine_tests_asan

tsan: $(SRC) $(HDRS)
	$(CLANGXX) $(FLAGS) -O1 -fsanitize=thread -fopenmp $(SRC) -o engine_tests_tsan
	OMP_TOOL_LIBRARIES=$(ARCHER) TSAN_OPTIONS="ignore_noninstrumented_modules=1 halt_on_error=1" \
	  ./engine_tests_tsan

clean:
	rm -f engine_tests_release engine_tests_asan engine_tests_tsan
