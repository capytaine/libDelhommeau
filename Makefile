all: run

DEBUG_COMPILE_OPTIONS=-g -Wall -Wextra -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
# COMPILE_OPTIONS=$(DEBUG_COMPILE_OPTIONS) -cpp #-fopenmp
COMPILE_OPTIONS=-cpp -fopenmp

LIBDIR=./lib

SRC=src/constants.f90\
src/Green_Rankine.f90\
src/Delhommeau_integrals.f90\
src/Green_wave.f90\
src/matrices.f90

OBJ=$(SRC:.f90=.o)

%.o: %.f90
	@mkdir -p $(LIBDIR)
	gfortran $(COMPILE_OPTIONS) -J$(LIBDIR) -c $< -o $@

# STATIC_LIB=$(LIBDIR)/libDelhommeau.a
# $(STATIC_LIB): $(OBJ)
# 	@mkdir -p $(LIBDIR)
# 	ar r $(STATIC_LIB) $(OBJ)
#
# DYNAMIC_LIB=$(LIBDIR)/libDelhommeau.so
# $(DYNAMIC_LIB): $(OBJ)
# 	@mkdir -p $(LIBDIR)
# 	gfortran -shared -fPIC $(COMPILE_OPTIONS) -o $(DYNAMIC_LIB) $(OBJ)

EXAMPLES_SRC=examples/minimal/minimal_example.f90
EXAMPLES_BIN=$(EXAMPLES_SRC:.f90=.bin)

BENCHMARKS_SRC=\
benchmarks/openmp/benchmark_omp.f90\
benchmarks/tabulations/benchmark_tabulation.f90

BENCHMARKS_BIN=$(BENCHMARKS_SRC:.f90=.bin)

%.bin: %.f90 $(OBJ)
	gfortran $(COMPILE_OPTIONS) -I$(LIBDIR) $^ -o $@

run :$(EXAMPLES_BIN)
	examples/minimal/minimal_example.bin
	# benchmarks/tabulations/benchmark_tabulation.bin

BENCHMARK_RESULTS_DIR = results/$(shell git rev-parse HEAD)
BENCHMARK_RESULT      = $(BENCHMARK_RESULTS_DIR)/omp.csv
$(BENCHMARK_RESULT): $(BENCHMARKS_BIN)
	benchmarks/openmp/benchmark_omp.bin
	mkdir -p $(BENCHMARK_RESULTS_DIR)
	mv benchmark_results.csv $(BENCHMARK_RESULT)
	benchmarks/openmp/read_output.py

benchmark: $(BENCHMARK_RESULT)

clean:
	rm -rf $(OBJ) $(LIBDIR) $(EXAMPLES_BIN) $(BENCHMARKS_BIN)

.PHONY: all run clean benchmark
