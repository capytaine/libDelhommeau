all: clean run

DEBUG_COMPILE_OPTIONS=-g -Wall -Wextra -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
# COMPILE_OPTIONS=$(DEBUG_COMPILE_OPTIONS) -cpp #-fopenmp
COMPILE_OPTIONS=-cpp -fopenmp

LIBDIR=./lib
BINDIR=./bin

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

#####################
#  Minimal example  #
#####################
MINIMAL_SRC=examples/minimal/minimal_example.f90
MINIMAL_BIN=$(BINDIR)/minimal
$(MINIMAL_BIN): $(MINIMAL_SRC) $(OBJ)
	@mkdir -p $(BINDIR)
	gfortran $(COMPILE_OPTIONS) $(MINIMAL_SRC) -o $(OBJ) -I$(LIBDIR) -o $(MINIMAL_BIN)

#######################
#  Benchmark example  #
#######################
BENCHMARK_SRC=examples/benchmark/benchmark.f90
BENCHMARK_BIN=$(BINDIR)/benchmark
$(BENCHMARK_BIN): $(BENCHMARK_SRC) $(OBJ)
	@mkdir -p $(BINDIR)
	gfortran $(BENCHMARK_SRC) -o $(OBJ) -I$(LIBDIR) -L$(LIBDIR) $(COMPILE_OPTIONS) -o $(BENCHMARK_BIN)

run: $(MINIMAL_BIN) $(BENCHMARK_BIN)
	$(MINIMAL_BIN)
	$(BENCHMARK_BIN)

clean:
	rm -rf $(OBJ) $(LIBDIR) $(BINDIR)

.PHONY: all run clean
