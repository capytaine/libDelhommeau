
SRC=src/constants.f90\
src/Green_Rankine.f90\
src/Delhommeau_integrals.f90\
src/Green_wave.f90\
src/matrices.f90

OBJ=$(SRC:.f90=.o)

LIBDIR=lib/
DYNAMIC_LIB=$(LIBDIR)/libDelhommeau.so
STATIC_LIB=$(LIBDIR)/libDelhommeau.a

$(STATIC_LIB): $(OBJ)
	@mkdir -p $(LIBDIR)
	ar r $(STATIC_LIB) $(OBJ)

$(DYNAMIC_LIB): $(OBJ)
	@mkdir -p $(LIBDIR)
	gfortran -shared -fPIC -o $(DYNAMIC_LIB) $(OBJ)

%.o: %.f90
	@mkdir -p $(LIBDIR)
	gfortran -cpp -c $< -J$(LIBDIR) -o $@

.PHONY: clean

clean:
	rm -rf $(OBJ) $(LIBDIR)
