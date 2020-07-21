
SRC=src/constants.f90\
src/Green_Rankine.f90\
src/Initialize_Green_wave.f90\
src/Green_wave.f90\
src/matrices.f90

OBJ=$(SRC:.f90=.o)

LIB=lib/libDelhommeau.so

MODDIR=./tmp/

$(LIB): $(OBJ)
	@mkdir -p lib/
	gfortran -shared -o $(LIB) $(OBJ)

%.o: %.f90
	@mkdir -p $(MODDIR)
	gfortran -cpp -c $< -J$(MODDIR) -o $@

.PHONY: clean

clean:
	rm -rf $(OBJ) $(MODDIR) $(LIB)
