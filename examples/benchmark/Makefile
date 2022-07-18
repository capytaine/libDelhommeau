BIN=test

run: $(BIN)
	./$(BIN)

$(BIN): test.f90
		gfortran $< -I../../lib -L../../lib -lDelhommeau -fopenmp -o $@ 

.PHONY: clean run

clean:
	rm -rf $(BIN)
