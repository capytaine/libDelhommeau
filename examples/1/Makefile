BIN=test

run: $(BIN)
	./$(BIN)

$(BIN): test.f90
		gfortran $< -I../../lib -L../../lib -lDelhommeau -o $@ 

.PHONY: clean run

clean:
	rm -rf $(BIN)
