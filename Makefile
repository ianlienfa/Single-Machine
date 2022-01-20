# INC = /Library/gurobi912/mac64/include/
CPP = g++
CARGS = -std=c++14 -m64 -g -arch x86_64 
# CPPLIB = -L/Library/gurobi912/mac64/lib -lgurobi_c++ -lgurobi91
CPPLIB = 
INC = 

# all: seq_relax_c++ permSet_outtree seq_lagRelax_c++ seq_c++
all: permSet_outtree

clean: 
	rm *_c++ permSet_outtree

%_c++: %_c++.cpp
	$(CPP) $(CARGS) -o $@ $< -I$(INC) $(CPPLIB) -lm

permSet_outtree: permSet_outtree.cpp
	$(CPP) $(CARGS) -o $@ $< -I$(INC) $(CPPLIB) -lm

.PHONY: run
run:
	@echo "[Makefile] Running: permSet_outtree"
	@echo "----------------------------------------"
	./permSet_outtree < input.in
	@echo
