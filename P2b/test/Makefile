CXX=mpicxx
INCDIR=../include
CXXFLAGS=-O3 -I$(INCDIR)

SRC =$(INCDIR)/tma.h
SRC+=$(wildcard $(INCDIR)/tma/*.h)
SRC+=$(wildcard $(INCDIR)/tma/*/*.h)

main.cpp: $(SRC)

test: main.cpp $(SRC)
	@$(CXX) $(CXXFLAGS) -o $@ $<
	
run: test
	@./test

clean:
	@rm -f test
