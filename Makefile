#Detecting OS and Architecture
UNAME_S := $(shell uname -s)
UNAME_P := $(shell uname -p)
UNAME_M := $(shell uname -m)

CCC = g++ -O3 -m64 -march=native

ifeq ($(UNAME_S), Linux)
	CCC+= -lrt
endif

CPP_EX = CVPCP

#Ejecuta instruccion CVPCP
all_cpp: $(CPP_EX)

# ------------------------------------------------------------
# Clean
clean :
	/bin/rm -rf $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp *.o

# ------------------------------------------------------------
# Run

CVPCP: CVPCP.o
	$(CCC) CVPCP.o -o CVPCP

CVPCP.o: main.cpp readInstance.h global.h IGVND.h getRSS.c getCPUTime.c
	$(CCC) -c main.cpp -o CVPCP.o

# Local Variables:
# mode: makefile
# End:
