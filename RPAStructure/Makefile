BUILDER = gfortran 
FFLAGS = -ffree-line-length-0 -Wall -O2
OPT = -fopenmp

#
F90 = $(BUILDER)

# Sources and utilities
SRC = 	asym_param.f90 asym_gauss.f90 asym.f90       \
	asym_io.f90 asym_cstes.f90 beta.f90 wi.f90 determinante.f90 main.f90 \
        asym_inm.f90

# define your objects
OBJS = 	asym_param.o asym_gauss.o asym.o       \
	asym_io.o asym_cstes.o beta.o wi.o determinante.o \
        asym_inm.o

OBJMAIN = main.o
MODS = 	asym_param.mod asym_cstes.mod asym_gauss.mod asym_io.mod asym.mod beta.mod wi.mod \
        determinante.mod asym_inm.mod main.mod

MAIN = main.so


# this rule defines how your python module is created


$(MAIN) : librpa.a $(OBJMAIN)
	f2py --fcompiler=$(F90) -I. -c -m main main.f90 librpa.a --f90flags=$(OPT) -lgomp
	f2py -h main.pyf main.f90 --overwrite-signature

librpa.a: $(OBJS)
	touch librpa.a
	\rm librpa.a
	ar rv librpa.a $(OBJS)



%.o: %.f90
	$(F90) $(FFLAGS) -o $@ -c $< $(OPT)




# Dependences

asym_cstes.o: asym_cstes.f90 
determinante.o: determinante.f90 asym_cstes.o
asym_io.o: asym_io.f90 asym_cstes.o
asym_param.o: asym_param.f90 asym_cstes.o asym_io.o
asym_inm.o : asym_cstes.o asym_param.o asym_io.o
asym_gauss.o: asym_gauss.f90 asym_cstes.o
beta.o : asym_gauss.o asym_cstes.o
wi.o : asym_gauss.o asym_cstes.o beta.o asym_param.o
asym.o: asym.f90 asym_cstes.o asym_gauss.o asym_io.o asym_param.o beta.o wi.o \
        determinante.o asym_inm.o
main.o: main.f90 asym.o asym_cstes.o asym_io.o
main.so: main.o main.f90

# house keeping

clean:
	rm -f *.a *.o *.mod *.d *.so *.pyf $(MAIN) $(MAIN2)

clean_all:
	rm -f *.a *.o *.mod *.so *.a *.pyf $(MAIN) $(MAIN2) *~ *.d set_pairing

