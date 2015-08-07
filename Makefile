# for gfortran
#

FLAGSALWAYS =  -ffree-form
FLAGSOPT= -O3
FLAGSDEBUG= -g

F90=gfortran

# choose debugging or optimization
FLAGS=  ${FLAGSALWAYS}  ${FLAGSOPT}  #  change the last to ${FLAGSOPT} or ${FLAGSDEBUG}


#------------------------------------

%.o: %.f90
	$(F90) $(FLAGS) -c $< -o $@

%.o: %.F90
	$(F90) $(FLAGS) -c $< -o $@


robust: ROBUST.o
	${F90} ${FLAGS} EXAMPLE.F90 -o EXAMPLE.exe 

clean:
	/bin/rm -f *.o *.mod
