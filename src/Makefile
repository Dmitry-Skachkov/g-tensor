Name= Diag1
F90= ifort
Library= dsyevj3
OBJS= $(Name).o $(Library).o


$(Name): ${OBJS}
	$(F90) -o ~/bin/$(Name)  ${OBJS}  

$(Name).o: $(Name).f90
	$(F90) -c $(Name).f90 -o $(Name).o 

$(Library).o: $(Library).f
	$(F90) -c $(Library).f -o $(Library).o 

clean: 
	rm -rf *.o *.mod ~/bin/$(Name)



