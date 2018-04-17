CC          = icpc
CLINKER     = icpc

#CFLAGS      =   -Wall -O3 -march=pentium
CFLAGS      =   -Wall -O3 -xHost

#CFLAGS      = -i-fast -lm  
#LIBS        = -lm
DEPEND= makedepend

SRC        = PolymerMD.c ran_uniform.c system.c boxmuller.c readinput.c initialization.c write.c potential.c verletlist.c force.c kinetic.c scale.c gr.c velocityverlet.c nvt.c nph.c npt.c store.c bondorder.c
OBJS       = PolymerMD.o ran_uniform.o system.o boxmuller.o readinput.o initialization.o write.o potential.o verletlist.o force.o kinetic.o scale.o gr.o velocityverlet.o nvt.o nph.o npt.o store.o bondorder.o
EXECS      = PolymerMD

default: PolymerMD

all: $(EXECS)

PolymerMD:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o PolymerMD $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

PolymerMD.o: system.h ran_uniform.h 
ran_uniform.o: system.h ran_uniform.h
boxmuller.o: system.h ran_uniform.h
system.o: system.h
readinput.o: system.h ran_uniform.h
initialization.o: system.h ran_uniform.h
write.o: system.h
potential.o: system.h
verletlist.o: system.h
force.o: system.h
kinetic.o: system.h
scale.o: system.h
gr.o: system.h
nvt.o: system.h
nph.o: system.h
npt.o: system.h
store.o: system.h
velocityverlet.o: system.h
bondorder.o: system.h
