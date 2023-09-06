# Compiles main.f

CFT = gfortran
LDR = gfortran
FFLAGS  = -ffixed-line-length-136 -fdefault-real-8 -Wall -fno-automatic -Ofast
GLIB = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
COMMAND = willchand

.f.o :
	$(CFT) $(FFLAGS) $*.f90 -c

MAIN = willchand.o

SUBS = \
brakzero.o hunt.o polint.o zbrak.o common.o zline.o


$(COMMAND): $(MAIN) $(SUBS)
	$(LDR) $(FFLAGS) -o $(COMMAND) $(MAIN) $(SUBS) 
