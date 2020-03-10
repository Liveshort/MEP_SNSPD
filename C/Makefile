CC     = gcc
FC     = gfortran
CFLAGS = -I -Wall -Wextra

DEPDIR = include
DEPDIR2= include_lapack
ODIR   = obj
LIBDIR = lib

_DEPS  = snspd.h thermal.h helper.h transmission.h yang.h yang_par.h waterfall_2s_res.h waterfall_3s_res.h types.h linalg.h
DEPS   = $(patsubst %,$(DEPDIR)/%,$(_DEPS))

_OBJ   = types.o linalg.o helper.o transmission.o thermal.o yang.o yang_par.o waterfall_2s_res.o waterfall_3s_res.o snspd.o simulation.o
OBJ    = $(patsubst %,$(ODIR)/%,$(_OBJ))

BLSLB  = $(LIBDIR)/librefblas.a
LPKLB  = $(LIBDIR)/liblapack.a
LPKELB = $(LIBDIR)/liblapacke.a
LIBS   = -lm $(LPKELB) $(LPKLB) $(BLSLB)

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -g -I $(DEPDIR) -I $(DEPDIR2) -std=c11 -c -fPIC -Ofast -o $@ $< $(CFLAGS)

all: simulation

simulation: $(OBJ) $(LIBS)
	$(FC) -g $(FFLAGS) $(LDFLAGS) -o bld/$@ $^

run:
	./bld/simulation ${args}

.PHONY: clean all run

clean:
	rm -f $(ODIR)/*.o bld/snspd bld/snspd.so bld/simulation
