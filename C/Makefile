CC     = gcc
CFLAGS = -I -Wall -Wextra

DEPDIR = include
DEPDIR2= include_lapack
ODIR   = obj
LIBDIR = lib

_DEPS  = snspd.h thermal.h helper.h yang.h types.h linalg.h lapacke_example_aux.h
DEPS   = $(patsubst %,$(DEPDIR)/%,$(_DEPS))

_OBJ   = types.o linalg.o helper.o thermal.o yang.o snspd.o lapacke_example_aux.o simulation.o
OBJ    = $(patsubst %,$(ODIR)/%,$(_OBJ))

BLSLB  = $(LIBDIR)/librefblas.a
LPKLB  = $(LIBDIR)/liblapack.a
LPKELB = $(LIBDIR)/liblapacke.a
LIBS   = -lm $(LPKELB) $(LPKLB) $(BLSLB)

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -I $(DEPDIR) -I $(DEPDIR2) -c -fPIC -Ofast -o $@ $< $(CFLAGS)

all: simulation

#snspd: $(OBJ)
#	$(CC) -shared -o bld/$@.so -fPIC $^ $(CFLAGS) $(LIBS)

#$(CC) -o bld/$@ $^ $(CFLAGS) $(LIBS)
simulation: $(OBJ) $(LIBS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o bld/$@ $^

run:
	./bld/simulation

.PHONY: clean all run

clean:
	rm -f $(ODIR)/*.o bld/snspd bld/snspd.so bld/simulation
