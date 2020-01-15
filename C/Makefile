CC     = gcc
CFLAGS = -I -Wall -Wextra

DEPDIR = include
DEPDIR2= include_lapack
ODIR   = obj
LIBDIR = lib

_DEPS  = snspd.h thermal.h helper.h yang.h yang_parallel.h two_stage_waterfall_res.h types.h linalg.h
DEPS   = $(patsubst %,$(DEPDIR)/%,$(_DEPS))

_OBJ   = types.o linalg.o helper.o thermal.o yang.o yang_parallel.o two_stage_waterfall_res.o snspd.o simulation.o
OBJ    = $(patsubst %,$(ODIR)/%,$(_OBJ))

BLSLB  = $(LIBDIR)/librefblas.a
LPKLB  = $(LIBDIR)/liblapack.a
LPKELB = $(LIBDIR)/liblapacke.a
LIBS   = -lm $(LPKELB) $(LPKLB) $(BLSLB)

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -g -I $(DEPDIR) -I $(DEPDIR2) -std=c11 -c -fPIC -Ofast -o $@ $< $(CFLAGS)

all: simulation

#snspd: $(OBJ)
#	$(CC) -shared -o bld/$@.so -fPIC $^ $(CFLAGS) $(LIBS)

#$(CC) -o bld/$@ $^ $(CFLAGS) $(LIBS)
simulation: $(OBJ) $(LIBS)
	$(FC) -g $(FFLAGS) $(LDFLAGS) -o bld/$@ $^

run:
	./bld/simulation ${args}

.PHONY: clean all run

clean:
	rm -f $(ODIR)/*.o bld/snspd bld/snspd.so bld/simulation
