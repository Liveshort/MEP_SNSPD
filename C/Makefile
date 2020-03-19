CC     = gcc
FC     = gfortran
CFLAGS = -I -Wall -Wextra

DEPDIR = include
OBJDIR = obj
LIBDIR = lib

CBLSLB = $(LIBDIR)/libcblas.a
BLSLB  = $(LIBDIR)/librefblas.a
LPKLB  = $(LIBDIR)/liblapack.a
LPKELB = $(LIBDIR)/liblapacke.a

OPNBLS = $(LIBDIR)/libopenblas.a

ifneq ("$(wildcard $(OPNBLS))","")
BLSDIR = include_openblas
LIBS   = $(OPNBLS) -lm -lpthread
else
BLSDIR = include_netlib
LIBS   = -lm $(LPKELB) $(LPKLB) $(BLSLB) $(CBLSLB)
endif

_DEPS  = snspd.h thermal.h helper.h transmission.h yang.h yang_par.h waterfall_2s_res.h waterfall_3s_res.h types.h linalg.h
DEPS   = $(patsubst %,$(DEPDIR)/%,$(_DEPS))

_OBJ   = types.o linalg.o helper.o transmission.o thermal.o yang.o yang_par.o waterfall_2s_res.o waterfall_3s_res.o snspd.o simulation.o
OBJ    = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

$(OBJDIR)/%.o: %.c $(DEPS)
	$(CC) -g -I $(DEPDIR) -I $(BLSDIR) -std=c11 -c -fPIC -Ofast -o $@ $< $(CFLAGS)

all: simulation

simulation: $(OBJ) $(LIBS)
	$(FC) -g $(FFLAGS) $(LDFLAGS) -o bld/$@ $^

run:
	./bld/simulation ${args}

.PHONY: clean all run

clean:
	rm -f $(OBJDIR)/*.o bld/snspd bld/snspd.so bld/simulation
