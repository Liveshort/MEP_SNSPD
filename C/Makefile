CC     = gcc
CFLAGS = -I -Wall -Wextra
DEPS   = snspd.h thermal.h helper.h yang.h types.h

_OBJ   = types.o thermal.o helper.o yang.o snspd.o simulation.o
OBJ    = $(patsubst %,$(ODIR)/%,$(_OBJ))

LIBS   = -lm

ODIR   = obj

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -fPIC -o $@ $< $(CFLAGS)

all: snspd simulation

snspd: $(OBJ)
	$(CC) -shared -o bld/$@.so $^ $(CFLAGS) $(LIBS)

simulation: $(OBJ)
	$(CC) -o bld/$@ $^ $(CFLAGS) $(LIBS)

run:
	./bld/simulation

.PHONY: clean all run

clean:
	rm -f $(ODIR)/*.o bld/snspd bld/snspd.so bld/simulation
