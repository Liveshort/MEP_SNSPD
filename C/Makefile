CC     = gcc
CFLAGS = -I -Wall -Wextra
DEPS   = snspd.h thermal.h helper.h bergen_standard.h types.h

_OBJ   = types.o thermal.o helper.o bergen_standard.o snspd.o simulation.o
OBJ    = $(patsubst %,$(ODIR)/%,$(_OBJ))

ODIR   = obj

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -fPIC -o $@ $< $(CFLAGS)

all: snspd simulation

snspd: $(OBJ)
	$(CC) -shared -o bld/$@.so $^ $(CFLAGS)

simulation: $(OBJ)
	$(CC) -o bld/$@ $^ $(CFLAGS)

run:
	./bld/simulation

.PHONY: clean all run

clean:
	rm -f $(ODIR)/*.o bld/snspd bld/snspd.so bld/simulation
