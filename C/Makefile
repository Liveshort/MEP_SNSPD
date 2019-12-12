CC     = gcc
CFLAGS = -I -Wall -Wextra
DEPS   = snspd.h thermal.h

_OBJ   = snspd.o thermal.o simulation.o
OBJ    = $(patsubst %,$(ODIR)/%,$(_OBJ))

ODIR   = obj
SDIR   = src

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -fPIC -o $@ $< $(CFLAGS)

all: snspd simulation

snspd: $(OBJ)
	$(CC) -shared -o $@.so $^ $(CFLAGS)

simulation: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean all

clean:
	rm -f $(ODIR)/*.o snspd snspd.so
