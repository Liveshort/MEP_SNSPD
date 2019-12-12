CC     = gcc
CFLAGS = -I -Wall -Wextra
DEPS   = snspd.h

_OBJ   = snspd.o
OBJ    = $(patsubst %,$(ODIR)/%,$(_OBJ))

ODIR   = obj

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -fPIC -o $@ $< $(CFLAGS)

snspd: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
	$(CC) -shared -o $@.so $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o snspd snspd.so
