TOPSRCDIR = ..

BLASLIB      = $(TOPSRCDIR)/librefblas.a
CBLASLIB     = $(TOPSRCDIR)/libcblas.a
LAPACKLIB    = $(TOPSRCDIR)/liblapack.a
LAPACKELIB   = $(TOPSRCDIR)/liblapacke.a

.SUFFIXES: .c .o
.c.o:
	$(CC) $(CFLAGS) -I. -I../include_lapack -c -o $@ $<

.PHONY: all
all: xexample_DGESV_rowmajor

LIBRARIES = $(LAPACKELIB) $(LAPACKLIB) $(BLASLIB)

# Double Precision Examples
xexample_DGESV_rowmajor: lapacke_dgesv_row.o lapacke_example_aux.o $(LIBRARIES)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^
	./$@

.PHONY: clean cleanobj cleanexe
clean: cleanobj cleanexe
cleanobj:
	rm -f *.o
cleanexe:
	rm -f x*
