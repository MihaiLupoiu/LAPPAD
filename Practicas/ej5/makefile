FC = gfortran
CFLAGS = -g
LFLAGS = -o
LIBS = ./../SPARSKIT2/MATGEN/MISC/zlatev.f ./../SPARSKIT2/ITSOL/iters.f ./../SPARSKIT2/BLASSM/matvec.f ./../SPARSKIT2/FORMATS/formats.f ./../SPARSKIT2/libskit.a /usr/local/lib/blas.a ./../SPARSKIT2/ITSOL/itaux.f 

all:
	$(FC) riters.f $(LIBS) $(CFLAGS) $(LFLAGS) ej5.exe
clean:
	rm *.exe
