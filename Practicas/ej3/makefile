FC = gfortran
CFLAGS = -g
LFLAGS = -o
LIBS = ./../SPARSKIT2/MATGEN/MISC/zlatev.f ./../SPARSKIT2/FORMATS/formats.f ./../SPARSKIT2/INFO/dinfo13.f ./../SPARSKIT2/libskit.a

all:
	$(FC) ej3.f $(LIBS) $(LFLAGS) ej3.exe
debug:
	$(FC) ej3.f $(LIBS) $(CFLAGS) $(LFLAGS) ej3.exe
clean:
	rm *.exe
