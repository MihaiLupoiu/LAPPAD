FC = gfortran
CFLAGS = -g
LFLAGS = -o
LIBS = ./../SPARSKIT2/MATGEN/MISC/zlatev.f ./../SPARSKIT2/FORMATS/formats.f ./../SPARSKIT2/INFO/dinfo13.f ./../SPARSKIT2/libskit.a

all:
	$(FC) ej2.f $(LIBS) $(LFLAGS) ej2.exe
debug:
	$(FC) ej2.f $(LIBS) $(CFLAGS) $(LFLAGS) ej2.exe
clean:
	rm *.exe
