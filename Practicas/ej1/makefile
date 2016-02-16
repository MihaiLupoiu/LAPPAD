FC = gfortran
CFLAGS = -g
LFLAGS = -o
LIBS = ./../SPARSKIT2/MATGEN/MISC/zlatev.f ./../SPARSKIT2/FORMATS/formats.f ./../SPARSKIT2/libskit.a

all:
	$(FC) ej1.f $(LIBS) $(CFLAGS) $(LFLAGS) ej1.exe
clean:
	rm *.exe