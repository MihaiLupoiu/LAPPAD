FC = gfortran
CFLAGS = -g
LFLAGS = -o
LIBS = ./../SPARSKIT2/MATGEN/MISC/zlatev.f ./../SPARSKIT2/FORMATS/formats.f ./../SPARSKIT2/INFO/dinfo13.f ./../SPARSKIT2/BLASSM/matvec.f ./../SPARSKIT2/libskit.a

all:
	$(FC) checkMatrix.f $(LIBS) $(LFLAGS) check.exe
	$(FC) ej4_1.f $(LIBS) $(LFLAGS) ej4_1.exe
	$(FC) ej4_2.f $(LIBS) $(LFLAGS) ej4_2.exe
	$(FC) ej4_3.f $(LIBS) $(LFLAGS) ej4_3.exe
	$(FC) ej4_4.f $(LIBS) $(LFLAGS) ej4_4.exe
	$(FC) ej4_5.f $(LIBS) $(LFLAGS) ej4_5.exe
debug:
	$(FC) ej4.f $(LIBS) $(CFLAGS) $(LFLAGS) ej4.exe
clean:
	rm *.exe