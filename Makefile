FC=gfortran
FFLAGS=-c
EXECUTABLE=optim.out
SRCDIR=src
OBJDIR=objs

$(EXECUTABLE): $(OBJDIR)/funcmod.o $(OBJDIR)/optmod.o $(OBJDIR)/main.o
	$(FC) $(OBJDIR)/funcmod.o $(OBJDIR)/optmod.o $(OBJDIR)/main.o -o $(EXECUTABLE)

$(OBJDIR)/funcmod.o: $(SRCDIR)/funcmod.f90
	$(FC) $(FFLAGS) $(SRCDIR)/funcmod.f90 -J$(OBJDIR) -o $(OBJDIR)/funcmod.o

$(OBJDIR)/optmod.o: $(SRCDIR)/optmod.f90
	$(FC) $(FFLAGS) $(SRCDIR)/optmod.f90 -J$(OBJDIR) -o $(OBJDIR)/optmod.o

$(OBJDIR)/main.o: $(SRCDIR)/main.f90
	$(FC) $(FFLAGS) $(SRCDIR)/main.f90 -J$(OBJDIR) -o $(OBJDIR)/main.o

clean:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod $(EXECUTABLE)

reset:
	rm output/report.txt; rm output/hyperparams.txt