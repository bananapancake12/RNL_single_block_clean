##################################################################
#                        Makefile for f90                         #
###################################################################

########################     compiler     #########################

F90 = ifort
F77 = ifort
#F90 = gfortran 
#F77 = gfortran 

########################  compiler flags  #########################

F90FLAGS= -c
#F90FLAGS = -c -warn -CB -debug extended
F77FLAGS= -c
LFLAGS =
#PREP = scalasca -inst
#PREP = skin

########################  objects alpha   #########################

INIT    = .
SRCDIR  = $(INIT)
OBJ     = $(INIT)
OBJDIR  = $(OBJ)
CALCDIR = $(INIT)

OBJECTS = $(OBJ)/declaration_unify.o\
          $(OBJ)/unify_spectra.o
          
########################      build       #########################

unify_spectra :printmsgA $(OBJECTS)
	@echo
	@echo Linking... 
	$(PREP) $(F90) -o $@ $(OBJECTS) $(LFLAGS)
	@echo 
	@echo unify_spectra built, congratulations.
	@echo 

########################     compile      #########################

$(OBJDIR)/declaration_unify.o : $(SRCDIR)/declaration_unify.f90 $(SRCDIR)/makefile
	@echo compiling declaration_unify.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/declaration_unify.f90

$(OBJDIR)/unify_spectra.o : $(SRCDIR)/unify_spectra.f90  $(SRCDIR)/makefile
	@echo compiling unify_spectra.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/unify_spectra.f90

########################      message     #########################

printmsgA :
	@echo
	@echo Building ...
	@echo F77 Compiler flags : $(F77FLAGS)
	@echo F90 Compiler flags : $(F90FLAGS)
	@echo Linker   flags : $(LFLAGS)
	@echo Prepend  flags : $(PREP)

########################      clean       #########################

clean:
	@find . \( -name '*.o' \) -exec rm {} \;

########################   end of file    #########################
