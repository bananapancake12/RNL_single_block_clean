##################################################################
#                        Makefile for f90                         #
###################################################################

#BGP_SYS = /bgsys/drivers/ppcfloor/comm

########################     compiler     #########################

#CSD3
F90 = mpif90
F77 = mpif90

#WSL
# F90 = mpiifx
# F77 = mpiifx
#F90 = $(BGP_SYS)/bin/mpixlf90
#F77 = $(BGP_SYS)/bin/mpixlf77

########################  compiler flags  #########################

#F90FLAGS= -c
#F90FLAGS= -c
#F90FLAGS= -c
#F90FLAGS = -c -warn -CB -debug extended
# Debug flags (default for most files)
DEBUG_F90FLAGS = -c -O0 -g -check bounds -traceback  #-check uninit -check pointers#-fpe0  -warn all
DEBUG_F77FLAGS = -c
# FFT / performance-critical flags (used for FFT sources)
FFT_F90FLAGS   = -c
FFT_F77FLAGS   = -c
# Default active flags (can be overridden on make command line)
F90FLAGS = $(DEBUG_F90FLAGS)
F77FLAGS = $(DEBUG_F77FLAGS)
LFLAGS =
#PREP = scorep
#PREP = skin

########################  objects alpha   #########################

INIT    = .
SRCDIR  = $(INIT)
OBJ     = $(INIT)
OBJDIR  = $(OBJ)
CALCDIR = $(INIT)

OBJECTS = $(OBJ)/declaration.o\
		  $(OBJ)/tridLU_3D.o\
		  $(OBJ)/transpose.o\
		  $(OBJ)/error.o\
		  $(OBJ)/record_out.o\
		  $(OBJ)/FOU3D.o\
		  $(OBJ)/littleharsh_mod.o\
		  $(OBJ)/init_mod.o\
		  $(OBJ)/start.o\
          $(OBJ)/stats.o\
          $(OBJ)/spectra.o\
          $(OBJ)/rft_buff.o\
          $(OBJ)/cft_buff.o\
          $(OBJ)/littleharsh.o
# 		  $(OBJ)/sl_stats2.o\
# 		  $(OBJ)/inst_sl_stats.o\

########################      build       #########################

littleharsh :printmsgA $(OBJECTS)
	@echo
	@echo Linking...
	$(PREP) $(F90) -o $@ $(OBJECTS) $(LFLAGS)
	@echo
	@echo littleharsh built, congratulations.
	@echo

########################     compile      #########################

$(OBJDIR)/declaration.o : $(SRCDIR)/declaration.f90 $(SRCDIR)/makefile
	@echo compiling declaration.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/declaration.f90

$(OBJDIR)/tridLU_3D.o : $(SRCDIR)/tridLU_3D.f90 $(SRCDIR)/makefile
	@echo compiling tridLU_3D.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/tridLU_3D.f90

$(OBJDIR)/transpose.o : $(SRCDIR)/transpose.f90 $(SRCDIR)/makefile
	@echo compiling transpose.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/transpose.f90

$(OBJDIR)/error.o : $(SRCDIR)/error.f90 $(SRCDIR)/makefile
	@echo compiling error.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/error.f90

$(OBJDIR)/FOU3D.o : $(SRCDIR)/FOU3D.f90 $(SRCDIR)/makefile
	@echo compiling FOU3D.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/FOU3D.f90

$(OBJDIR)/littleharsh_mod.o : $(SRCDIR)/littleharsh_mod.f90 $(OBJDIR)/declaration.o
	@echo compiling littleharsh_mod.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(OBJDIR) $(SRCDIR)/littleharsh_mod.f90

$(OBJDIR)/init_mod.o : $(SRCDIR)/init_mod.f90 $(OBJDIR)/declaration.o
	@echo compiling init_mod.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(OBJDIR) $(SRCDIR)/init_mod.f90

$(OBJDIR)/start.o : $(SRCDIR)/start.f90 $(SRCDIR)/makefile
	@echo compiling start.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/start.f90

$(OBJDIR)/stats.o : $(SRCDIR)/stats.f90 $(SRCDIR)/makefile
	@echo compiling stats.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/stats.f90

# $(OBJDIR)/sl_stats2.o : $(SRCDIR)/sl_stats2.f90 $(SRCDIR)/makefile
# 	@echo compiling sl_stats2.f90
# 	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/sl_stats2.f90

# $(OBJDIR)/inst_sl_stats.o : $(SRCDIR)/inst_sl_stats.f90 $(SRCDIR)/makefile
# 	@echo compiling inst_sl_stats.f90
# 	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/inst_sl_stats.f90

$(OBJDIR)/spectra.o : $(SRCDIR)/spectra.f90 $(SRCDIR)/makefile
	@echo compiling spectra.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/spectra.f90


$(OBJDIR)/rft_buff.o : $(SRCDIR)/rft_buff.f $(SRCDIR)/makefile
	@echo compiling rft_buff.f
	@cd $(OBJDIR); $(PREP) $(F77) $(FFT_F77FLAGS) -I$(SRCDIR) $(SRCDIR)/rft_buff.f -o rft_buff.o

$(OBJDIR)/cft_buff.o : $(SRCDIR)/cft_buff.f $(SRCDIR)/makefile
	@echo compiling cft_buff.f
	@cd $(OBJDIR); $(PREP) $(F77) $(FFT_F77FLAGS) -I$(SRCDIR) $(SRCDIR)/cft_buff.f -o cft_buff.o

$(OBJDIR)/record_out.o : $(SRCDIR)/record_out.f90 $(SRCDIR)/makefile
	@echo compiling record_out.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/record_out.f90

$(OBJDIR)/littleharsh.o : $(SRCDIR)/littleharsh.f90  $(SRCDIR)/makefile
	@echo compiling littleharsh.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/littleharsh.f90

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
	@find . \( -name '*.o' -o -name '*.mod' -o -name '*.smod' \) -exec rm {} \;

########################   end of file    #########################