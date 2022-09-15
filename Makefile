# Makefile for sunset_code
#
#
# ========================== OPTIONS ==============================================================
# -------------------------------------------------------------------------------------------------
# thermo     Isothermal (0) or thermal (1) flow                                        (default: 1)
# react      Reacting (1) or inert (0) flow                                            (default: 0)
# restart    Start from initial conditions (0) or restart file (1)                     (default: 0)
# hardinf    Non-reflecting inflow (0) or hard inflow (1)                              (default: 0)
# multispec  Single (0) or multispecies (1) flow                                       (default: 0)
# mpi        Shared only (0) or distributed-shared (1) acceleration                    (default: 0)          
# wisot      Adiabatic/prescribed heat flux (0) or isothermal (1) walls                (default: 0)
# dim3       Two (0) or three (1) dimensional simulation                               (default: 0)
# pgrad      Drive the flow with a pressure gradient and P.I.D control                 (default: 0)
# tdtp       Temperature dependent transport properties (1) or constant (0)            (default: 0)
# yout       Output the complete composition (1) or don't (0)                          (default: 1)
# -------------------------------------------------------------------------------------------------
#
# EXAMPLE USAGE:
# make thermo=0 react=0 mpi=1 etc...
#
# Choose compiler depending on whether mpi
ifeq ($(mpi),1)
FC := mpifort
LD := mpifort
else
FC := gfortran
LD := gfortran
endif

# Set compiler flags based on make options
CFLAGS := -Wall -O3 -g -m64
FFLAGS := -fopenmp -fbounds-check -ffpe-trap=zero -O3 -Wall -g -J./obj -I./obj -m64

# Isothermal or not. tdtp can only happen if not(isoT)
ifeq ($(thermo), 0)
FFLAGS += -DisoT
else
ifeq ($(tdtp),1)
FFLAGS += -Dtdtp
endif
endif
# Multi-species?
ifeq ($(multispec), 1)
FFLAGS += -Dms
endif
# Reacting, and if so, force multispeces
ifeq ($(react), 1)
FFLAGS += -Dreact
ifeq ($multispec),0)
FFLAGS += -Dms
endif
endif
# Restart from dump file.
ifeq ($(restart), 1)
FFLAGS += -Drestart
endif
# Inflow boundary types
ifeq ($(hardinf), 1)
FFLAGS += -Dhardinf
endif
# Multiprocessor? (use mpi?)
ifeq ($(mpi),1)
FFLAGS += -Dmp
endif
# Wall boundary types
ifneq ($(wisot),0)
FFLAGS += -Dwall_isoT
endif
# Three dimensional?
ifeq ($(dim3),1)
FFLAGS += -Ddim3
endif
# Flow driven by pressure gradient?
ifeq ($(pgrad),1)
FFLAGS += -Dpgrad
endif
# Output full chemical composition?
ifneq ($(yout),0)
FFLAGS += -Doutput_composition
endif
LDFLAGS := -fopenmp -m64 -lopenblas 

# Identify directories
SUB_DIRS := common base
SRC_DIR  := $(addprefix source/,$(SUB_DIRS))

# identify object files
#parameters come first, as almost everything depends on them.
OBJ_FILES := obj/kind_parameters.o obj/common_parameter.o obj/common_vars.o
OBJ_FILES += obj/rbfs.o obj/mirror_boundaries.o obj/derivatives.o 
OBJ_FILES += obj/mpi_transfers.o obj/thermodynamics.o
OBJ_FILES += obj/neighbours.o obj/output.o obj/setup.o
OBJ_FILES += obj/labf.o obj/fd.o obj/chemistry.o
OBJ_FILES += obj/characteristic_boundaries.o obj/rhs.o
OBJ_FILES += obj/step.o
OBJ_FILES += $(foreach sdir,$(SRC_DIR),$(patsubst $(sdir)/%.F90,obj/%.o,$(wildcard $(sdir)/*.F90)))
OBJ_FILES += $(foreach sdir,$(SRC_DIR),$(patsubst $(sdir)/%.cpp,obj/%.o,$(wildcard $(sdir)/*.cpp)))

HDEPS := $(OBJ_FILES:.o=.d)

vpath %.F90 $(SRC_DIR)
vpath %.cpp $(SRC_DIR)

#-------
default: sunset
sunset: $(OBJ_FILES)
	$(LD) -o $@ $^ $(LDFLAGS)

obj/%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<

#-include $(HDEPS)

clean:
	rm -vf ./obj/*.o
	rm -vf ./obj/*.mod
	rm -vf ./sunset
	rm -rfv fort.*	
	rm -vf ./data_out/layer*
	rm -vf ./data_out/flame*
	rm -vf ./data_out/time.out
	rm -vf ./data_out/statistics/*.out
	rm -vf ./paraview_files/LAYER*

