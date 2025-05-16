# This is a convenience thing, on ARCHER it will default to GNU whereas locally to local (different mpi wrappers are used)
ifdef CRAYOS_VERSION
.DEFAULT_GOAL :=GNU
else
.DEFAULT_GOAL :=local
endif

FTN=ftn
FFLAGS=-O3 -J . -DDEF_MODEL=MODEL_MONC -DMODEL_MONC=4

local: FTN=mpif90
local: GNU

debug: FFLAGS = -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -J . -DDEF_MODEL=MODEL_MONC -DMODEL_MONC=4
debug: local

GNU: casim
Cray: casim
Intel: casim

%.o: %.F90
	$(FTN) -c -o $@ $< $(FFLAGS)

casim: src/variable_precision.o src/which_mode_to_use.o src/type_process.o src/type_aerosol.o src/thresholds.o src/special.o src/qsat_casim_func.o src/mphys_parameters.o src/process_routines.o src/mphys_switches.o src/mphys_die.o src/mphys_constants.o src/sweepout_rate.o src/passive_fields.o src/ventfac.o src/preconditioning.o src/lookup.o src/m3_incs.o src/sum_procs.o src/distributions.o src/snow_autoconversion.o src/gauss_4A_func.o src/derived_constants.o src/breakup.o src/autoconversion.o src/aggregation.o src/aerosol_routines.o src/sedimentation.o src/mphys_tidy.o src/ice_nucleation.o src/ice_multiplication.o src/ice_melting.o src/ice_deposition.o src/homogeneous_freezing.o src/graupel_wetgrowth.o src/graupel_embryo.o src/evaporation.o src/adjust_deposition.o src/activation.o src/ice_accretion.o src/condensation.o src/accretion.o src/micro_main.o src/initialize.o src/generic_diagnostic_variables.o src.cloud_fraction_dummy.o src/casim_parent.o src/casim_stph.o
	mkdir -p build
	mv src/*.o build/.
	mv *.mod build/.

clean:
	rm build/*
