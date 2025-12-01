.SUFFIXES: .o .f90

FC=gfortran -Wall -pedantic -fbounds-check -fcheck=all -O0 -g -pg

SRC_modules = modmain.f90 modrand.f90
SRC_read1 = read_var.f90 read_input.f90 read_lat.f90 read_gvec.f90 read_sym.f90
SRC_read2 =  read_vzz.f90 read_fermi.f90 read_gap.f90 read_elkeig.f90 read_rand.f90
SRC_init = init.f90 fermi.f90 findkpt.f90 r3frac.f90 
SRC_eps = eps00_Gamma.f90 eps00_q.f90 eps_gg.f90 eps_gg_file.f90 eps_gg1.f90 eps_gg1_file.f90 inv_lfe.f90 istogram.f90 combine.f90
SRC_rand = montecarlo.f90 set_mesh_and_el.f90 fft_interpolation.f90 cft.f90 random_kp.f90 grid.f90 eps_gg1_rand.f90 eps_gg1_file_rand.f90 symmetrization.f90 istogram_rand.f90 eps_gg_rand.f90 eps_gg_file_rand.f90 combine_rand.f90 interpolation.f90
# SRC_sc = response_functions.f90 
SRC_main = bcsdiel.f90 

SRC = $(SRC_modules) $(SRC_main) $(SRC_read1) $(SRC_read2) $(SRC_init) $(SRC_eps) $(SRC_rand) $(SRC_sc)

EXE = BCSdiel
OBJ = $(SRC:.f90=.o)

.f90.o:
	$(FC) -c $<

all:    $(OBJ)
	$(FC) -o $(EXE) $(OBJ) -lblas -llapack

clean: 
	 rm -f *.o *.mod
