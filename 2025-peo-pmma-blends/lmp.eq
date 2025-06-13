# VARIABLES
variable        data_name      index 	data.220
variable        settings_name  index 	sys.settings
variable        avg_freq       index    1000
variable        coords_freq    index    1000
variable        thermo_freq    index    1000
variable        dump4avg       index    100
variable        Tf	       index    220
variable        vseed1      index 1
variable        vseed2      index 2

#===========================================================
# GENERAL PROCEDURES
#===========================================================
units		real	# g/mol, angstroms, fs, kcal/mol, K, atm, charge*angstrom
dimension	3	# 3 dimensional simulation
newton		off	# use Newton's 3rd law
boundary	p p p	# periodic boundary conditions 
atom_style	full    # molecular + charge

#===========================================================
# FORCE FIELD DEFINITION
#===========================================================
special_bonds lj/coul 0.0 0.0 0.5
pair_style     lj/cut/coul/long 12.0 12.0 # inner outer_LJ outer_Coul (cutoff values, see LAMMPS Doc)
bond_style     hybrid harmonic            
angle_style    hybrid harmonic           
dihedral_style hybrid fourier 
improper_style harmonic 
pair_modify     shift yes mix geometric
neigh_modify delay 0 every 1 check yes page 500000 one 50000

#===========================================================
# SETUP SIMULATIONS
#===========================================================
# READ IN COEFFICIENTS/COORDINATES/TOPOLOGY
read_data ${data_name} 
include   ${settings_name}
kspace_style   pppm 1e-4            # long-range electrostatics sum method

# SET RUN PARAMETERS
timestep	1.0		# fs
run_style	verlet 		# Velocity-Verlet integrator

# SET OUTPUTS
thermo_style    custom step temp vol density etotal pe ebond eangle edihed ecoul elong evdwl enthalpy press
thermo_modify   format float %14.6f
thermo ${thermo_freq}

# DECLARE RELEVANT OUTPUT VARIABLES
variable        my_step equal   step
variable        my_temp equal   temp
variable        my_rho  equal   density
variable        my_pe   equal   pe
variable        my_ebon equal   ebond
variable        my_eang equal   eangle
variable        my_edih equal   edihed
variable        my_evdw equal   evdwl
variable        my_eel  equal   (ecoul+elong)
variable        my_ent  equal   enthalpy
variable        my_P    equal   press
variable        my_vol  equal   vol

fix  averages all ave/time ${dump4avg} $(v_avg_freq/v_dump4avg) ${avg_freq} v_my_temp v_my_rho v_my_vol v_my_pe v_my_edih v_my_evdw v_my_eel v_my_ent v_my_P file thermo.avg

# CREATE COORDINATE DUMPS 
dump crds all custom ${coords_freq} all.lammpstrj id type mol x y z
dump_modify crds sort id

#===========================================================
# TRANSITION TO NOSE-HOOVER FOR ADDITIONAL NVT COOLING
#===========================================================
fix dynamics all npt temp ${Tf} ${Tf} 100.0 iso 1.0 1.0 2000.0
fix mom all momentum 100 linear 1 1 1 angular
run             10000000
write_data      data.${Tf}.2 pair ii
unfix dynamics
