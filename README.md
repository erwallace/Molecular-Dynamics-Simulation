# Molecular-Dynamics-Simulation
A MS simulation using a Lennard-Jones potential and reduced units to determine the temperature dependence of mean squared deviation, velocity autocorrelation functions and diffusion constants. 
The run_md function in md.py was initially written as part of a university project, however as you can see I have since expanded this to fully automate all of the data analysis and plotting. The project report has also been included to give context to the calculations.

# Running Instructions
1. The first step is to find the relaxation period (t_relax) of the velocity autocorrelation function (VCF) at the lowest temperature at which you wish to simulate (as this will be the longest). Simulation parameters (initial_condtions.parameters()) need to be set accordingly. A simulation of approx. 1s should be more than sufficient.
2. The t_relax can then be found by running plotting.vcf_relaxation_period(), it is the time taken for the VCF to decay to a constant value.
3. Length of simulation will now need to be adjusted for subsequent calculations using the formula: t_eq + (n_vcf_measurements-1)*t_relax + duration. (t_eq = eqilibration period, n_vcf_measurements = number of VCF measurements at a specific temperature, t_relax = relaxation period of vcf, duration = length of each VCF measurement as a muliplier of t_relax).
4. Define variables: temperatures, n_vcf_measurments, t_relax in top_level.py and run the script. This will run the MD simulation(s), calculate all VCFs and plot them for you.

N.B. on my laptop, a simulation for 1 second with 172 particles took 130 seconds. I reccomend performing a similar calculation on your own system to approximate how long a given simualtion will take. MD simulations scale linearly with time and quadratically with the number of particles.
