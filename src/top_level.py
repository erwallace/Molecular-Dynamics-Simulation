import calculations as calc
import plotting
import time
import output

simulation_start = time.strftime("%H:%M:%S", time.localtime())
print(f'\nThis simulation started at {simulation_start} \n')
t_start = time.perf_counter()

# define variables
temps = [0.2, 0.4, 0.6, 0.8]  # units of LJs ε
n_vcf_measurements = 20
relax_time = 0.3

# run md sims for all temps and output the vcf of each to a csv
calc.multiple_vcfs(temps, n_vcf_measurements, relax_time, duration=20)

# calculate values for plotting from raw vcf values
calc.vcf_pre_plot(temps, n_vcf_measurements, relax_time)

# plot vcf graphs
plotting.vcf_plot(temps, normalised=True) 
plotting.vcf_plot(temps, normalised=False)
plotting.vcf_plot(temps, normalised=True, x_axis_length=2*relax_time, additional_name='(short_x_axis)') 
plotting.vcf_plot(temps, normalised=False, x_axis_length=2*relax_time, additional_name='(short_x_axis)') 
plotting.vcf_integral(temps)

# output variables used in this simulation
output.simulation_parameters(temps, n_vcf_measurements, relax_time)

t_finish = time.perf_counter()
t_total = t_finish - t_start
minutes = t_total // 60
seconds = t_total % 60

print(f'Set of simulations completed after {minutes:.0f} mins {seconds:.0f} secs \n')


