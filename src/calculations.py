import md
import numpy as np
import pandas as pd
import initial_conditions as ic

params = ic.parameters()
dt = params[0]
time = params[3]
duration_time = time
duration_steps = int(time/dt)
# these are called and set as global variables 


def msd(r_):
    ''' this function calculates the mean squared deviation of the simulation 
    for each step

    Inputs
    ------
    r_: np.array(N+1, NP, 3)
        positions of all particles at each time step (without periodic boundary 
        conditions)

    Outputs
    -------
    _msd: np.array(N) 
        the mean squared deviation for each time step
    '''

    # Introduce a delay AND how many steps it will go calculate over.
    # Has to be longer than the relaxation period of the system. 

    params = ic.parameters()
    NP = params[2]
    N = int(params[3]/params[0])

    msd_ = np.zeros(N)
    dr = r_ - r_[0]

    for n in range(N):
        r_sum = 0
        for i in range(NP):
            r_sum += np.linalg.norm(dr[n, i, :])
        msd_[n] = r_sum 

    return msd_   


def vcf(v, t_start=0):
    ''' this function calculates the velocity autocorrelation function of the
    simulation for each step from an initial time for a specific duration 
    (the default is 3x the relaxation period).

    Inputs
    ------
    v: np.array(N+1, NP, 3)
        velocities of all particles at each time step
    t_start: float
        time at which the vcf begins being calculated. Must be greater than the
        equilibration time.

    Outputs
    -------
    _vcf: np.array(N) 
        the velocity autocorrelation function for each time step
    '''

    params = ic.parameters()
    dt = params[0]
    NP = params[2]

    step_start = int(t_start/dt)

    global duration_steps

    vcf_ = np.zeros(duration_steps)
    dv = np.zeros((NP,3))

    for n in range(step_start, step_start+duration_steps):
        v_sum = 0
        for i in range(NP):
            dv[i] = np.dot(v[n,i],v[step_start,i])
            v_sum += np.linalg.norm(dv[i])**2
        vcf_[n-step_start] = v_sum 

    return vcf_


def output_to_csv(array, vcf_or_msd, temp):
    ''' this function exports vcf values to a csv file for all repeats of a 
    simulation at a specific temperature

    Inputs
    ------
    array: np.array((repeats,N))
        each row of this array contains all of the vcf values for a single 
        repeat at a given temperature
    vcf_or_msd: str
        string 'vcf' or 'msd' used to name the output file correctly
    temp: float
        the temperature of the simulation run, also to be included in the file 
        name

    Outputs
    -------
    {vcf_or_msd}_T{temp}.csv: csv
        csv file containing all repeats of a simulation at a specific 
        temperature
    '''

    params = ic.parameters()
    dt = params[0]
    N = int(params[3]/params[0])
    
    global duration_steps
    
    t_axis = [i*dt for i in range(duration_steps)]

    data = {'time / s': t_axis}

    for i in range(len(array)):
        data[f'T={temp}_{i}'] = array[i]
    '''
    try replacing with:
    for vcf in array:
        data[f'T={temp}_{i}'] = vcf
    '''

    df = pd.DataFrame(data)
    df.to_csv(f'{vcf_or_msd}_T{temp}.csv')

    print(f'{vcf_or_msd}_T{temp}.csv created \n')


def multiple_vcfs(temps, repeats, time_to_relax, duration=3):
    ''' this function runs md simulations and outputs the vcf results to a csv.

    Ensure that the length of the simulation is long enough to include all 
    measurement that you want to make (i.e. t_eq + (repeats-1)*time_to_relax + 
    duration)

    Inputs
    ------
    temps: lst
        list of temperatures (floats) measured in LJs epsilon
    repeats: int
        number of vcf measurements made on a md run (at a specific temperature)
    time_to_relax: float
        the time required for the vcf to decay to zero in seconds
    duration: float
        a multiplier for time_to_relax to determine the duration for which the 
        vcf is measured for each repeat
    
    Outputs
    -------
    {vcf_or_msd}_T{temp}.csv: csv
        csv file containing all vcf vlues for all repeats for a simulation at 
        a specific temperature
    '''    

    params = ic.parameters()
    dt = params[0]
    time = params[3]
    t_eq = params[6]

    steps_to_relax = int(time_to_relax/dt)
    
    global duration_steps, duration_time
    duration_steps = duration * steps_to_relax
    duration_time = duration_steps*dt


    time_needed = float(t_eq + (repeats-1)*time_to_relax + duration_time)
    # print(t_eq, (repeats-1)*time_to_relax, duration_time)
    if time < time_needed :
        raise IndexError(f'Length of simulation (as defined in initial_conditions.parameters) is not long enough for all vcf calculations in calculations.multiple_vcfs. Increase length of simulation to {time_needed} seconds')

    # this creates a list of timings to start each vcf repeat. 
    # The first timing is the end of the equilibration period with each subsequent time separated by 1 relaxation time 
    initial_timings = [t_eq+(i*time_to_relax) for i in range(repeats)]

    for i in range(len(temps)):
        
        lst_vcfs = np.zeros((repeats,duration_steps))
        md_out = md.run_md(ic.lattice_density(0.8), ic.random_v(), temps[i])
        v = md_out[2]

        for j in range(repeats):
            
            lst_vcfs[j] = vcf(v, initial_timings[j])

        output_to_csv(lst_vcfs, 'vcf', temps[i])

    '''
    try replacing with:
    for temp in temps:
        md_out = md.run_md(ic.lattice_density(0.8), ic.random_v(), temp)
    '''


def vcf_pre_plot(temps, repeats, time_to_relax):
    ''' this function expands the current cvs to include all vcf derived values
    necessary for plotting (average of vcf values, baseline corrected values 
    and normalised values).

    Inputs
    ------
    {vcf_or_msd}_T{temp}.csv: csv
        csv file containing only times and raw vcf values for all repeats of a 
        simulation at a specific temperature
    
    Outputs
    -------
    {vcf_or_msd}_T{temp}.csv: csv
        csv file containing average of vcf values, baseline corrected values 
        and normalised values, ready to be plotted
    '''

    params = ic.parameters()
    dt = params[0]
    steps_to_relax = int(time_to_relax / dt)

    for i in range(len(temps)):

        df = pd.read_csv(f'vcf_T{temps[i]}.csv')

        # average over each of the repeats column
        df['average'] = df.loc[ : , f'T={temps[i]}_0':f'T={temps[i]}_{repeats-1}'].mean(axis='columns')
        # baseline correction: should be after relaxation period
        df['baseline_corr'] = df['average'] - df['average'].iloc[steps_to_relax : ].mean()
        # normalisation
        df['normalised'] = df['baseline_corr'] / df['baseline_corr'].max()
        # cumulative integral
        df['cum_integral'] = df['baseline_corr'].iloc[1:].cumsum() * dt

        # side note: it is equivalent to baseline correct and then average OR average and then baseline correct

        df.to_csv(f'vcf_T{temps[i]}.csv')
