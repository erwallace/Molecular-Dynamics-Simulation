import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

import initial_conditions as ic
import md
import calculations as calc

# matplotlib presets
mpl.rcParams['axes.labelsize'] = 15
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['axes.linewidth'] = 1.2


def vcf_relaxation_period(temps):
    ''' this function produces a plot of vcf vs. time to indicate the relaxation 
    time of the vcf for a range of temperatures at a given set of conditions. 
    The relaxation time is approx. the amount of time between collisions and is 
    signified by the amount of time required for the vcf to decay to a constant 
    value (excluding noise).

    If the relaxation times vary for differing temperatures then the largest 
    value should be chosen.

    Inputs
    ------
    temps: lst
        list of temperatures (floats) measured in LJs epsilon
    
    Outputs
    -------
    vcf_relaxation_period.pdf: plot
        this is a preliminary plot to determine the relaxation time of the vcf
    '''

    params = ic.parameters()
    dt = params[0]
    N = int(params[3]/params[0])

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 3.75))
    ax.set_xlabel('Time / s')
    ax.set_ylabel('vcf')

    t_axis = [i*dt for i in range(N)]
    
    for i in range(len(temps)):
        md_out = md.run_md(ic.lattice_density(0.8), ic.random_v(), temps[i])
        v = md_out[2]
        vcf_ = calc.vcf(v)
        # no baseline correction (need to know relaxation time first), just see when it decays to constant value
        ax.plot(t_axis[1:], vcf_[1:], label=f'T={temps[i]}')  # excludes t=0.0s data

    '''
    try replacing with:
    for temp in temps:
        md_out = md.run_md(ic.lattice_density(0.8), ic.random_v(), temp)

        ax.plot(t_axis[1:], vcf_[1:], label=f'T={temp}')
    '''

    ax.legend()
    plt.savefig('vcf_relaxation_period.pdf')


def vcf_plot(temps, normalised=True, x_axis_length=0, additional_name=''):
    ''' this function produces a plot of (un)normalised vcf vs. time for a 
    variety of temperatures.

    Inputs
    ------
    temps: lst
        list of temperatures (floats) measured in LJs epsilon
    normalised: bool
        gives a choice of whether you the plot produced will have a normalised 
        (to 1) vcf values
    x_axis_length: float
        length of the plots x-axis in seconds. If not specified then all data
        will plotted
    additional_name: string
        allows any additional notes to be added to the name of the output file
    
    Outputs
    -------
    vcf_normalised.pdf or vcf_unnormalised.pdf: plot
        plot of vcf vs. time for a range of temperatures
    '''

    params = ic.parameters()
    dt = params[0]

    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.set_xlabel('Time / s')
    ax.set_ylabel('vcf')
    ax.axhline(0, color='black')

    for i in range(len(temps)):
        
        df = pd.read_csv(f'vcf_T{temps[i]}.csv')

        if x_axis_length != 0: 
            x_axis_steps = int(x_axis_length/dt)
        else:
            x_axis_steps = df.shape[0] 

        if normalised is True:
            df.iloc[1:x_axis_steps].plot(x='time / s', y='normalised', kind='line', label=f'T={temps[i]}', ax=ax) # ax=ax plots on same plot

            ax.legend()
            plt.savefig(f'vcf_normalised_{additional_name}.pdf')

        else:
            df.iloc[1:x_axis_steps].plot(x='time / s', y='baseline_corr', kind='line', label=f'T={temps[i]}', ax=ax)

            ax.legend()
            plt.savefig(f'vcf_unnormalised_{additional_name}.pdf')

def vcf_integral(temps):
    ''' this function produces a plot of the cumulative integral of the vcf vs.
    time for a variety of temperatures. The cumulative integral is calculated
    from the unnormalised vcf vs. time plot. 

    The value of each line tends towards the diffusion constant at that
    temperature as time tends to infinity.

    Inputs
    ------
    temps: lst
        list of temperatures (floats) measured in LJs epsilon
    
    Outputs
    -------
    vcf_cum_integral.pdf: plot
        plot of cumulative integral vs. time for a range of temperatures
    '''

    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.set_xlabel('Time / s')
    ax.set_ylabel('Cumulative Integral')

    for i in range(len(temps)):
        
        df = pd.read_csv(f'vcf_T{temps[i]}.csv')

        df.iloc[1:].plot(x='time / s', y='cum_integral', kind='line', label=f'T={temps[i]}', ax=ax) # ax=ax plots on same plot

        ax.legend()
        plt.savefig('vcf_cum_integral.pdf')
