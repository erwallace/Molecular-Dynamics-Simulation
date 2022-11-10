import random
import numpy as np

def parameters():
    ''' this function contains all parameters required for the run_md function,
    allowing consistnecy between functions as all additional functions will
    initially call parameters to set their parameter values.
    
    Inputs
    ------
    none

    Outputs
    -------
    tuple containing the following:
    dt: float 
        time step between cycles of the MD in seconds
    L: int 
        length of the cube within which the MD runs
    NP: int 
        number of particles in the simulation
    time: float 
        length of the simulation in seconds 
    s: float 
        sigma value for Lennard-Jones potential
    e: float 
        epsilon value for Lennard-Jones potential
    t_eq: int
        length of equilibration period of simulation in seconds (before any 
        mesurements are taken)
    '''
    
    dt = 0.005  # parameters[0]: time step is 5 ms
    L = 6       # parameters[1]: box length
    NP = 172    # parameters[2]: number of particles
    time = 12.7 # parameters[3]: length of simulation/s
    s = 1.0     # parameters[4]: sigma for LJ potential
    e = 1.0     # parameters[5]: epsilon for LJ potential
    t_eq = 1.0  # parameters[6]: equilibration period/s

    return dt, L, NP, time, s, e, t_eq

def lattice_density(density):
    ''' this function creates a lattice regular 3 dimensional lattice and then 
    removes particles at random until the density is equal to the function 
    input.

    WARNING: this fucntion should not currently be used as you should not be 
    specifying density, L and NP as they are linked via density=NP/(L**3), 
    assuming mass = 1

    Inputs
    ------
    density: float
        must take a value 0 < dnesity <= 1. This is the desired density of the 
        initial particle configuration

    Outputs
    -------
    r_lattice: np.array(NP,3)
        initial positions of all particles
    '''

    # import parameters
    params = parameters()
    L = params[1]
    NP = params[2]

    # create a cubic lattice with length L and density 1
    axis = np.arange(1, L+1)
    r_lattice = np.meshgrid(axis, axis, axis) 
    r_lattice = np.array(r_lattice).T.reshape(-1, 3)

    n_particles = np.count_nonzero(r_lattice)/3
    n_par_actual = int(n_particles*density)

    # random list of particles to be removed form lattice
    to_rmv = [random.randint(0, NP) for i in range(int(n_particles-n_par_actual))]
    to_rmv.sort(reverse=True)

    # remove particles from lattice
    for i in range(len(to_rmv)):
        r_lattice = np.delete(r_lattice, i, 0)

    # check density
    n_particles = np.count_nonzero(r_lattice)/3
    density_actual = n_particles/(L**3)

    if  int(n_particles) != NP:
        raise ValueError(f'The number of particles in the initial lattice ({int(n_particles)}) does not correspond to the correct total number of particles (NP={NP})')
    else:
        print(f'Initial configuration density is {density_actual:.3f}. This corresponds to {int(n_particles)} particles.')
    
    return r_lattice

def random_r():
    
    ''' this function randomly assigns initial particle positions within a 
    cube of length L.

    Inputs
    ------
    none

    Outputs
    -------
    r_random: np.array(NP,3)
        ...     
    '''

    #import parameters
    params = parameters()
    L = params[1]
    NP = params[2]
    
    r_random = np.random.sample(NP*3) # selects without replacement (to prevent r_ij = 0)
    r_random = L*np.reshape(r_random, (NP, 3))

    return r_random

def random_v():
    ''' this function randomly generates (unscaled) initial velocities of all
    particles.

    Inputs
    ------
    none

    Outputs
    -------
    v_random: np.array(NP,3)
        ...
    '''
    
    # import parameters
    params = parameters()
    NP = params[2]

    vavg = [0]
    v = np.zeros((NP, 3))

    for i in range(NP):
        v[i] = np.random.uniform(-0.2,0.2,3)
        vavg += v[i]

    vavg = vavg/NP

    for i in range(NP): # this ensures that the velocity of the COM is 0
        v[i] += -vavg

    return v
