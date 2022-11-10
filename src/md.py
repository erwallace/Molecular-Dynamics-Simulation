import time
import numpy as np
import initial_conditions as ic
from pickle import FALSE, TRUE

np.random.seed(2)

def tail_correction():
    ''' this function gives an error correction for when potentials have a cut 
    off, as used in periodic boundary conditions

    Inputs
    ------
    none

    Outputs
    -------
    float 
        value of the tail correction
    '''
    
    # import parameters
    params = ic.parameters()
    L = params[1]
    NP = params[2]
    s = params[4]
    e = params[5]
    
    return 6.4*NP*e*np.pi*((1/5)*s**3*(L/2)**-5 - (1/11)*s**9*(L/2)**-11)

def run_md(initial_r, initial_v, T_0, animation = FALSE):
    ''' [explanation]
    
    Inputs
    ------
    initial_r: np.array(NP,3) of floats
        list of vectors in a np.array containing the initail positions of all particles.
    initial_v: np.array(NP,3) of floats 
        list of vectors in a np.array containing the initail velocities of all particles.
    T_0: float
        initial temperatrue of the simulation in units of e/kT
    animation: bool 
        NOT CURRENTLY SET UP. if TRUE an animation of the simulation will play once calucations
        are complete. 

    Outputs
    -------
    tuple containing the following

    r: np.array(N+1, NP, 3)
        positions of all particles at each time step (within periodic boundary conditions)
    r_: np.array(N+1, NP, 3)
        positions of all particles at each time step (without periodic boundary conditions)
    v: np.array(N+1, NP, 3)
        velocities of all particles at each time step
    f: np.array(N+1, NP, 3)
        forces acting on all particles at each time step
    E: np.array(N)
        total energy at each time step
    K: np.array(N)
        kinetic energy at each time step
    V: np.array(N)
        potential energy at each time step
    '''

    t_start = time.perf_counter()

    # import parameters
    params = ic.parameters()
    dt =    params[0]
    L =     int(params[1])
    NP =    params[2]
    t =     params[3]
    s =     params[4]
    e =     params[5]

    N = int(t / dt)

    # Allocating arrays for 2D problem: first axis - time. second axis - particle's number. third - coordinate
    v = np.zeros((N+1, NP, 3))
    r = np.zeros((N+1, NP, 3)) # uses periodic boundary conditions
    r_ = np.zeros((N+1, NP, 3)) # no periodic boundary conditions (needed for r_acf)
    f = np.zeros((N+1, NP, 3))

    E = np.zeros(N)
    K = np.zeros(N)
    V = np.zeros(N)

    ## Initial postions/velocities
    v[0] = initial_v

    # scale velocities to T_0
    K[0] = 0.5*(np.linalg.norm(v[0])**2)
    T_n = K[0] / ((3/2)*NP)
    v[0] = v[0]*np.sqrt(T_0/T_n)

    r[0] = initial_r
    r_[0] = r[0]

    def compute_forces(n):
        '''The function computes forces on each particle at time step n'''
        
        tail_corr = tail_correction()
        
        for i in range(NP):
            for j in range(NP):
                if i != j:
                    rij = r[n,i] - r[n,j]
                
                    # minimum image convention
                    for k in range(3):
                        if abs(rij[k]) > (L/2) and rij[k] > 0:
                            rij[k] = rij[k] - L
                        elif abs(rij[k]) > (L/2) and rij[k] < 0:
                            rij[k] = rij[k] + L

                    rij_abs = np.linalg.norm(rij)
                    f[n,i] += 6*e*( 2*(s**12)/(rij_abs**13) - (s**6)/(rij_abs**7) ) * (rij/(rij_abs))
                    V[n] += 2*e*(((s/rij_abs)**12)-((s/rij_abs)**6)) + tail_corr #half usual V as it will be double counted

    print('Running dynamics')

    ## Run dynamics:
    for n in range(N):
        compute_forces(n)

        # record thermodynamic quantities
        for i in range(NP):
            K[n] += 0.5*(np.linalg.norm(v[n,i])**2)
        E[n] = K[n] + V[n]
        
        # advance time step
        v[n+1] = v[n] + f[n] * dt #time is scaled by mass
        r[n+1] = r[n] + v[n+1] * dt
        r_[n+1] = r_[n] + v[n+1] * dt
        r[n+1] = r[n+1] % L # periodic boundary conditions - perfect tiling
        
        T_n = K[n] / ((3/2)*NP)
        v[n+1] = v[n+1]*np.sqrt(T_0/T_n) # scale velocities to give initial temp

    t_finish = time.perf_counter()
    t_total = t_finish - t_start
    minutes = t_total // 60
    seconds = t_total % 60

    print(f'Dynamics for T={T_0} finished after {minutes:.0f} mins {seconds:.0f} secs')

    return r, r_, v, f, E, K, V, T_0
