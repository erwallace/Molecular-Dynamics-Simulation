import initial_conditions as ic

def simulation_parameters(temps, repeats, time_to_relax):
    '''
    This function creates a file and outputs all variables used in the simulation for 
    future refernce.
    
    Inputs
    ------
    temps: lst
        list of temperatures (floats) measured in LJs epsilon
    repeats: int
        number of vcf measurements made on a md run (at a specific temperature)
    time_to_relax: float
        the time required for the vcf to decay to zero in seconds
    
    Outputs
    -------
    simulation_parameter.txt: txt file
        txt file containing all parameters used in this simulation.
    '''
    params = ic.parameters()

    with open('simulation_parameters.txt', 'w') as outfile:
        #outfile.write(f'')
        print(f'dt = {params[0]}',
              f'L = {params[1]}',
              f'NP = {params[2]}',
              f'time = {params[3]}', 
              f's = {params[4]}',
              f'e = {params[5]}',
              f't_eq = {params[6]}',
              f'temperatures = {temps}',
              f'no. vcf measurements per temp= {repeats}',
              f'relaxation period = {time_to_relax}', 
              file=outfile, sep='\n')
    
    outfile.close()
        
        