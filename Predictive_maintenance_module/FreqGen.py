"""
FreqGen is a module aiming to generate sets of frequency components that can be observed in induction machines.
The scope is to recieve minimal parameters and to return lists of possible frequencies.
"""

def grid_frequencies(
        fgrid:float,
        N:int
        )->list:
    """
    Returns a list of N grid harmonics
    fgrid: Grid frequency
    N: Number of harmonics to consider
    """
    out = list()
    for i in range(1,N+1):
        out.append(i*fgrid)
    return out

def broken_bar_frequencies( 
    fgrid,N,s    
    )->list:    
    """
    Returns a list of 2*N frequency components
    associated to the sidebands of broken bars
    s: slip (per unit)
    N: Number of harmonics to consider
    """
    out = list()
    for i in range(1,N):
        x = s*i
        frecm = (1-2*x)*50
        frecp = (1+2*x)*50
        out.append(frecm)
        out.append(frecp)
    return out

# def 