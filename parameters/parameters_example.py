'''
Simple example illustrating the use of the parameter toolbox.

See https://pypi.org/project/parameters for detailed information.

(Tom Tetzlaff, t.tetzlaff@fz-juelich.de, 2024)

'''

import parameters
import copy
import hashlib
from pprint import pformat

#############################################################################
def parameter_set_list(P):
    '''
    Generates list of parameters sets from a parameter space, 
    and adds a unique identifier ('label') to each parameters set (e.g., for data file names).    
    
    Parameters:
    -----------
    
    P : parameters.ParameterSpace
        Parameter space.
    
    Returns:
    --------

    l : list(dict) 
        List of parameter sets.
    
    '''

    l=[]

    for z in P.iter_inner(copy=True):
        l+=[dict(z)]
                 
        ## add md5 checksum as label of parameter set 
        l[-1]['label'] = hashlib.md5(pformat(l[-1]).encode('utf-8')).hexdigest() 
        
    return l

#############################################################################
def example():

    P = parameters.ParameterSpace({})

    P['a'] = 3.0
    P['b'] = parameters.ParameterRange([1,2,3])
    P['c'] = parameters.ParameterRange(['bla', 'blu', 'blubb'])

    PL = parameter_set_list(P)

    return P,PL

#############################################################################

P,PL = example()

print()
print('Example parameter space:')
print(P)
print()
print('List of parameter sets:')
for par in PL:
    print(par) ## at this point, one could for example submit jobs to some queuing system for each parameter combinations defined in the parameter space
print()
