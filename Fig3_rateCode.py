"""
Python script for reproduce the rate code task from the Clopath et al. 2010
publication (Fig. 4 a). Network consists of ten, recurrent connected neurons.
Every neuron receives input from one neuron with a Poisson distributed
activity as stimulating input.
Poisson neurons firing with different frequencies (2Hz, 4Hz, 6Hz, ..)
as described in the original publication.
In the original publication, the resulting weights are averaged over 100s.
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from ANNarchy import *
setup(dt=1.0,seed=23456)
from network import *
from cmap import myCmap


# Global parameters
duration = 1000 #ms

# Populations
"""
The 'PoissonPopulation' object is used to create a population with a Poisson distributed
activity with 10 neurons and a initial firing rate of 100 Hz.
The population of the ten AdEx neurons with the neuron model is created with the
'Population' object.
"""
poisPop = PoissonPopulation(geometry=10,rates=100.0)
pop_Ten = Population(geometry=10,neuron=AdExNeuron, name="pop_Ten")

# Projections
"""
Every neuron in the 'PoissonPopulation' is with one neuron in the neuron
population connected. Every spike of a Poisson neuron leads to a spike
of the corresponding AdEx neuron.
"""
projInp_Ten = Projection(
    pre = poisPop,
    post= pop_Ten,
    target='Exc'
).connect_one_to_one(weights = 30.0)

# Create a Projection for the recurrent connections
projTen_Ten = Projection(
    pre= pop_Ten,
    post= pop_Ten,
    target= 'Exc',
    synapse= ffSyn
).connect_all_to_all(weights = 0.25,allow_self_connections=True)
# Set network parameter
projTen_Ten.wMax= 0.55#0.25

projTen_Ten.transmit=0.0
projTen_Ten.set_fix = 0.0 # use the homeostatic mechanisms in the LTD term

def run():
    # Compile the network with ANNarchy
    compile()
    # Number of repetitions
    repeats = 20
    w = np.zeros((repeats,10,10)) #
    print('Start rate code experiment')
    for r in range(repeats):
        # Repeat in 100 times to get the 100s as in Clopath et al. 2010
        for i in range(100):
            # Set the different firint rates
            poisPop.rates = np.linspace(20,2,10)
            # Simulate the network for 1000 ms
            simulate(duration)
            # Save the changes in the weights
            w[r,:,:] += projTen_Ten.w
            # Reset the network for the next 1000 ms
            reset()
        w[r,:,:] = w[r,:,:]/100


    img = np.ones((10,10))
    w = np.mean(w,axis=0)

    """
    Adapt the output like in Clopath et al. 2010.
    Depending on the connection, set another number to get another color.
    prepare a matrix of weights with different values for the different
    connections as mentioned in the Clopath et al., 2010 publication
    weak connections (< (2/3 of max. weight)) == 0
    strong unidirectional (> (2/3 of max. weight)) connections == 1.0
    strong bidirectional (> (2/3 of max. weight)) connections == 2.0
    """

    # Weak connection (under 2/3 of maximum weight) have the value = 0.0
    maxima = (np.nanmax(w)*2./3.)
    idx_b = np.where(w < maxima)
    img[idx_b] = 0.0

    # Strong biidirectional connections (> 2/3 of maximum weight) = 2.0
    idx_r = np.asarray(np.where(w >=maxima))
    for i in range(len(idx_r[0])):
        ix = (idx_r[0,i],idx_r[1,i])
        for j in range(len(idx_r[0])):
            ix2 = (idx_r[0,j],idx_r[1,j])
            if ix2 == (ix[1],ix[0]):
                img[ix[0],ix[1]] = 2.0
                img[ix[1],ix[0]] = 2.0

    # Strong unidirectional connections (> 2/3 of maximum weight) are every else

    # Set the main diagonal to nan
    for i in range(10):
        w[i,i] = np.nan
        img[i,i] = np.nan

    # Start plotting
    plt.figure()
    plt.imshow(img.T,interpolation='none',cmap= myCmap(),vmin=0,vmax=2)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('Neuron Post',fontsize=20)
    plt.ylabel('Neuron Pre',fontsize=20)
    plt.savefig('Fig3_rateCode.png',bbox_inches='tight')
    plt.show()
    print("done")

if __name__ == "__main__":
    run()