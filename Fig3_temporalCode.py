"""
Python script for reproduce the rate code task from the Clopath et al. 2010
publication (Fig. 4 b). Network consists of ten, recurrent connected neurons.
Every neuron receives input from one extra neuron as input to force to spike.
Use the SpikeSourceArray of ANNarchy to determine the spiking time points.
Every extern neuron spikes at an other time point.
In the original publication, the resulting weights are averaged over 100s.
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from ANNarchy import *
setup(dt=1.0,seed=20000) # 23456
from network import *
from cmap import myCmap


# Global parameters
duration = 200 #ms
# Time points for spikes
spike_times =[[1],[2],[3],[4],[5],[6],[7],[8],[9],[10]]

# Populations
"""
Use the SpikeSourceArray of ANNarchy to control the spiking time points of
the AdEx neurons."""
inpPop = SpikeSourceArray(spike_times=spike_times)
pop_Ten = Population(geometry=10,neuron=AdExNeuron, name="pop_Ten")

# Projections
"""
Create the one-to-one projections from the SpikeSourceArray population
to the AdEx population. If one neuron in the SpikeSourceArray population spikes,
1 ms later the corresponding AdEx neuron spikes.
"""
projInp_Ten = Projection(
    pre = inpPop,
    post=pop_Ten,
    target='Exc'
).connect_one_to_one(weights = 30.0)

# Create the projection for the recurrent connections
projTen_Ten = Projection(
    pre=pop_Ten,
    post=pop_Ten,
    target='Exc',
    synapse=ffSyn
).connect_all_to_all(weights = 0.5,allow_self_connections=True)
projTen_Ten.wMax= 0.75

projTen_Ten.set_fix = 0.0 # use the homeostatic mechanisms in the LTD term

def run():

    compile()
    # Repeat the experiments 1000 times, that the weights can be stable
    for i in range(1000):
        # Define the time points for spikes
        spkT_N1 = [0+(i*duration)]
        spkT_N2 = [20+(i*duration)]
        spkT_N3 = [40+(i*duration)]
        spkT_N4 = [60+(i*duration)]
        spkT_N5 = [80+(i*duration)]
        spkT_N6 = [100+(i*duration)]
        spkT_N7 = [120+(i*duration)]
        spkT_N8 = [140+(i*duration)]
        spkT_N9 = [160+(i*duration)]
        spkT_N10 = [180+(i*duration)]

        inpPop.spike_times=[spkT_N1,spkT_N2,spkT_N3,spkT_N4,spkT_N5,spkT_N6,spkT_N7,spkT_N8,spkT_N9,spkT_N10]

        simulate(duration)
    w = projTen_Ten.w

    img = np.ones((10,10))

    """
    Adapt the output like in Clopath et al. 2010.
    Depending on the connection, set another number to get another color.

    prepare a matrix of weights with different values for the different
    connections as mentioned in the Clopath et al., 2010 publication
    weak connections (< (2/3 of max. weight)) == 0
    strong unidirectional (> (2/3 of max. weight)) connections == 1.0
    strong bidirectional (> (2/3 of max. weight)) connections == 2.0
    """

    maxima = (np.nanmax(w)*2./3.)
    idx = np.where(w < maxima)
    img[idx[0],idx[1]] = 0.0

    # Strong biidirectional connections (> 2/3 of maximum weight) = 2.0
    idx_r = np.asarray(np.where(w >=maxima))
    for i in range(len(idx_r[0])):
        ix = (idx_r[0,i],idx_r[1,i])
        for j in range(len(idx_r[0])):
            ix2 = (idx_r[0,j],idx_r[1,j])
            if ix2 == (ix[1],ix[0]):
                img[ix[0],ix[1]] = 2.0
                img[ix[1],ix[0]] = 2.0
    # Set selfconnection weights to nan, because they not exist
    for i in range(10):
        w[i][i] = np.nan
        img[i,i]= np.nan

    # Start plotting
    plt.figure()
    plt.imshow(img.T,interpolation='none',cmap=myCmap(),vmin=0,vmax=2)
    plt.xlabel('Neuron Post',fontsize=20)
    plt.ylabel('Neuron Pre',fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig('Fig3_temporalCode.png',bbox_inches='tight')
    # plt.show()
    print("done")

if __name__ == "__main__":
    run()
