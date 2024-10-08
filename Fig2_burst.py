"""
Python script to reproduce the burst spiking experiments (Fig.3 in Clopath et al. (2010)).

1. Experiment -- Weight change depending on the number of postsynaptic spikes:
 One presynaptic spike is paired with one,two or three postsynaptic spikes.
 The first postsynaptic spike arrives +/-10 ms after/befor the presynaptic spike.
 As mentioned in the original publication, the postsynaptic neuron have a firing rate of 50 Hz
 (20 ms between every postsynaptic spike ). The time delay refers to the time point of the presynaptic spike
 and the first postsynaptic spike.

2. Experiment -- Weight change depending on postsynaptic spiking frequency:
 One presynaptic spikes follows three postsynaptic spikes.
 The frequency for the postsynaptic spikes changes from 20 to 100 Hz.
 Between the presynaptic spike and the first postsynaptic spike is a time delay of +/- 10 ms.

3. Experiment -- Weight change depending on time lag between pre- and postsynaptic spike:
 Three postsynaptic spikes (with a frequency of 50 Hz) are paired with one presynaptic spike.
 The time lag between the occurrence of the first of the three postsynaptic spike to the occurrence of the one
 presynaptic spikes varies from -100 to 60 ms.
 The time between the postsynaptic spikes are 20 ms.

See Fig. 3 in Clopath et al. (2010)
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from ANNarchy import *
setup(dt=1)
from network import *

# Global parameters
duration = 240 #ms
d_t = 10 # time between pre- and first postsynaptic spike
t_1 = 110 # time point of the presynaptic spike

# Init the weights for the three tasks
initW1 = 0.0025
initW2 = 0.0025
initW3 = 0.01

# Time points
spike_times1 = [t_1]
spike_times2 = [t_1]

# Populations
"""
To control the spike timings of the AdEx neurons, two additional input populations
are used. The spike timing of the SpikeSourceArray can be determined with a
list of time points. """
# SpikeSourceArray populations to determine the spike timings of the  pre-
# and post neuron
inpPop1 = SpikeSourceArray(spike_times=spike_times1)
inpPop2 = SpikeSourceArray(spike_times=spike_times2)

# Populations for the neurons
popN1 = Population(geometry=1,neuron=AdExNeuron, name="N1")
popN2 = Population(geometry=1,neuron=AdExNeuron, name="N2")

# Projections
"""
Define simple projection from the input SpikeSourceArray populations
to the neuron populations.
If the neuron in the input population spikes,
1 ms later a spike in the connected AdEx neuron population is triggered.
"""

# Simple projection to propagate the spike of the input to the  presynaptic neuron
projST1_V1 = Projection(
    pre=inpPop1,
    post=popN1,
    target='Exc'
).connect_one_to_one(weights = 30.0)

# Simple projection to propagate the spike of the input to the postsynaptic neuron
projST2_V1 = Projection(
    pre=inpPop2,
    post=popN2,
    target='Exc'
).connect_one_to_one(weights = 30.0)

# Connection between the two neuron population to observe the weight changes
projN1_N2 = Projection(
    pre=popN1,
    post=popN2,
    target='Exc',
    synapse=ffSyn
).connect_one_to_one(weights = 0.01)

projN1_N2.set_fix = 1.0 # use a fix apmlitude for the LTD term

def run():
    # compile command to create the ANNarchy network
    compile()

    # Monitor objects
    dendrite = projN1_N2.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW'])

    """
    First experiment:
    Increase the number of postsynaptic spikes from 1 to 3.
    """
    dWSpk_pos = np.zeros(3)
    deltaWSpk_pos = np.zeros(3)

    # Loop to record weight change for one to three postsynaptic spikes after
    # the presynaptic one
    for i in range(3):
        projN1_N2.w = initW1 # Set the initial weight
        # Set the spike times for the second input neuron to
        # Determine the spike times of the postsynaptic neuron
        inpPop2.spike_times = np.linspace(t_1+d_t,t_1+d_t+20*(i),i+1).tolist()
        simulate(duration)
        d_w = m_d.get('w')
        dWSpk_pos[i] = d_w[-1]
        delta_w = m_d.get('deltaW')
        deltaWSpk_pos[i] = np.sum(delta_w)
        # Reset all variables
        reset()

    dWSpk_neg = np.zeros(3)
    deltaWSpk_neg = np.zeros(3)

    # Loop to record weight change for one to three postsynaptic spikes before the presynaptic one
    for i in range(3):
        projN1_N2.w = initW1
        inpPop2.spike_times = np.linspace(t_1-d_t,t_1-(d_t+20*(i)),i+1).tolist()
        simulate(duration)
        d_w = m_d.get('w')
        dWSpk_neg[i] = d_w[-1]
        delta_w = m_d.get('deltaW')
        deltaWSpk_neg[i] = np.sum(delta_w)
        reset()

    """
    Second experiment:
    Increase the postsynaptic firing from 20 to 100.
    """
    n_freq = 10
    rates = np.linspace(20,100,n_freq)
    dWBurst_pos = np.zeros(len(rates))
    deltaWBurst_pos = np.zeros(len(rates))

    # Loop to record weight change by changing the postsynaptic firing rate and
    # the first postsynaptic spike 10 ms after the presynaptic spike
    for i in range(n_freq):
        d_t2 = (1000./rates[i])
        projN1_N2.w = initW2
        inpPop2.spike_times = np.linspace(t_1+d_t,t_1+d_t+(d_t2*2),3).tolist()
        simulate(duration)
        d_w = m_d.get('w')
        dWBurst_pos[i] = d_w[-1]
        delta_w = m_d.get('deltaW')
        deltaWBurst_pos[i] = np.sum(delta_w)
        reset()


    dWBurst_neg = np.zeros(len(rates))
    deltaWBurst_neg = np.zeros(len(rates))
    # Loop to record weight change by changing the postsynaptic firing rate and
    # the first postsynaptic spike 10 ms before the presynaptic spike
    for i in range(n_freq):
        d_t2 = (1000./rates[i])
        projN1_N2.w = initW2
        inpPop2.spike_times = np.linspace(t_1-d_t,t_1-(d_t+d_t2*2),3).tolist()
        simulate(duration)
        d_w = m_d.get('w')
        dWBurst_neg[i] = d_w[-1]
        delta_w = m_d.get('deltaW')
        deltaWBurst_neg[i] = np.sum(delta_w)
        reset()

    """
    Third experiment:
    Change the delay between pre- and postsynaptic spikes from -100 to 60 ms.
    """
    lags = np.linspace(-100,60,33)
    dWLag_pos = np.zeros(len(lags))
    deltaWLag_pos = np.zeros(len(lags))
    # Loop over different delays between the presynaptic and the first postsynaptic spike
    for i in range(len(lags)):
        d_t2 = 20
        projN1_N2.w = initW3
        inpPop2.spike_times = np.linspace(t_1+lags[i],t_1+lags[i]+(d_t2*2),3).tolist()
        simulate(duration)
        d_w = m_d.get('w')
        dWLag_pos[i] = d_w[-1]
        delta_w = m_d.get('deltaW')
        deltaWLag_pos[i] = np.sum(delta_w)
        reset()

    # Start plotting
    fig = plt.figure(figsize=(12,10))
    gs=GridSpec(6,4)
    ax0= plt.subplot(gs[0:4,0:2])
    ax0.plot(np.clip(deltaWSpk_pos/initW1*100,0,250),'x',color='black',lw=5,ms=20)
    ax0.plot( (2*initW1 +deltaWSpk_neg)/initW1*100,'x',color='black',lw=5,ms=20)
    ax0.hlines(100,-0.2,3,colors='k')
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.xaxis.set_ticks_position('bottom')
    ax0.yaxis.set_ticks_position('left')
    plt.xticks(np.linspace(0,2,3),np.linspace(1,3,3),fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(-0.2,2.2)
    plt.xlabel('Number of Spikes',fontsize=25)
    plt.ylabel('Normalized weight (%)',fontsize=25)
    plt.ylim(0.0,250)

    ax1= plt.subplot(gs[0:4,2:4])
    ax1.plot(np.clip(deltaWBurst_pos/initW2*100,0,250),'--',color='black',lw=5,ms=15)
    ax1.plot((2*initW2+deltaWBurst_neg)/initW2*100,'--',color='black',lw=5,ms=15)
    ax1.hlines(100,0,n_freq,colors='k')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    plt.xticks(np.linspace(0,n_freq-1,5),np.linspace(20,100,5),fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylim(0.0,250)
    plt.xlabel('Frequency (Hz)',fontsize=25)

    ax2 = plt.subplot(gs[5:6,:])
    ax2.plot((initW3+deltaWLag_pos)/initW3*100,'-',color='black',lw=5)
    ax2.hlines(100,0,33,colors='k')
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    plt.ylim(0.0,300)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('Time lag (ms)',fontsize=25)
    fig.savefig('Fig2_burst.png',bbox_inches='tight')
    # plt.show()
    print("done")

if __name__ == "__main__":
    run()
