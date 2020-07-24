# cython: embedsignature=True
from cpython.exc cimport PyErr_CheckSignals
from libcpp.vector cimport vector
from libcpp.map cimport map, pair
from libcpp cimport bool
from math import ceil
import numpy as np
import sys
cimport numpy as np

import ANNarchy
from ANNarchy.core.cython_ext.Connector cimport LILConnectivity as LIL
from ANNarchy.core.cython_ext.Connector cimport CSRConnectivity, CSRConnectivityPre1st

cdef extern from "ANNarchy.h":

    # User-defined functions


    # User-defined constants


    # Data structures

    # Export Population 0 (pop0)
    cdef struct PopStruct0 :
        # Number of neurons
        int get_size()
        void set_size(int)
        # Maximum delay in steps
        int get_max_delay()
        void set_max_delay(int)
        void update_max_delay(int)
        # Activate/deactivate the population
        bool is_active()
        void set_active(bool)
        # Reset the population
        void reset()


        # Local parameter rates
        vector[double] get_rates()
        double get_single_rates(int rk)
        void set_rates(vector[double])
        void set_single_rates(int, double)

        # Local variable p
        vector[double] get_p()
        double get_single_p(int rk)
        void set_p(vector[double])
        void set_single_p(int, double)

        # Local variable r
        vector[double] get_r()
        double get_single_r(int rk)
        void set_r(vector[double])
        void set_single_r(int, double)




        # Compute firing rate
        void compute_firing_rate(double window)


        # memory management
        long int size_in_bytes()
        void clear()

    # Export Population 1 (pop1)
    cdef struct PopStruct1 :
        # Number of neurons
        int get_size()
        void set_size(int)
        # Maximum delay in steps
        int get_max_delay()
        void set_max_delay(int)
        void update_max_delay(int)
        # Activate/deactivate the population
        bool is_active()
        void set_active(bool)
        # Reset the population
        void reset()


        # Global parameter EL
        double  get_EL()
        void set_EL(double)

        # Global parameter VTrest
        double  get_VTrest()
        void set_VTrest(double)

        # Global parameter taux
        double  get_taux()
        void set_taux(double)

        # Local variable g_vm
        vector[double] get_g_vm()
        double get_single_g_vm(int rk)
        void set_g_vm(vector[double])
        void set_single_g_vm(int, double)

        # Local variable Spike
        vector[double] get_Spike()
        double get_single_Spike(int rk)
        void set_Spike(vector[double])
        void set_single_Spike(int, double)

        # Local variable Reset
        vector[double] get_Reset()
        double get_single_Reset(int rk)
        void set_Reset(vector[double])
        void set_single_Reset(int, double)

        # Local variable xtrace
        vector[double] get_xtrace()
        double get_single_xtrace(int rk)
        void set_xtrace(vector[double])
        void set_single_xtrace(int, double)

        # Local variable state
        vector[double] get_state()
        double get_single_state(int rk)
        void set_state(vector[double])
        void set_single_state(int, double)

        # Local variable r
        vector[double] get_r()
        double get_single_r(int rk)
        void set_r(vector[double])
        void set_single_r(int, double)




        # Compute firing rate
        void compute_firing_rate(double window)


        # memory management
        long int size_in_bytes()
        void clear()

    # Export Population 2 (pop2)
    cdef struct PopStruct2 :
        # Number of neurons
        int get_size()
        void set_size(int)
        # Maximum delay in steps
        int get_max_delay()
        void set_max_delay(int)
        void update_max_delay(int)
        # Activate/deactivate the population
        bool is_active()
        void set_active(bool)
        # Reset the population
        void reset()


        # Global parameter gL
        double  get_gL()
        void set_gL(double)

        # Global parameter DeltaT
        double  get_DeltaT()
        void set_DeltaT(double)

        # Global parameter tauw
        double  get_tauw()
        void set_tauw(double)

        # Global parameter a
        double  get_a()
        void set_a(double)

        # Global parameter b
        double  get_b()
        void set_b(double)

        # Global parameter EL
        double  get_EL()
        void set_EL(double)

        # Global parameter C
        double  get_C()
        void set_C(double)

        # Global parameter tauz
        double  get_tauz()
        void set_tauz(double)

        # Global parameter tauVT
        double  get_tauVT()
        void set_tauVT(double)

        # Global parameter Isp
        double  get_Isp()
        void set_Isp(double)

        # Global parameter VTMax
        double  get_VTMax()
        void set_VTMax(double)

        # Global parameter VTrest
        double  get_VTrest()
        void set_VTrest(double)

        # Global parameter taux
        double  get_taux()
        void set_taux(double)

        # Global parameter tauLTD
        double  get_tauLTD()
        void set_tauLTD(double)

        # Global parameter tauLTP
        double  get_tauLTP()
        void set_tauLTP(double)

        # Global parameter taumean
        double  get_taumean()
        void set_taumean(double)

        # Global parameter tau_gExc
        double  get_tau_gExc()
        void set_tau_gExc(double)

        # Local parameter inter_vm
        vector[double] get_inter_vm()
        double get_single_inter_vm(int rk)
        void set_inter_vm(vector[double])
        void set_single_inter_vm(int, double)

        # Local variable vm
        vector[double] get_vm()
        double get_single_vm(int rk)
        void set_vm(vector[double])
        void set_single_vm(int, double)

        # Local variable vmean
        vector[double] get_vmean()
        double get_single_vmean(int rk)
        void set_vmean(vector[double])
        void set_single_vmean(int, double)

        # Local variable umeanLTD
        vector[double] get_umeanLTD()
        double get_single_umeanLTD(int rk)
        void set_umeanLTD(vector[double])
        void set_single_umeanLTD(int, double)

        # Local variable umeanLTP
        vector[double] get_umeanLTP()
        double get_single_umeanLTP(int rk)
        void set_umeanLTP(vector[double])
        void set_single_umeanLTP(int, double)

        # Local variable xtrace
        vector[double] get_xtrace()
        double get_single_xtrace(int rk)
        void set_xtrace(vector[double])
        void set_single_xtrace(int, double)

        # Local variable wad
        vector[double] get_wad()
        double get_single_wad(int rk)
        void set_wad(vector[double])
        void set_single_wad(int, double)

        # Local variable z
        vector[double] get_z()
        double get_single_z(int rk)
        void set_z(vector[double])
        void set_single_z(int, double)

        # Local variable VT
        vector[double] get_VT()
        double get_single_VT(int rk)
        void set_VT(vector[double])
        void set_single_VT(int, double)

        # Local variable g_Exc
        vector[double] get_g_Exc()
        double get_single_g_Exc(int rk)
        void set_g_Exc(vector[double])
        void set_single_g_Exc(int, double)

        # Local variable state
        vector[double] get_state()
        double get_single_state(int rk)
        void set_state(vector[double])
        void set_single_state(int, double)

        # Local variable Spike
        vector[double] get_Spike()
        double get_single_Spike(int rk)
        void set_Spike(vector[double])
        void set_single_Spike(int, double)

        # Local variable resetvar
        vector[double] get_resetvar()
        double get_single_resetvar(int rk)
        void set_resetvar(vector[double])
        void set_single_resetvar(int, double)

        # Local variable r
        vector[double] get_r()
        double get_single_r(int rk)
        void set_r(vector[double])
        void set_single_r(int, double)




        # Compute firing rate
        void compute_firing_rate(double window)


        # memory management
        long int size_in_bytes()
        void clear()


    # Export Projection 0
    cdef struct ProjStruct0 :
        # Flags
        bool _transmission
        bool _plasticity
        bool _update
        int _update_period
        long _update_offset
        # Size
        int get_size()
        int nb_synapses(int)
        void set_size(int)


        # LIL Connectivity
        vector[int] get_post_rank()
        vector[vector[int]] get_pre_rank()
        void set_post_rank(vector[int])
        void set_pre_rank(vector[vector[int]])
        void inverse_connectivity_matrix()

        # Local variable w
        double w








        # memory management
        long int size_in_bytes()
        void clear()

    # Export Projection 1
    cdef struct ProjStruct1 :
        # Flags
        bool _transmission
        bool _plasticity
        bool _update
        int _update_period
        long _update_offset
        # Size
        int get_size()
        int nb_synapses(int)
        void set_size(int)


        # LIL Connectivity
        vector[int] get_post_rank()
        vector[vector[int]] get_pre_rank()
        void set_post_rank(vector[int])
        void set_pre_rank(vector[vector[int]])
        void inverse_connectivity_matrix()

        # Local variable w
        vector[vector[double]] get_w()
        vector[double] get_dendrite_w(int)
        double get_synapse_w(int, int)
        void set_w(vector[vector[double]])
        void set_dendrite_w(int, vector[double])
        void set_synapse_w(int, int, double)




        # Global parameter vmean_fix
        double get_vmean_fix()
        void set_vmean_fix(double)

        # Global parameter urefsquare
        double get_urefsquare()
        void set_urefsquare(double)

        # Global parameter thetaLTD
        double get_thetaLTD()
        void set_thetaLTD(double)

        # Global parameter thetaLTP
        double get_thetaLTP()
        void set_thetaLTP(double)

        # Global parameter aLTD
        double get_aLTD()
        void set_aLTD(double)

        # Global parameter aLTP
        double get_aLTP()
        void set_aLTP(double)

        # Global parameter wMin
        double get_wMin()
        void set_wMin(double)

        # Global parameter wMax
        double get_wMax()
        void set_wMax(double)

        # Global parameter transmit
        double get_transmit()
        void set_transmit(double)

        # Global parameter set_fix
        double get_set_fix()
        void set_set_fix(double)

        # Local variable ltdTerm_fix
        vector[vector[double]] get_ltdTerm_fix()
        vector[double] get_dendrite_ltdTerm_fix(int)
        double get_synapse_ltdTerm_fix(int, int)
        void set_ltdTerm_fix(vector[vector[double]])
        void set_dendrite_ltdTerm_fix(int, vector[double])
        void set_synapse_ltdTerm_fix(int, int, double)

        # Local variable ltdTerm
        vector[vector[double]] get_ltdTerm()
        vector[double] get_dendrite_ltdTerm(int)
        double get_synapse_ltdTerm(int, int)
        void set_ltdTerm(vector[vector[double]])
        void set_dendrite_ltdTerm(int, vector[double])
        void set_synapse_ltdTerm(int, int, double)

        # Local variable ltpTerm
        vector[vector[double]] get_ltpTerm()
        vector[double] get_dendrite_ltpTerm(int)
        double get_synapse_ltpTerm(int, int)
        void set_ltpTerm(vector[vector[double]])
        void set_dendrite_ltpTerm(int, vector[double])
        void set_synapse_ltpTerm(int, int, double)

        # Local variable deltaW
        vector[vector[double]] get_deltaW()
        vector[double] get_dendrite_deltaW(int)
        double get_synapse_deltaW(int, int)
        void set_deltaW(vector[vector[double]])
        void set_dendrite_deltaW(int, vector[double])
        void set_synapse_deltaW(int, int, double)





        # memory management
        long int size_in_bytes()
        void clear()


    # Monitors
    cdef cppclass Monitor:
        vector[int] ranks
        int period_
        int period_offset_
        long offset_

    void addRecorder(Monitor*)
    void removeRecorder(Monitor*)

    # Population 0 (pop0) : Monitor
    cdef cppclass PopRecorder0 (Monitor):
        PopRecorder0(vector[int], int, int, long) except +
        long int size_in_bytes()
        void clear()

        vector[vector[double]] rates
        bool record_rates

        vector[vector[double]] p
        bool record_p

        vector[vector[double]] r
        bool record_r

        map[int, vector[long]] spike
        bool record_spike
        void clear_spike()

    # Population 1 (pop1) : Monitor
    cdef cppclass PopRecorder1 (Monitor):
        PopRecorder1(vector[int], int, int, long) except +
        long int size_in_bytes()
        void clear()

        vector[double] EL
        bool record_EL

        vector[double] VTrest
        bool record_VTrest

        vector[double] taux
        bool record_taux

        vector[vector[double]] g_vm
        bool record_g_vm

        vector[vector[double]] Spike
        bool record_Spike

        vector[vector[double]] Reset
        bool record_Reset

        vector[vector[double]] xtrace
        bool record_xtrace

        vector[vector[double]] state
        bool record_state

        vector[vector[double]] r
        bool record_r

        map[int, vector[long]] spike
        bool record_spike
        void clear_spike()

    # Population 2 (pop2) : Monitor
    cdef cppclass PopRecorder2 (Monitor):
        PopRecorder2(vector[int], int, int, long) except +
        long int size_in_bytes()
        void clear()

        vector[double] gL
        bool record_gL

        vector[double] DeltaT
        bool record_DeltaT

        vector[double] tauw
        bool record_tauw

        vector[double] a
        bool record_a

        vector[double] b
        bool record_b

        vector[double] EL
        bool record_EL

        vector[double] C
        bool record_C

        vector[double] tauz
        bool record_tauz

        vector[double] tauVT
        bool record_tauVT

        vector[double] Isp
        bool record_Isp

        vector[double] VTMax
        bool record_VTMax

        vector[double] VTrest
        bool record_VTrest

        vector[double] taux
        bool record_taux

        vector[double] tauLTD
        bool record_tauLTD

        vector[double] tauLTP
        bool record_tauLTP

        vector[double] taumean
        bool record_taumean

        vector[double] tau_gExc
        bool record_tau_gExc

        vector[vector[double]] inter_vm
        bool record_inter_vm

        vector[vector[double]] vm
        bool record_vm

        vector[vector[double]] vmean
        bool record_vmean

        vector[vector[double]] umeanLTD
        bool record_umeanLTD

        vector[vector[double]] umeanLTP
        bool record_umeanLTP

        vector[vector[double]] xtrace
        bool record_xtrace

        vector[vector[double]] wad
        bool record_wad

        vector[vector[double]] z
        bool record_z

        vector[vector[double]] VT
        bool record_VT

        vector[vector[double]] g_Exc
        bool record_g_Exc

        vector[vector[double]] state
        bool record_state

        vector[vector[double]] Spike
        bool record_Spike

        vector[vector[double]] resetvar
        bool record_resetvar

        vector[vector[double]] r
        bool record_r

        map[int, vector[long]] spike
        bool record_spike
        void clear_spike()

    # Projection 0 : Monitor
    cdef cppclass ProjRecorder0 (Monitor):
        ProjRecorder0(vector[int], int, int, long) except +

        vector[double] w
        bool record_w

    # Projection 1 : Monitor
    cdef cppclass ProjRecorder1 (Monitor):
        ProjRecorder1(vector[int], int, int, long) except +

        vector[double] vmean_fix
        bool record_vmean_fix

        vector[double] urefsquare
        bool record_urefsquare

        vector[double] thetaLTD
        bool record_thetaLTD

        vector[double] thetaLTP
        bool record_thetaLTP

        vector[double] aLTD
        bool record_aLTD

        vector[double] aLTP
        bool record_aLTP

        vector[double] wMin
        bool record_wMin

        vector[double] wMax
        bool record_wMax

        vector[double] transmit
        bool record_transmit

        vector[double] set_fix
        bool record_set_fix

        vector[vector[vector[double]]] ltdTerm_fix
        bool record_ltdTerm_fix

        vector[vector[vector[double]]] ltdTerm
        bool record_ltdTerm

        vector[vector[vector[double]]] ltpTerm
        bool record_ltpTerm

        vector[vector[vector[double]]] deltaW
        bool record_deltaW

        vector[vector[vector[double]]] w
        bool record_w


    # Instances

    PopStruct0 pop0
    PopStruct1 pop1
    PopStruct2 pop2

    ProjStruct0 proj0
    ProjStruct1 proj1

    # Methods
    void initialize(double, long)
    void init_rng_dist()
    void setSeed(long)
    void run(int nbSteps) nogil
    int run_until(int steps, vector[int] populations, bool or_and)
    void step()

    # Time
    long getTime()
    void setTime(long)

    # dt
    double getDt()
    void setDt(double dt_)


    # Number of threads
    void setNumberThreads(int)


# Population wrappers

# Wrapper for population 0 (pop0)
cdef class pop0_wrapper :

    def __cinit__(self, size, max_delay):

        pop0.set_size(size)
        pop0.set_max_delay(max_delay)
    # Number of neurons
    property size:
        def __get__(self):
            return pop0.get_size()
    # Reset the population
    def reset(self):
        pop0.reset()
    # Set the maximum delay of outgoing projections
    def set_max_delay(self, val):
        pop0.set_max_delay(val)
    # Updates the maximum delay of outgoing projections and rebuilds the arrays
    def update_max_delay(self, val):
        pop0.update_max_delay(val)
    # Allows the population to compute
    def activate(self, bool val):
        pop0.set_active(val)


    # Local parameter rates
    cpdef np.ndarray get_rates(self):
        return np.array(pop0.get_rates())
    cpdef set_rates(self, np.ndarray value):
        pop0.set_rates( value )
    cpdef double get_single_rates(self, int rank):
        return pop0.get_single_rates(rank)
    cpdef set_single_rates(self, int rank, value):
        pop0.set_single_rates(rank, value)

    # Local variable p
    cpdef np.ndarray get_p(self):
        return np.array(pop0.get_p())
    cpdef set_p(self, np.ndarray value):
        pop0.set_p( value )
    cpdef double get_single_p(self, int rank):
        return pop0.get_single_p(rank)
    cpdef set_single_p(self, int rank, value):
        pop0.set_single_p(rank, value)

    # Local variable r
    cpdef np.ndarray get_r(self):
        return np.array(pop0.get_r())
    cpdef set_r(self, np.ndarray value):
        pop0.set_r( value )
    cpdef double get_single_r(self, int rank):
        return pop0.get_single_r(rank)
    cpdef set_single_r(self, int rank, value):
        pop0.set_single_r(rank, value)





    # Compute firing rate
    cpdef compute_firing_rate(self, double window):
        pop0.compute_firing_rate(window)


    # memory management
    def size_in_bytes(self):
        return pop0.size_in_bytes()

    def clear(self):
        return pop0.clear()

# Wrapper for population 1 (pop1)
cdef class pop1_wrapper :

    def __cinit__(self, size, max_delay):

        pop1.set_size(size)
        pop1.set_max_delay(max_delay)
    # Number of neurons
    property size:
        def __get__(self):
            return pop1.get_size()
    # Reset the population
    def reset(self):
        pop1.reset()
    # Set the maximum delay of outgoing projections
    def set_max_delay(self, val):
        pop1.set_max_delay(val)
    # Updates the maximum delay of outgoing projections and rebuilds the arrays
    def update_max_delay(self, val):
        pop1.update_max_delay(val)
    # Allows the population to compute
    def activate(self, bool val):
        pop1.set_active(val)


    # Global parameter EL
    cpdef double get_EL(self):
        return pop1.get_EL()
    cpdef set_EL(self, double value):
        pop1.set_EL(value)

    # Global parameter VTrest
    cpdef double get_VTrest(self):
        return pop1.get_VTrest()
    cpdef set_VTrest(self, double value):
        pop1.set_VTrest(value)

    # Global parameter taux
    cpdef double get_taux(self):
        return pop1.get_taux()
    cpdef set_taux(self, double value):
        pop1.set_taux(value)

    # Local variable g_vm
    cpdef np.ndarray get_g_vm(self):
        return np.array(pop1.get_g_vm())
    cpdef set_g_vm(self, np.ndarray value):
        pop1.set_g_vm( value )
    cpdef double get_single_g_vm(self, int rank):
        return pop1.get_single_g_vm(rank)
    cpdef set_single_g_vm(self, int rank, value):
        pop1.set_single_g_vm(rank, value)

    # Local variable Spike
    cpdef np.ndarray get_Spike(self):
        return np.array(pop1.get_Spike())
    cpdef set_Spike(self, np.ndarray value):
        pop1.set_Spike( value )
    cpdef double get_single_Spike(self, int rank):
        return pop1.get_single_Spike(rank)
    cpdef set_single_Spike(self, int rank, value):
        pop1.set_single_Spike(rank, value)

    # Local variable Reset
    cpdef np.ndarray get_Reset(self):
        return np.array(pop1.get_Reset())
    cpdef set_Reset(self, np.ndarray value):
        pop1.set_Reset( value )
    cpdef double get_single_Reset(self, int rank):
        return pop1.get_single_Reset(rank)
    cpdef set_single_Reset(self, int rank, value):
        pop1.set_single_Reset(rank, value)

    # Local variable xtrace
    cpdef np.ndarray get_xtrace(self):
        return np.array(pop1.get_xtrace())
    cpdef set_xtrace(self, np.ndarray value):
        pop1.set_xtrace( value )
    cpdef double get_single_xtrace(self, int rank):
        return pop1.get_single_xtrace(rank)
    cpdef set_single_xtrace(self, int rank, value):
        pop1.set_single_xtrace(rank, value)

    # Local variable state
    cpdef np.ndarray get_state(self):
        return np.array(pop1.get_state())
    cpdef set_state(self, np.ndarray value):
        pop1.set_state( value )
    cpdef double get_single_state(self, int rank):
        return pop1.get_single_state(rank)
    cpdef set_single_state(self, int rank, value):
        pop1.set_single_state(rank, value)

    # Local variable r
    cpdef np.ndarray get_r(self):
        return np.array(pop1.get_r())
    cpdef set_r(self, np.ndarray value):
        pop1.set_r( value )
    cpdef double get_single_r(self, int rank):
        return pop1.get_single_r(rank)
    cpdef set_single_r(self, int rank, value):
        pop1.set_single_r(rank, value)





    # Compute firing rate
    cpdef compute_firing_rate(self, double window):
        pop1.compute_firing_rate(window)


    # memory management
    def size_in_bytes(self):
        return pop1.size_in_bytes()

    def clear(self):
        return pop1.clear()

# Wrapper for population 2 (pop2)
cdef class pop2_wrapper :

    def __cinit__(self, size, max_delay):

        pop2.set_size(size)
        pop2.set_max_delay(max_delay)
    # Number of neurons
    property size:
        def __get__(self):
            return pop2.get_size()
    # Reset the population
    def reset(self):
        pop2.reset()
    # Set the maximum delay of outgoing projections
    def set_max_delay(self, val):
        pop2.set_max_delay(val)
    # Updates the maximum delay of outgoing projections and rebuilds the arrays
    def update_max_delay(self, val):
        pop2.update_max_delay(val)
    # Allows the population to compute
    def activate(self, bool val):
        pop2.set_active(val)


    # Global parameter gL
    cpdef double get_gL(self):
        return pop2.get_gL()
    cpdef set_gL(self, double value):
        pop2.set_gL(value)

    # Global parameter DeltaT
    cpdef double get_DeltaT(self):
        return pop2.get_DeltaT()
    cpdef set_DeltaT(self, double value):
        pop2.set_DeltaT(value)

    # Global parameter tauw
    cpdef double get_tauw(self):
        return pop2.get_tauw()
    cpdef set_tauw(self, double value):
        pop2.set_tauw(value)

    # Global parameter a
    cpdef double get_a(self):
        return pop2.get_a()
    cpdef set_a(self, double value):
        pop2.set_a(value)

    # Global parameter b
    cpdef double get_b(self):
        return pop2.get_b()
    cpdef set_b(self, double value):
        pop2.set_b(value)

    # Global parameter EL
    cpdef double get_EL(self):
        return pop2.get_EL()
    cpdef set_EL(self, double value):
        pop2.set_EL(value)

    # Global parameter C
    cpdef double get_C(self):
        return pop2.get_C()
    cpdef set_C(self, double value):
        pop2.set_C(value)

    # Global parameter tauz
    cpdef double get_tauz(self):
        return pop2.get_tauz()
    cpdef set_tauz(self, double value):
        pop2.set_tauz(value)

    # Global parameter tauVT
    cpdef double get_tauVT(self):
        return pop2.get_tauVT()
    cpdef set_tauVT(self, double value):
        pop2.set_tauVT(value)

    # Global parameter Isp
    cpdef double get_Isp(self):
        return pop2.get_Isp()
    cpdef set_Isp(self, double value):
        pop2.set_Isp(value)

    # Global parameter VTMax
    cpdef double get_VTMax(self):
        return pop2.get_VTMax()
    cpdef set_VTMax(self, double value):
        pop2.set_VTMax(value)

    # Global parameter VTrest
    cpdef double get_VTrest(self):
        return pop2.get_VTrest()
    cpdef set_VTrest(self, double value):
        pop2.set_VTrest(value)

    # Global parameter taux
    cpdef double get_taux(self):
        return pop2.get_taux()
    cpdef set_taux(self, double value):
        pop2.set_taux(value)

    # Global parameter tauLTD
    cpdef double get_tauLTD(self):
        return pop2.get_tauLTD()
    cpdef set_tauLTD(self, double value):
        pop2.set_tauLTD(value)

    # Global parameter tauLTP
    cpdef double get_tauLTP(self):
        return pop2.get_tauLTP()
    cpdef set_tauLTP(self, double value):
        pop2.set_tauLTP(value)

    # Global parameter taumean
    cpdef double get_taumean(self):
        return pop2.get_taumean()
    cpdef set_taumean(self, double value):
        pop2.set_taumean(value)

    # Global parameter tau_gExc
    cpdef double get_tau_gExc(self):
        return pop2.get_tau_gExc()
    cpdef set_tau_gExc(self, double value):
        pop2.set_tau_gExc(value)

    # Local parameter inter_vm
    cpdef np.ndarray get_inter_vm(self):
        return np.array(pop2.get_inter_vm())
    cpdef set_inter_vm(self, np.ndarray value):
        pop2.set_inter_vm( value )
    cpdef double get_single_inter_vm(self, int rank):
        return pop2.get_single_inter_vm(rank)
    cpdef set_single_inter_vm(self, int rank, value):
        pop2.set_single_inter_vm(rank, value)

    # Local variable vm
    cpdef np.ndarray get_vm(self):
        return np.array(pop2.get_vm())
    cpdef set_vm(self, np.ndarray value):
        pop2.set_vm( value )
    cpdef double get_single_vm(self, int rank):
        return pop2.get_single_vm(rank)
    cpdef set_single_vm(self, int rank, value):
        pop2.set_single_vm(rank, value)

    # Local variable vmean
    cpdef np.ndarray get_vmean(self):
        return np.array(pop2.get_vmean())
    cpdef set_vmean(self, np.ndarray value):
        pop2.set_vmean( value )
    cpdef double get_single_vmean(self, int rank):
        return pop2.get_single_vmean(rank)
    cpdef set_single_vmean(self, int rank, value):
        pop2.set_single_vmean(rank, value)

    # Local variable umeanLTD
    cpdef np.ndarray get_umeanLTD(self):
        return np.array(pop2.get_umeanLTD())
    cpdef set_umeanLTD(self, np.ndarray value):
        pop2.set_umeanLTD( value )
    cpdef double get_single_umeanLTD(self, int rank):
        return pop2.get_single_umeanLTD(rank)
    cpdef set_single_umeanLTD(self, int rank, value):
        pop2.set_single_umeanLTD(rank, value)

    # Local variable umeanLTP
    cpdef np.ndarray get_umeanLTP(self):
        return np.array(pop2.get_umeanLTP())
    cpdef set_umeanLTP(self, np.ndarray value):
        pop2.set_umeanLTP( value )
    cpdef double get_single_umeanLTP(self, int rank):
        return pop2.get_single_umeanLTP(rank)
    cpdef set_single_umeanLTP(self, int rank, value):
        pop2.set_single_umeanLTP(rank, value)

    # Local variable xtrace
    cpdef np.ndarray get_xtrace(self):
        return np.array(pop2.get_xtrace())
    cpdef set_xtrace(self, np.ndarray value):
        pop2.set_xtrace( value )
    cpdef double get_single_xtrace(self, int rank):
        return pop2.get_single_xtrace(rank)
    cpdef set_single_xtrace(self, int rank, value):
        pop2.set_single_xtrace(rank, value)

    # Local variable wad
    cpdef np.ndarray get_wad(self):
        return np.array(pop2.get_wad())
    cpdef set_wad(self, np.ndarray value):
        pop2.set_wad( value )
    cpdef double get_single_wad(self, int rank):
        return pop2.get_single_wad(rank)
    cpdef set_single_wad(self, int rank, value):
        pop2.set_single_wad(rank, value)

    # Local variable z
    cpdef np.ndarray get_z(self):
        return np.array(pop2.get_z())
    cpdef set_z(self, np.ndarray value):
        pop2.set_z( value )
    cpdef double get_single_z(self, int rank):
        return pop2.get_single_z(rank)
    cpdef set_single_z(self, int rank, value):
        pop2.set_single_z(rank, value)

    # Local variable VT
    cpdef np.ndarray get_VT(self):
        return np.array(pop2.get_VT())
    cpdef set_VT(self, np.ndarray value):
        pop2.set_VT( value )
    cpdef double get_single_VT(self, int rank):
        return pop2.get_single_VT(rank)
    cpdef set_single_VT(self, int rank, value):
        pop2.set_single_VT(rank, value)

    # Local variable g_Exc
    cpdef np.ndarray get_g_Exc(self):
        return np.array(pop2.get_g_Exc())
    cpdef set_g_Exc(self, np.ndarray value):
        pop2.set_g_Exc( value )
    cpdef double get_single_g_Exc(self, int rank):
        return pop2.get_single_g_Exc(rank)
    cpdef set_single_g_Exc(self, int rank, value):
        pop2.set_single_g_Exc(rank, value)

    # Local variable state
    cpdef np.ndarray get_state(self):
        return np.array(pop2.get_state())
    cpdef set_state(self, np.ndarray value):
        pop2.set_state( value )
    cpdef double get_single_state(self, int rank):
        return pop2.get_single_state(rank)
    cpdef set_single_state(self, int rank, value):
        pop2.set_single_state(rank, value)

    # Local variable Spike
    cpdef np.ndarray get_Spike(self):
        return np.array(pop2.get_Spike())
    cpdef set_Spike(self, np.ndarray value):
        pop2.set_Spike( value )
    cpdef double get_single_Spike(self, int rank):
        return pop2.get_single_Spike(rank)
    cpdef set_single_Spike(self, int rank, value):
        pop2.set_single_Spike(rank, value)

    # Local variable resetvar
    cpdef np.ndarray get_resetvar(self):
        return np.array(pop2.get_resetvar())
    cpdef set_resetvar(self, np.ndarray value):
        pop2.set_resetvar( value )
    cpdef double get_single_resetvar(self, int rank):
        return pop2.get_single_resetvar(rank)
    cpdef set_single_resetvar(self, int rank, value):
        pop2.set_single_resetvar(rank, value)

    # Local variable r
    cpdef np.ndarray get_r(self):
        return np.array(pop2.get_r())
    cpdef set_r(self, np.ndarray value):
        pop2.set_r( value )
    cpdef double get_single_r(self, int rank):
        return pop2.get_single_r(rank)
    cpdef set_single_r(self, int rank, value):
        pop2.set_single_r(rank, value)





    # Compute firing rate
    cpdef compute_firing_rate(self, double window):
        pop2.compute_firing_rate(window)


    # memory management
    def size_in_bytes(self):
        return pop2.size_in_bytes()

    def clear(self):
        return pop2.clear()


# Projection wrappers

# Wrapper for projection 0
cdef class proj0_wrapper :

    def __cinit__(self, synapses):

        cdef LIL syn = synapses
        cdef int size = syn.size
        cdef int nb_post = syn.post_rank.size()
        proj0.set_size( size )
        proj0.set_post_rank( syn.post_rank )
        proj0.set_pre_rank( syn.pre_rank )

        # Use only the first weight
        proj0.w = syn.w[0][0]




    property size:
        def __get__(self):
            return proj0.get_size()

    def nb_synapses(self, int n):
        return proj0.nb_synapses(n)

    # Transmission flag
    def _get_transmission(self):
        return proj0._transmission
    def _set_transmission(self, bool l):
        proj0._transmission = l

    # Update flag
    def _get_update(self):
        return proj0._update
    def _set_update(self, bool l):
        proj0._update = l

    # Plasticity flag
    def _get_plasticity(self):
        return proj0._plasticity
    def _set_plasticity(self, bool l):
        proj0._plasticity = l

    # Update period
    def _get_update_period(self):
        return proj0._update_period
    def _set_update_period(self, int l):
        proj0._update_period = l

    # Update offset
    def _get_update_offset(self):
        return proj0._update_offset
    def _set_update_offset(self, long l):
        proj0._update_offset = l


    # Connectivity
    def post_rank(self):
        return proj0.get_post_rank()
    def set_post_rank(self, val):
        proj0.set_post_rank(val)
        proj0.inverse_connectivity_matrix()
    def pre_rank(self, int n):
        return proj0.get_pre_rank()[n]
    def pre_rank_all(self):
        return proj0.get_pre_rank()
    def set_pre_rank(self, val):
        proj0.set_pre_rank(val)
        proj0.inverse_connectivity_matrix()

    # Local variable w
    def get_w(self):
        return proj0.w
    def set_w(self, value):
        proj0.w = value
    def get_dendrite_w(self, int rank):
        return proj0.w
    def set_dendrite_w(self, int rank, double value):
        proj0.w = value
    def get_synapse_w(self, int rank_post, int rank_pre):
        return proj0.w
    def set_synapse_w(self, int rank_post, int rank_pre, double value):
        proj0.w = value







    # memory management
    def size_in_bytes(self):
        return proj0.size_in_bytes()

    def clear(self):
        return proj0.clear()

# Wrapper for projection 1
cdef class proj1_wrapper :

    def __cinit__(self, synapses):

        cdef LIL syn = synapses
        cdef int size = syn.size
        cdef int nb_post = syn.post_rank.size()
        proj1.set_size( size )
        proj1.set_post_rank( syn.post_rank )
        proj1.set_pre_rank( syn.pre_rank )

        proj1.set_w(syn.w)




    property size:
        def __get__(self):
            return proj1.get_size()

    def nb_synapses(self, int n):
        return proj1.nb_synapses(n)

    # Transmission flag
    def _get_transmission(self):
        return proj1._transmission
    def _set_transmission(self, bool l):
        proj1._transmission = l

    # Update flag
    def _get_update(self):
        return proj1._update
    def _set_update(self, bool l):
        proj1._update = l

    # Plasticity flag
    def _get_plasticity(self):
        return proj1._plasticity
    def _set_plasticity(self, bool l):
        proj1._plasticity = l

    # Update period
    def _get_update_period(self):
        return proj1._update_period
    def _set_update_period(self, int l):
        proj1._update_period = l

    # Update offset
    def _get_update_offset(self):
        return proj1._update_offset
    def _set_update_offset(self, long l):
        proj1._update_offset = l


    # Connectivity
    def post_rank(self):
        return proj1.get_post_rank()
    def set_post_rank(self, val):
        proj1.set_post_rank(val)
        proj1.inverse_connectivity_matrix()
    def pre_rank(self, int n):
        return proj1.get_pre_rank()[n]
    def pre_rank_all(self):
        return proj1.get_pre_rank()
    def set_pre_rank(self, val):
        proj1.set_pre_rank(val)
        proj1.inverse_connectivity_matrix()

    # Local variable w
    def get_w(self):
        return proj1.get_w()
    def set_w(self, value):
        proj1.set_w( value )
    def get_dendrite_w(self, int rank):
        return proj1.get_dendrite_w(rank)
    def set_dendrite_w(self, int rank, vector[double] value):
        proj1.set_dendrite_w(rank, value)
    def get_synapse_w(self, int rank_post, int rank_pre):
        return proj1.get_synapse_w(rank_post, rank_pre)
    def set_synapse_w(self, int rank_post, int rank_pre, double value):
        proj1.set_synapse_w(rank_post, rank_pre, value)



    # Global parameter vmean_fix
    def get_vmean_fix(self):
        return proj1.get_vmean_fix()
    def set_vmean_fix(self, value):
        proj1.set_vmean_fix(value)

    # Global parameter urefsquare
    def get_urefsquare(self):
        return proj1.get_urefsquare()
    def set_urefsquare(self, value):
        proj1.set_urefsquare(value)

    # Global parameter thetaLTD
    def get_thetaLTD(self):
        return proj1.get_thetaLTD()
    def set_thetaLTD(self, value):
        proj1.set_thetaLTD(value)

    # Global parameter thetaLTP
    def get_thetaLTP(self):
        return proj1.get_thetaLTP()
    def set_thetaLTP(self, value):
        proj1.set_thetaLTP(value)

    # Global parameter aLTD
    def get_aLTD(self):
        return proj1.get_aLTD()
    def set_aLTD(self, value):
        proj1.set_aLTD(value)

    # Global parameter aLTP
    def get_aLTP(self):
        return proj1.get_aLTP()
    def set_aLTP(self, value):
        proj1.set_aLTP(value)

    # Global parameter wMin
    def get_wMin(self):
        return proj1.get_wMin()
    def set_wMin(self, value):
        proj1.set_wMin(value)

    # Global parameter wMax
    def get_wMax(self):
        return proj1.get_wMax()
    def set_wMax(self, value):
        proj1.set_wMax(value)

    # Global parameter transmit
    def get_transmit(self):
        return proj1.get_transmit()
    def set_transmit(self, value):
        proj1.set_transmit(value)

    # Global parameter set_fix
    def get_set_fix(self):
        return proj1.get_set_fix()
    def set_set_fix(self, value):
        proj1.set_set_fix(value)

    # Local variable ltdTerm_fix
    def get_ltdTerm_fix(self):
        return proj1.get_ltdTerm_fix()
    def set_ltdTerm_fix(self, value):
        proj1.set_ltdTerm_fix( value )
    def get_dendrite_ltdTerm_fix(self, int rank):
        return proj1.get_dendrite_ltdTerm_fix(rank)
    def set_dendrite_ltdTerm_fix(self, int rank, vector[double] value):
        proj1.set_dendrite_ltdTerm_fix(rank, value)
    def get_synapse_ltdTerm_fix(self, int rank_post, int rank_pre):
        return proj1.get_synapse_ltdTerm_fix(rank_post, rank_pre)
    def set_synapse_ltdTerm_fix(self, int rank_post, int rank_pre, double value):
        proj1.set_synapse_ltdTerm_fix(rank_post, rank_pre, value)

    # Local variable ltdTerm
    def get_ltdTerm(self):
        return proj1.get_ltdTerm()
    def set_ltdTerm(self, value):
        proj1.set_ltdTerm( value )
    def get_dendrite_ltdTerm(self, int rank):
        return proj1.get_dendrite_ltdTerm(rank)
    def set_dendrite_ltdTerm(self, int rank, vector[double] value):
        proj1.set_dendrite_ltdTerm(rank, value)
    def get_synapse_ltdTerm(self, int rank_post, int rank_pre):
        return proj1.get_synapse_ltdTerm(rank_post, rank_pre)
    def set_synapse_ltdTerm(self, int rank_post, int rank_pre, double value):
        proj1.set_synapse_ltdTerm(rank_post, rank_pre, value)

    # Local variable ltpTerm
    def get_ltpTerm(self):
        return proj1.get_ltpTerm()
    def set_ltpTerm(self, value):
        proj1.set_ltpTerm( value )
    def get_dendrite_ltpTerm(self, int rank):
        return proj1.get_dendrite_ltpTerm(rank)
    def set_dendrite_ltpTerm(self, int rank, vector[double] value):
        proj1.set_dendrite_ltpTerm(rank, value)
    def get_synapse_ltpTerm(self, int rank_post, int rank_pre):
        return proj1.get_synapse_ltpTerm(rank_post, rank_pre)
    def set_synapse_ltpTerm(self, int rank_post, int rank_pre, double value):
        proj1.set_synapse_ltpTerm(rank_post, rank_pre, value)

    # Local variable deltaW
    def get_deltaW(self):
        return proj1.get_deltaW()
    def set_deltaW(self, value):
        proj1.set_deltaW( value )
    def get_dendrite_deltaW(self, int rank):
        return proj1.get_dendrite_deltaW(rank)
    def set_dendrite_deltaW(self, int rank, vector[double] value):
        proj1.set_dendrite_deltaW(rank, value)
    def get_synapse_deltaW(self, int rank_post, int rank_pre):
        return proj1.get_synapse_deltaW(rank_post, rank_pre)
    def set_synapse_deltaW(self, int rank_post, int rank_pre, double value):
        proj1.set_synapse_deltaW(rank_post, rank_pre, value)





    # memory management
    def size_in_bytes(self):
        return proj1.size_in_bytes()

    def clear(self):
        return proj1.clear()


# Monitor wrappers
cdef class Monitor_wrapper:
    cdef Monitor *thisptr
    def __cinit__(self, list ranks, int period, int period_offset, long offset):
        pass
    property ranks:
        def __get__(self): return self.thisptr.ranks
        def __set__(self, val): self.thisptr.ranks = val
    property period:
        def __get__(self): return self.thisptr.period_
        def __set__(self, val): self.thisptr.period_ = val
    property offset:
        def __get__(self): return self.thisptr.offset_
        def __set__(self, val): self.thisptr.offset_ = val
    property period_offset:
        def __get__(self): return self.thisptr.period_offset_
        def __set__(self, val): self.thisptr.period_offset_ = val

def add_recorder(Monitor_wrapper recorder):
    addRecorder(recorder.thisptr)
def remove_recorder(Monitor_wrapper recorder):
    removeRecorder(recorder.thisptr)


# Population Monitor wrapper
cdef class PopRecorder0_wrapper(Monitor_wrapper):
    def __cinit__(self, list ranks, int period, period_offset, long offset):
        self.thisptr = new PopRecorder0(ranks, period, period_offset, offset)

    def size_in_bytes(self):
        return (<PopRecorder0 *>self.thisptr).size_in_bytes()

    def clear(self):
        (<PopRecorder0 *>self.thisptr).clear()


    property rates:
        def __get__(self): return (<PopRecorder0 *>self.thisptr).rates
        def __set__(self, val): (<PopRecorder0 *>self.thisptr).rates = val
    property record_rates:
        def __get__(self): return (<PopRecorder0 *>self.thisptr).record_rates
        def __set__(self, val): (<PopRecorder0 *>self.thisptr).record_rates = val
    def clear_rates(self):
        (<PopRecorder0 *>self.thisptr).rates.clear()

    property p:
        def __get__(self): return (<PopRecorder0 *>self.thisptr).p
        def __set__(self, val): (<PopRecorder0 *>self.thisptr).p = val
    property record_p:
        def __get__(self): return (<PopRecorder0 *>self.thisptr).record_p
        def __set__(self, val): (<PopRecorder0 *>self.thisptr).record_p = val
    def clear_p(self):
        (<PopRecorder0 *>self.thisptr).p.clear()

    property r:
        def __get__(self): return (<PopRecorder0 *>self.thisptr).r
        def __set__(self, val): (<PopRecorder0 *>self.thisptr).r = val
    property record_r:
        def __get__(self): return (<PopRecorder0 *>self.thisptr).record_r
        def __set__(self, val): (<PopRecorder0 *>self.thisptr).record_r = val
    def clear_r(self):
        (<PopRecorder0 *>self.thisptr).r.clear()

    property spike:
        def __get__(self): return (<PopRecorder0 *>self.thisptr).spike
        def __set__(self, val): (<PopRecorder0 *>self.thisptr).spike = val
    property record_spike:
        def __get__(self): return (<PopRecorder0 *>self.thisptr).record_spike
        def __set__(self, val): (<PopRecorder0 *>self.thisptr).record_spike = val
    def clear_spike(self):
        (<PopRecorder0 *>self.thisptr).clear_spike()

# Population Monitor wrapper
cdef class PopRecorder1_wrapper(Monitor_wrapper):
    def __cinit__(self, list ranks, int period, period_offset, long offset):
        self.thisptr = new PopRecorder1(ranks, period, period_offset, offset)

    def size_in_bytes(self):
        return (<PopRecorder1 *>self.thisptr).size_in_bytes()

    def clear(self):
        (<PopRecorder1 *>self.thisptr).clear()


    property EL:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).EL
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).EL = val
    property record_EL:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).record_EL
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).record_EL = val
    def clear_EL(self):
        (<PopRecorder1 *>self.thisptr).EL.clear()

    property VTrest:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).VTrest
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).VTrest = val
    property record_VTrest:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).record_VTrest
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).record_VTrest = val
    def clear_VTrest(self):
        (<PopRecorder1 *>self.thisptr).VTrest.clear()

    property taux:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).taux
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).taux = val
    property record_taux:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).record_taux
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).record_taux = val
    def clear_taux(self):
        (<PopRecorder1 *>self.thisptr).taux.clear()

    property g_vm:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).g_vm
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).g_vm = val
    property record_g_vm:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).record_g_vm
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).record_g_vm = val
    def clear_g_vm(self):
        (<PopRecorder1 *>self.thisptr).g_vm.clear()

    property Spike:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).Spike
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).Spike = val
    property record_Spike:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).record_Spike
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).record_Spike = val
    def clear_Spike(self):
        (<PopRecorder1 *>self.thisptr).Spike.clear()

    property Reset:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).Reset
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).Reset = val
    property record_Reset:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).record_Reset
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).record_Reset = val
    def clear_Reset(self):
        (<PopRecorder1 *>self.thisptr).Reset.clear()

    property xtrace:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).xtrace
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).xtrace = val
    property record_xtrace:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).record_xtrace
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).record_xtrace = val
    def clear_xtrace(self):
        (<PopRecorder1 *>self.thisptr).xtrace.clear()

    property state:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).state
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).state = val
    property record_state:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).record_state
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).record_state = val
    def clear_state(self):
        (<PopRecorder1 *>self.thisptr).state.clear()

    property r:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).r
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).r = val
    property record_r:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).record_r
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).record_r = val
    def clear_r(self):
        (<PopRecorder1 *>self.thisptr).r.clear()

    property spike:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).spike
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).spike = val
    property record_spike:
        def __get__(self): return (<PopRecorder1 *>self.thisptr).record_spike
        def __set__(self, val): (<PopRecorder1 *>self.thisptr).record_spike = val
    def clear_spike(self):
        (<PopRecorder1 *>self.thisptr).clear_spike()

# Population Monitor wrapper
cdef class PopRecorder2_wrapper(Monitor_wrapper):
    def __cinit__(self, list ranks, int period, period_offset, long offset):
        self.thisptr = new PopRecorder2(ranks, period, period_offset, offset)

    def size_in_bytes(self):
        return (<PopRecorder2 *>self.thisptr).size_in_bytes()

    def clear(self):
        (<PopRecorder2 *>self.thisptr).clear()


    property gL:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).gL
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).gL = val
    property record_gL:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_gL
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_gL = val
    def clear_gL(self):
        (<PopRecorder2 *>self.thisptr).gL.clear()

    property DeltaT:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).DeltaT
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).DeltaT = val
    property record_DeltaT:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_DeltaT
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_DeltaT = val
    def clear_DeltaT(self):
        (<PopRecorder2 *>self.thisptr).DeltaT.clear()

    property tauw:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).tauw
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).tauw = val
    property record_tauw:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_tauw
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_tauw = val
    def clear_tauw(self):
        (<PopRecorder2 *>self.thisptr).tauw.clear()

    property a:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).a
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).a = val
    property record_a:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_a
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_a = val
    def clear_a(self):
        (<PopRecorder2 *>self.thisptr).a.clear()

    property b:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).b
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).b = val
    property record_b:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_b
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_b = val
    def clear_b(self):
        (<PopRecorder2 *>self.thisptr).b.clear()

    property EL:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).EL
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).EL = val
    property record_EL:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_EL
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_EL = val
    def clear_EL(self):
        (<PopRecorder2 *>self.thisptr).EL.clear()

    property C:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).C
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).C = val
    property record_C:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_C
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_C = val
    def clear_C(self):
        (<PopRecorder2 *>self.thisptr).C.clear()

    property tauz:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).tauz
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).tauz = val
    property record_tauz:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_tauz
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_tauz = val
    def clear_tauz(self):
        (<PopRecorder2 *>self.thisptr).tauz.clear()

    property tauVT:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).tauVT
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).tauVT = val
    property record_tauVT:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_tauVT
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_tauVT = val
    def clear_tauVT(self):
        (<PopRecorder2 *>self.thisptr).tauVT.clear()

    property Isp:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).Isp
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).Isp = val
    property record_Isp:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_Isp
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_Isp = val
    def clear_Isp(self):
        (<PopRecorder2 *>self.thisptr).Isp.clear()

    property VTMax:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).VTMax
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).VTMax = val
    property record_VTMax:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_VTMax
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_VTMax = val
    def clear_VTMax(self):
        (<PopRecorder2 *>self.thisptr).VTMax.clear()

    property VTrest:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).VTrest
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).VTrest = val
    property record_VTrest:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_VTrest
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_VTrest = val
    def clear_VTrest(self):
        (<PopRecorder2 *>self.thisptr).VTrest.clear()

    property taux:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).taux
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).taux = val
    property record_taux:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_taux
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_taux = val
    def clear_taux(self):
        (<PopRecorder2 *>self.thisptr).taux.clear()

    property tauLTD:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).tauLTD
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).tauLTD = val
    property record_tauLTD:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_tauLTD
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_tauLTD = val
    def clear_tauLTD(self):
        (<PopRecorder2 *>self.thisptr).tauLTD.clear()

    property tauLTP:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).tauLTP
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).tauLTP = val
    property record_tauLTP:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_tauLTP
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_tauLTP = val
    def clear_tauLTP(self):
        (<PopRecorder2 *>self.thisptr).tauLTP.clear()

    property taumean:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).taumean
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).taumean = val
    property record_taumean:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_taumean
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_taumean = val
    def clear_taumean(self):
        (<PopRecorder2 *>self.thisptr).taumean.clear()

    property tau_gExc:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).tau_gExc
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).tau_gExc = val
    property record_tau_gExc:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_tau_gExc
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_tau_gExc = val
    def clear_tau_gExc(self):
        (<PopRecorder2 *>self.thisptr).tau_gExc.clear()

    property inter_vm:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).inter_vm
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).inter_vm = val
    property record_inter_vm:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_inter_vm
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_inter_vm = val
    def clear_inter_vm(self):
        (<PopRecorder2 *>self.thisptr).inter_vm.clear()

    property vm:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).vm
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).vm = val
    property record_vm:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_vm
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_vm = val
    def clear_vm(self):
        (<PopRecorder2 *>self.thisptr).vm.clear()

    property vmean:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).vmean
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).vmean = val
    property record_vmean:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_vmean
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_vmean = val
    def clear_vmean(self):
        (<PopRecorder2 *>self.thisptr).vmean.clear()

    property umeanLTD:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).umeanLTD
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).umeanLTD = val
    property record_umeanLTD:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_umeanLTD
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_umeanLTD = val
    def clear_umeanLTD(self):
        (<PopRecorder2 *>self.thisptr).umeanLTD.clear()

    property umeanLTP:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).umeanLTP
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).umeanLTP = val
    property record_umeanLTP:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_umeanLTP
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_umeanLTP = val
    def clear_umeanLTP(self):
        (<PopRecorder2 *>self.thisptr).umeanLTP.clear()

    property xtrace:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).xtrace
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).xtrace = val
    property record_xtrace:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_xtrace
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_xtrace = val
    def clear_xtrace(self):
        (<PopRecorder2 *>self.thisptr).xtrace.clear()

    property wad:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).wad
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).wad = val
    property record_wad:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_wad
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_wad = val
    def clear_wad(self):
        (<PopRecorder2 *>self.thisptr).wad.clear()

    property z:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).z
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).z = val
    property record_z:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_z
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_z = val
    def clear_z(self):
        (<PopRecorder2 *>self.thisptr).z.clear()

    property VT:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).VT
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).VT = val
    property record_VT:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_VT
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_VT = val
    def clear_VT(self):
        (<PopRecorder2 *>self.thisptr).VT.clear()

    property g_Exc:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).g_Exc
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).g_Exc = val
    property record_g_Exc:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_g_Exc
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_g_Exc = val
    def clear_g_Exc(self):
        (<PopRecorder2 *>self.thisptr).g_Exc.clear()

    property state:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).state
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).state = val
    property record_state:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_state
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_state = val
    def clear_state(self):
        (<PopRecorder2 *>self.thisptr).state.clear()

    property Spike:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).Spike
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).Spike = val
    property record_Spike:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_Spike
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_Spike = val
    def clear_Spike(self):
        (<PopRecorder2 *>self.thisptr).Spike.clear()

    property resetvar:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).resetvar
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).resetvar = val
    property record_resetvar:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_resetvar
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_resetvar = val
    def clear_resetvar(self):
        (<PopRecorder2 *>self.thisptr).resetvar.clear()

    property r:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).r
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).r = val
    property record_r:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_r
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_r = val
    def clear_r(self):
        (<PopRecorder2 *>self.thisptr).r.clear()

    property spike:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).spike
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).spike = val
    property record_spike:
        def __get__(self): return (<PopRecorder2 *>self.thisptr).record_spike
        def __set__(self, val): (<PopRecorder2 *>self.thisptr).record_spike = val
    def clear_spike(self):
        (<PopRecorder2 *>self.thisptr).clear_spike()

# Projection Monitor wrapper
cdef class ProjRecorder0_wrapper(Monitor_wrapper):
    def __cinit__(self, list ranks, int period, int period_offset, long offset):
        self.thisptr = new ProjRecorder0(ranks, period, period_offset, offset)

    property w:
        def __get__(self): return (<ProjRecorder0 *>self.thisptr).w
        def __set__(self, val): (<ProjRecorder0 *>self.thisptr).w = val
    property record_w:
        def __get__(self): return (<ProjRecorder0 *>self.thisptr).record_w
        def __set__(self, val): (<ProjRecorder0 *>self.thisptr).record_w = val
    def clear_w(self):
        (<ProjRecorder0 *>self.thisptr).w.clear()

# Projection Monitor wrapper
cdef class ProjRecorder1_wrapper(Monitor_wrapper):
    def __cinit__(self, list ranks, int period, int period_offset, long offset):
        self.thisptr = new ProjRecorder1(ranks, period, period_offset, offset)

    property vmean_fix:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).vmean_fix
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).vmean_fix = val
    property record_vmean_fix:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_vmean_fix
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_vmean_fix = val
    def clear_vmean_fix(self):
        (<ProjRecorder1 *>self.thisptr).vmean_fix.clear()

    property urefsquare:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).urefsquare
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).urefsquare = val
    property record_urefsquare:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_urefsquare
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_urefsquare = val
    def clear_urefsquare(self):
        (<ProjRecorder1 *>self.thisptr).urefsquare.clear()

    property thetaLTD:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).thetaLTD
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).thetaLTD = val
    property record_thetaLTD:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_thetaLTD
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_thetaLTD = val
    def clear_thetaLTD(self):
        (<ProjRecorder1 *>self.thisptr).thetaLTD.clear()

    property thetaLTP:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).thetaLTP
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).thetaLTP = val
    property record_thetaLTP:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_thetaLTP
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_thetaLTP = val
    def clear_thetaLTP(self):
        (<ProjRecorder1 *>self.thisptr).thetaLTP.clear()

    property aLTD:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).aLTD
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).aLTD = val
    property record_aLTD:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_aLTD
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_aLTD = val
    def clear_aLTD(self):
        (<ProjRecorder1 *>self.thisptr).aLTD.clear()

    property aLTP:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).aLTP
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).aLTP = val
    property record_aLTP:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_aLTP
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_aLTP = val
    def clear_aLTP(self):
        (<ProjRecorder1 *>self.thisptr).aLTP.clear()

    property wMin:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).wMin
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).wMin = val
    property record_wMin:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_wMin
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_wMin = val
    def clear_wMin(self):
        (<ProjRecorder1 *>self.thisptr).wMin.clear()

    property wMax:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).wMax
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).wMax = val
    property record_wMax:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_wMax
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_wMax = val
    def clear_wMax(self):
        (<ProjRecorder1 *>self.thisptr).wMax.clear()

    property transmit:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).transmit
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).transmit = val
    property record_transmit:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_transmit
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_transmit = val
    def clear_transmit(self):
        (<ProjRecorder1 *>self.thisptr).transmit.clear()

    property set_fix:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).set_fix
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).set_fix = val
    property record_set_fix:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_set_fix
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_set_fix = val
    def clear_set_fix(self):
        (<ProjRecorder1 *>self.thisptr).set_fix.clear()

    property ltdTerm_fix:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).ltdTerm_fix
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).ltdTerm_fix = val
    property record_ltdTerm_fix:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_ltdTerm_fix
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_ltdTerm_fix = val
    def clear_ltdTerm_fix(self):
        (<ProjRecorder1 *>self.thisptr).ltdTerm_fix.clear()

    property ltdTerm:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).ltdTerm
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).ltdTerm = val
    property record_ltdTerm:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_ltdTerm
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_ltdTerm = val
    def clear_ltdTerm(self):
        (<ProjRecorder1 *>self.thisptr).ltdTerm.clear()

    property ltpTerm:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).ltpTerm
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).ltpTerm = val
    property record_ltpTerm:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_ltpTerm
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_ltpTerm = val
    def clear_ltpTerm(self):
        (<ProjRecorder1 *>self.thisptr).ltpTerm.clear()

    property deltaW:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).deltaW
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).deltaW = val
    property record_deltaW:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_deltaW
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_deltaW = val
    def clear_deltaW(self):
        (<ProjRecorder1 *>self.thisptr).deltaW.clear()

    property w:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).w
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).w = val
    property record_w:
        def __get__(self): return (<ProjRecorder1 *>self.thisptr).record_w
        def __set__(self, val): (<ProjRecorder1 *>self.thisptr).record_w = val
    def clear_w(self):
        (<ProjRecorder1 *>self.thisptr).w.clear()


# User-defined functions


# User-defined constants


# Initialize the network
def pyx_create(double dt, long seed):
    initialize(dt, seed)

def pyx_init_rng_dist():
    init_rng_dist()

# Simple progressbar on the command line
def progress(count, total, status=''):
    """
    Prints a progress bar on the command line.

    adapted from: https://gist.github.com/vladignatyev/06860ec2040cb497f0f3

    Modification: The original code set the '\r' at the end, so the bar disappears when finished.
    I moved it to the front, so the last status remains.
    """
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('\r[%s] %s%s ...%s' % (bar, percents, '%', status))
    sys.stdout.flush()

# Simulation for the given number of steps
def pyx_run(int nb_steps, progress_bar):
    cdef int nb, rest
    cdef int batch = 1000
    if nb_steps < batch:
        with nogil:
            run(nb_steps)
    else:
        nb = int(nb_steps/batch)
        rest = nb_steps % batch
        for i in range(nb):
            with nogil:
                run(batch)
            PyErr_CheckSignals()
            if nb > 1 and progress_bar:
                progress(i+1, nb, 'simulate()')
        if rest > 0:
            run(rest)

        if (progress_bar):
            print('\n')

# Simulation for the given number of steps except if a criterion is reached
def pyx_run_until(int nb_steps, list populations, bool mode):
    cdef int nb
    nb = run_until(nb_steps, populations, mode)
    return nb

# Simulate for one step
def pyx_step():
    step()

# Access time
def set_time(t):
    setTime(t)
def get_time():
    return getTime()

# Access dt
def set_dt(double dt):
    setDt(dt)
def get_dt():
    return getDt()


# Set number of threads
def set_number_threads(int n):
    setNumberThreads(n)


# Set seed
def set_seed(long seed):
    setSeed(seed)
