/*
 *  ANNarchy-version: 4.6.9.1
 */
#pragma once
#include "ANNarchy.h"
#include <random>


extern double dt;
extern long int t;
extern std::mt19937 rng;


///////////////////////////////////////////////////////////////
// Main Structure for the population of id 1 (pop1)
///////////////////////////////////////////////////////////////
struct PopStruct1{

    int size; // Number of neurons
    bool _active; // Allows to shut down the whole population
    int max_delay; // Maximum number of steps to store for delayed synaptic transmission

    // Access functions used by cython wrapper
    int get_size() { return size; }
    void set_size(int s) { size  = s; }
    int get_max_delay() { return max_delay; }
    void set_max_delay(int d) { max_delay  = d; }
    bool is_active() { return _active; }
    void set_active(bool val) { _active = val; }



    // Structures for managing spikes
    std::vector<long int> last_spike;
    std::vector<int> spiked;

    // Neuron specific parameters and variables

    // Global parameter EL
    double  EL ;

    // Global parameter VTrest
    double  VTrest ;

    // Global parameter taux
    double  taux ;

    // Local variable g_vm
    std::vector< double > g_vm;

    // Local variable Spike
    std::vector< double > Spike;

    // Local variable Reset
    std::vector< double > Reset;

    // Local variable xtrace
    std::vector< double > xtrace;

    // Local variable state
    std::vector< double > state;

    // Local variable r
    std::vector< double > r;

    // Global operations

    // Random numbers



    // Mean Firing rate
    std::vector< std::queue<long int> > _spike_history;
    long int _mean_fr_window;
    double _mean_fr_rate;
    void compute_firing_rate(double window){
        if(window>0.0){
            _mean_fr_window = int(window/dt);
            _mean_fr_rate = 1000./window;
        }
    };


    // Access methods to the parameters and variables

    // Global parameter EL
    double get_EL() { return EL; }
    void set_EL(double val) { EL = val; }

    // Global parameter VTrest
    double get_VTrest() { return VTrest; }
    void set_VTrest(double val) { VTrest = val; }

    // Global parameter taux
    double get_taux() { return taux; }
    void set_taux(double val) { taux = val; }

    // Local variable g_vm
    std::vector< double > get_g_vm() { return g_vm; }
    double get_single_g_vm(int rk) { return g_vm[rk]; }
    void set_g_vm(std::vector< double > val) { g_vm = val; }
    void set_single_g_vm(int rk, double val) { g_vm[rk] = val; }

    // Local variable Spike
    std::vector< double > get_Spike() { return Spike; }
    double get_single_Spike(int rk) { return Spike[rk]; }
    void set_Spike(std::vector< double > val) { Spike = val; }
    void set_single_Spike(int rk, double val) { Spike[rk] = val; }

    // Local variable Reset
    std::vector< double > get_Reset() { return Reset; }
    double get_single_Reset(int rk) { return Reset[rk]; }
    void set_Reset(std::vector< double > val) { Reset = val; }
    void set_single_Reset(int rk, double val) { Reset[rk] = val; }

    // Local variable xtrace
    std::vector< double > get_xtrace() { return xtrace; }
    double get_single_xtrace(int rk) { return xtrace[rk]; }
    void set_xtrace(std::vector< double > val) { xtrace = val; }
    void set_single_xtrace(int rk, double val) { xtrace[rk] = val; }

    // Local variable state
    std::vector< double > get_state() { return state; }
    double get_single_state(int rk) { return state[rk]; }
    void set_state(std::vector< double > val) { state = val; }
    void set_single_state(int rk, double val) { state[rk] = val; }

    // Local variable r
    std::vector< double > get_r() { return r; }
    double get_single_r(int rk) { return r[rk]; }
    void set_r(std::vector< double > val) { r = val; }
    void set_single_r(int rk, double val) { r[rk] = val; }



    // Method called to initialize the data structures
    void init_population() {
        _active = true;

        // Global parameter EL
        EL = 0.0;

        // Global parameter VTrest
        VTrest = 0.0;

        // Global parameter taux
        taux = 0.0;

        // Local variable g_vm
        g_vm = std::vector<double>(size, 0.0);

        // Local variable Spike
        Spike = std::vector<double>(size, 0.0);

        // Local variable Reset
        Reset = std::vector<double>(size, 0.0);

        // Local variable xtrace
        xtrace = std::vector<double>(size, 0.0);

        // Local variable state
        state = std::vector<double>(size, 0.0);

        // Local variable r
        r = std::vector<double>(size, 0.0);


        // Spiking variables
        spiked = std::vector<int>(0, 0);
        last_spike = std::vector<long int>(size, -10000L);



        // Mean Firing Rate
        _spike_history = std::vector< std::queue<long int> >(size, std::queue<long int>());
        _mean_fr_window = 0;
        _mean_fr_rate = 1.0;


    }

    // Method called to reset the population
    void reset() {

        spiked.clear();
        last_spike.clear();
        last_spike = std::vector<long int>(size, -10000L);



    }

    // Init rng dist
    void init_rng_dist() {

    }

    // Method to draw new random numbers
    void update_rng() {

    }

    // Method to update global operations on the population (min/max/mean...)
    void update_global_ops() {

    }

    // Method to enqueue output variables in case outgoing projections have non-zero delay
    void update_delay() {

    }

    // Method to dynamically change the size of the queue for delayed variables
    void update_max_delay(int value) {

    }

    // Main method to update neural variables
    void update() {

        if( _active ) {
            spiked.clear();

            // Updating local variables
            #pragma omp simd
            for(int i = 0; i < size; i++){

                // dg_vm/dt = EL/1000
                double _g_vm = (1.0/1000.0)*EL;

                // dg_vm/dt = EL/1000
                g_vm[i] += dt*_g_vm ;
                if(g_vm[i] < EL)
                    g_vm[i] = EL;


                // Spike = if state == 1: 1.0 else: 0.0
                Spike[i] = (state[i] == 1 ? 1.0 : 0.0);


                // dReset/dt = if state == 1: +1 else: -Reset
                double _Reset = (state[i] == 1 ? 1 : -Reset[i]);

                // dxtrace/dt = if state == 1: +1/taux else: -xtrace/taux
                double _xtrace = (state[i] == 1 ? 1.0/taux : (-xtrace[i])/taux);

                // dReset/dt = if state == 1: +1 else: -Reset
                Reset[i] += dt*_Reset ;


                // dxtrace/dt = if state == 1: +1/taux else: -xtrace/taux
                xtrace[i] += dt*_xtrace ;


                // state = if state >0: -1 else: 0
                state[i] = (state[i] > 0 ? -1 : 0);


            }
        } // active



        if( _active ) {
            for (int i = 0; i < size; i++) {
                // Spike emission
                if(g_vm[i] > VTrest){ // Condition is met
                    // Reset variables

                    g_vm[i] = EL;

                    state[i] = 1;

                    // Store the spike

                    {
                    spiked.push_back(i);
                    }
                    last_spike[i] = t;

                    // Refractory period


                    // Update the mean firing rate
                    if(_mean_fr_window> 0)
                        _spike_history[i].push(t);

                }

                // Update the mean firing rate
                if(_mean_fr_window> 0){
                    while((_spike_history[i].size() != 0)&&(_spike_history[i].front() <= t - _mean_fr_window)){
                        _spike_history[i].pop(); // Suppress spikes outside the window
                    }
                    r[i] = _mean_fr_rate * double(_spike_history[i].size());
                }



            }
        } // active

    }



    // Memory management: track the memory consumption
    long int size_in_bytes() {
        long int size_in_bytes = 0;
        // Parameters
        size_in_bytes += sizeof(double);	// EL
        size_in_bytes += sizeof(double);	// VTrest
        size_in_bytes += sizeof(double);	// taux
        // Variables
        size_in_bytes += sizeof(double) * g_vm.capacity();	// g_vm
        size_in_bytes += sizeof(double) * Spike.capacity();	// Spike
        size_in_bytes += sizeof(double) * Reset.capacity();	// Reset
        size_in_bytes += sizeof(double) * xtrace.capacity();	// xtrace
        size_in_bytes += sizeof(double) * state.capacity();	// state
        size_in_bytes += sizeof(double) * r.capacity();	// r

        return size_in_bytes;
    }

    // Memory management: track the memory consumption
    void clear() {
        // Variables
        g_vm.clear();
        g_vm.shrink_to_fit();
        Spike.clear();
        Spike.shrink_to_fit();
        Reset.clear();
        Reset.shrink_to_fit();
        xtrace.clear();
        xtrace.shrink_to_fit();
        state.clear();
        state.shrink_to_fit();
        r.clear();
        r.shrink_to_fit();

    }
};

