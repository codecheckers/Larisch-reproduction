/*
 *  ANNarchy-version: 4.6.8.1
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

    // Custom local parameter spike_times
    // std::vector< double > r ;
    std::vector< std::vector< long int > > spike_times ;
    std::vector< long int >  next_spike ;
    std::vector< int > idx_next_spike;
    long int _t;

    // Recompute the spike times
    void recompute_spike_times(){
        std::fill(next_spike.begin(), next_spike.end(), -10000);
        std::fill(idx_next_spike.begin(), idx_next_spike.end(), 0);
        for(int i=0; i< size; i++){
            if(!spike_times[i].empty()){
                int idx = 0;
                // Find the first spike time which is not in the past
                while(spike_times[i][idx] < _t){
                    idx++;
                }
                // Set the next spike
                if(idx < spike_times[i].size())
                    next_spike[i] = spike_times[i][idx];
                else
                    next_spike[i] = -10000;
            }
        }
    }


    // Access methods to the parameters and variables

    // Local variable r
    std::vector< double > get_r() { return r; }
    double get_single_r(int rk) { return r[rk]; }
    void set_r(std::vector< double > val) { r = val; }
    void set_single_r(int rk, double val) { r[rk] = val; }



    // Method called to initialize the data structures
    void init_population() {
        _active = true;

        // Local variable r
        r = std::vector<double>(size, 0.0);


        // Spiking variables
        spiked = std::vector<int>(0, 0);
        last_spike = std::vector<long int>(size, -10000L);



        // Mean Firing Rate
        _spike_history = std::vector< std::queue<long int> >(size, std::queue<long int>());
        _mean_fr_window = 0;
        _mean_fr_rate = 1.0;

        _t = 0;
        next_spike = std::vector<long int>(size, -10000);
        idx_next_spike = std::vector<int>(size, 0);
        this->recompute_spike_times();


    }

    // Method called to reset the population
    void reset() {

        spiked.clear();
        last_spike.clear();
        last_spike = std::vector<long int>(size, -10000L);



        _t = 0;
        this->recompute_spike_times();

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

        if(_active){
            spiked.clear();
            for(int i = 0; i < size; i++){
                // Emit spike
                if( _t == next_spike[i] ){
                    last_spike[i] = _t;
                    /*
                    while(++idx_next_spike[i]< spike_times[i].size()){
                        if(spike_times[i][idx_next_spike[i]] > _t)
                            break;
                    }
                    */
                    idx_next_spike[i]++ ;
                    if(idx_next_spike[i] < spike_times[i].size()){
                        next_spike[i] = spike_times[i][idx_next_spike[i]];
                    }
                    spiked.push_back(i);
                }
            }
            _t++;
        }

    }

    

    // Memory management: track the memory consumption
    long int size_in_bytes() {
        long int size_in_bytes = 0;
        // Parameters
        // Variables
        size_in_bytes += sizeof(double) * r.capacity();	// r
        
        return size_in_bytes;
    }

    // Memory management: track the memory consumption
    void clear() {
        // Variables
        r.clear();
        r.shrink_to_fit();
        
    }
};
