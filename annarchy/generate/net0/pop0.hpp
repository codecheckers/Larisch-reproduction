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
// Main Structure for the population of id 0 (pop0)
///////////////////////////////////////////////////////////////
struct PopStruct0{

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

    // Local parameter rates
    std::vector< double > rates;

    // Local variable p
    std::vector< double > p;

    // Local variable r
    std::vector< double > r;

    // Global operations

    // Random numbers
    std::vector<double> rand_0;
    std::uniform_real_distribution< double > dist_rand_0;



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

    // Local parameter rates
    std::vector< double > get_rates() { return rates; }
    double get_single_rates(int rk) { return rates[rk]; }
    void set_rates(std::vector< double > val) { rates = val; }
    void set_single_rates(int rk, double val) { rates[rk] = val; }

    // Local variable p
    std::vector< double > get_p() { return p; }
    double get_single_p(int rk) { return p[rk]; }
    void set_p(std::vector< double > val) { p = val; }
    void set_single_p(int rk, double val) { p[rk] = val; }

    // Local variable r
    std::vector< double > get_r() { return r; }
    double get_single_r(int rk) { return r[rk]; }
    void set_r(std::vector< double > val) { r = val; }
    void set_single_r(int rk, double val) { r[rk] = val; }



    // Method called to initialize the data structures
    void init_population() {
        _active = true;

        // Local parameter rates
        rates = std::vector<double>(size, 0.0);

        // Local variable p
        p = std::vector<double>(size, 0.0);

        // Local variable r
        r = std::vector<double>(size, 0.0);

        rand_0 = std::vector<double>(size, 0.0);


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

        dist_rand_0 = std::uniform_real_distribution< double >(0.0, 1.0);

    }

    // Method to draw new random numbers
    void update_rng() {

        if (_active){

            for(int i = 0; i < size; i++) {

                rand_0[i] = dist_rand_0(rng);

            }
        }

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

                // p = Uniform(0.0, 1.0) * 1000.0 / dt
                p[i] = 1000.0*rand_0[i]/dt;


            }
        } // active



        if( _active ) {
            for (int i = 0; i < size; i++) {
                // Spike emission
                if(p[i] < rates[i]){ // Condition is met
                    // Reset variables

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
        size_in_bytes += sizeof(double) * rates.capacity();	// rates
        // Variables
        size_in_bytes += sizeof(double) * p.capacity();	// p
        size_in_bytes += sizeof(double) * r.capacity();	// r

        return size_in_bytes;
    }

    // Memory management: track the memory consumption
    void clear() {
        // Variables
        p.clear();
        p.shrink_to_fit();
        r.clear();
        r.shrink_to_fit();

    }
};

