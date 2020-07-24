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
// Main Structure for the population of id 3 (N2)
///////////////////////////////////////////////////////////////
struct PopStruct3{

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

    // Global parameter gL
    double  gL ;

    // Global parameter DeltaT
    double  DeltaT ;

    // Global parameter tauw
    double  tauw ;

    // Global parameter a
    double  a ;

    // Global parameter b
    double  b ;

    // Global parameter EL
    double  EL ;

    // Global parameter C
    double  C ;

    // Global parameter tauz
    double  tauz ;

    // Global parameter tauVT
    double  tauVT ;

    // Global parameter Isp
    double  Isp ;

    // Global parameter VTMax
    double  VTMax ;

    // Global parameter VTrest
    double  VTrest ;

    // Global parameter taux
    double  taux ;

    // Global parameter tauLTD
    double  tauLTD ;

    // Global parameter tauLTP
    double  tauLTP ;

    // Global parameter taumean
    double  taumean ;

    // Global parameter tau_gExc
    double  tau_gExc ;

    // Local parameter inter_vm
    std::vector< double > inter_vm;

    // Local variable vm
    std::vector< double > vm;

    // Local variable vmean
    std::vector< double > vmean;

    // Local variable umeanLTD
    std::vector< double > umeanLTD;

    // Local variable umeanLTP
    std::vector< double > umeanLTP;

    // Local variable xtrace
    std::vector< double > xtrace;

    // Local variable wad
    std::vector< double > wad;

    // Local variable z
    std::vector< double > z;

    // Local variable VT
    std::vector< double > VT;

    // Local variable g_Exc
    std::vector< double > g_Exc;

    // Local variable state
    std::vector< double > state;

    // Local variable Spike
    std::vector< double > Spike;

    // Local variable resetvar
    std::vector< double > resetvar;

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

    // Global parameter gL
    double get_gL() { return gL; }
    void set_gL(double val) { gL = val; }

    // Global parameter DeltaT
    double get_DeltaT() { return DeltaT; }
    void set_DeltaT(double val) { DeltaT = val; }

    // Global parameter tauw
    double get_tauw() { return tauw; }
    void set_tauw(double val) { tauw = val; }

    // Global parameter a
    double get_a() { return a; }
    void set_a(double val) { a = val; }

    // Global parameter b
    double get_b() { return b; }
    void set_b(double val) { b = val; }

    // Global parameter EL
    double get_EL() { return EL; }
    void set_EL(double val) { EL = val; }

    // Global parameter C
    double get_C() { return C; }
    void set_C(double val) { C = val; }

    // Global parameter tauz
    double get_tauz() { return tauz; }
    void set_tauz(double val) { tauz = val; }

    // Global parameter tauVT
    double get_tauVT() { return tauVT; }
    void set_tauVT(double val) { tauVT = val; }

    // Global parameter Isp
    double get_Isp() { return Isp; }
    void set_Isp(double val) { Isp = val; }

    // Global parameter VTMax
    double get_VTMax() { return VTMax; }
    void set_VTMax(double val) { VTMax = val; }

    // Global parameter VTrest
    double get_VTrest() { return VTrest; }
    void set_VTrest(double val) { VTrest = val; }

    // Global parameter taux
    double get_taux() { return taux; }
    void set_taux(double val) { taux = val; }

    // Global parameter tauLTD
    double get_tauLTD() { return tauLTD; }
    void set_tauLTD(double val) { tauLTD = val; }

    // Global parameter tauLTP
    double get_tauLTP() { return tauLTP; }
    void set_tauLTP(double val) { tauLTP = val; }

    // Global parameter taumean
    double get_taumean() { return taumean; }
    void set_taumean(double val) { taumean = val; }

    // Global parameter tau_gExc
    double get_tau_gExc() { return tau_gExc; }
    void set_tau_gExc(double val) { tau_gExc = val; }

    // Local parameter inter_vm
    std::vector< double > get_inter_vm() { return inter_vm; }
    double get_single_inter_vm(int rk) { return inter_vm[rk]; }
    void set_inter_vm(std::vector< double > val) { inter_vm = val; }
    void set_single_inter_vm(int rk, double val) { inter_vm[rk] = val; }

    // Local variable vm
    std::vector< double > get_vm() { return vm; }
    double get_single_vm(int rk) { return vm[rk]; }
    void set_vm(std::vector< double > val) { vm = val; }
    void set_single_vm(int rk, double val) { vm[rk] = val; }

    // Local variable vmean
    std::vector< double > get_vmean() { return vmean; }
    double get_single_vmean(int rk) { return vmean[rk]; }
    void set_vmean(std::vector< double > val) { vmean = val; }
    void set_single_vmean(int rk, double val) { vmean[rk] = val; }

    // Local variable umeanLTD
    std::vector< double > get_umeanLTD() { return umeanLTD; }
    double get_single_umeanLTD(int rk) { return umeanLTD[rk]; }
    void set_umeanLTD(std::vector< double > val) { umeanLTD = val; }
    void set_single_umeanLTD(int rk, double val) { umeanLTD[rk] = val; }

    // Local variable umeanLTP
    std::vector< double > get_umeanLTP() { return umeanLTP; }
    double get_single_umeanLTP(int rk) { return umeanLTP[rk]; }
    void set_umeanLTP(std::vector< double > val) { umeanLTP = val; }
    void set_single_umeanLTP(int rk, double val) { umeanLTP[rk] = val; }

    // Local variable xtrace
    std::vector< double > get_xtrace() { return xtrace; }
    double get_single_xtrace(int rk) { return xtrace[rk]; }
    void set_xtrace(std::vector< double > val) { xtrace = val; }
    void set_single_xtrace(int rk, double val) { xtrace[rk] = val; }

    // Local variable wad
    std::vector< double > get_wad() { return wad; }
    double get_single_wad(int rk) { return wad[rk]; }
    void set_wad(std::vector< double > val) { wad = val; }
    void set_single_wad(int rk, double val) { wad[rk] = val; }

    // Local variable z
    std::vector< double > get_z() { return z; }
    double get_single_z(int rk) { return z[rk]; }
    void set_z(std::vector< double > val) { z = val; }
    void set_single_z(int rk, double val) { z[rk] = val; }

    // Local variable VT
    std::vector< double > get_VT() { return VT; }
    double get_single_VT(int rk) { return VT[rk]; }
    void set_VT(std::vector< double > val) { VT = val; }
    void set_single_VT(int rk, double val) { VT[rk] = val; }

    // Local variable g_Exc
    std::vector< double > get_g_Exc() { return g_Exc; }
    double get_single_g_Exc(int rk) { return g_Exc[rk]; }
    void set_g_Exc(std::vector< double > val) { g_Exc = val; }
    void set_single_g_Exc(int rk, double val) { g_Exc[rk] = val; }

    // Local variable state
    std::vector< double > get_state() { return state; }
    double get_single_state(int rk) { return state[rk]; }
    void set_state(std::vector< double > val) { state = val; }
    void set_single_state(int rk, double val) { state[rk] = val; }

    // Local variable Spike
    std::vector< double > get_Spike() { return Spike; }
    double get_single_Spike(int rk) { return Spike[rk]; }
    void set_Spike(std::vector< double > val) { Spike = val; }
    void set_single_Spike(int rk, double val) { Spike[rk] = val; }

    // Local variable resetvar
    std::vector< double > get_resetvar() { return resetvar; }
    double get_single_resetvar(int rk) { return resetvar[rk]; }
    void set_resetvar(std::vector< double > val) { resetvar = val; }
    void set_single_resetvar(int rk, double val) { resetvar[rk] = val; }

    // Local variable r
    std::vector< double > get_r() { return r; }
    double get_single_r(int rk) { return r[rk]; }
    void set_r(std::vector< double > val) { r = val; }
    void set_single_r(int rk, double val) { r[rk] = val; }



    // Method called to initialize the data structures
    void init_population() {
        _active = true;

        // Global parameter gL
        gL = 0.0;

        // Global parameter DeltaT
        DeltaT = 0.0;

        // Global parameter tauw
        tauw = 0.0;

        // Global parameter a
        a = 0.0;

        // Global parameter b
        b = 0.0;

        // Global parameter EL
        EL = 0.0;

        // Global parameter C
        C = 0.0;

        // Global parameter tauz
        tauz = 0.0;

        // Global parameter tauVT
        tauVT = 0.0;

        // Global parameter Isp
        Isp = 0.0;

        // Global parameter VTMax
        VTMax = 0.0;

        // Global parameter VTrest
        VTrest = 0.0;

        // Global parameter taux
        taux = 0.0;

        // Global parameter tauLTD
        tauLTD = 0.0;

        // Global parameter tauLTP
        tauLTP = 0.0;

        // Global parameter taumean
        taumean = 0.0;

        // Global parameter tau_gExc
        tau_gExc = 0.0;

        // Local parameter inter_vm
        inter_vm = std::vector<double>(size, 0.0);

        // Local variable vm
        vm = std::vector<double>(size, 0.0);

        // Local variable vmean
        vmean = std::vector<double>(size, 0.0);

        // Local variable umeanLTD
        umeanLTD = std::vector<double>(size, 0.0);

        // Local variable umeanLTP
        umeanLTP = std::vector<double>(size, 0.0);

        // Local variable xtrace
        xtrace = std::vector<double>(size, 0.0);

        // Local variable wad
        wad = std::vector<double>(size, 0.0);

        // Local variable z
        z = std::vector<double>(size, 0.0);

        // Local variable VT
        VT = std::vector<double>(size, 0.0);

        // Local variable g_Exc
        g_Exc = std::vector<double>(size, 0.0);

        // Local variable state
        state = std::vector<double>(size, 0.0);

        // Local variable Spike
        Spike = std::vector<double>(size, 0.0);

        // Local variable resetvar
        resetvar = std::vector<double>(size, 0.0);

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

                // dvm/dt = if state>=2:+3.462 else: if state==1:-vm + inter_vm +1/C*( - wad+b + z)+g_Exc else:1/C * ( -gL * (vm - EL) + gL * DeltaT * exp((vm - VT) / DeltaT) - wad + z ) + g_Exc
                double _vm = (state[i] >= 2 ? 3.4620000000000002 : (state[i] == 1 ? g_Exc[i] + inter_vm[i] - vm[i] + (b - wad[i] + z[i])/C : g_Exc[i] + (DeltaT*gL*exp((-VT[i] + vm[i])/DeltaT) + (-gL)*(-EL + vm[i]) - wad[i] + z[i])/C));

                // dvmean/dt = (pos(vm - EL)**2 - vmean)/taumean
                double _vmean = (-vmean[i] + pow(positive(-EL + vm[i]), 2))/taumean;

                // dumeanLTD/dt = (vm - umeanLTD)/tauLTD
                double _umeanLTD = (-umeanLTD[i] + vm[i])/tauLTD;

                // dumeanLTP/dt = (vm - umeanLTP)/tauLTP
                double _umeanLTP = (-umeanLTP[i] + vm[i])/tauLTP;

                // dxtrace /dt = (- xtrace )/taux
                double _xtrace = -xtrace[i]/taux;

                // dwad/dt = if state==1:+b else: (a * (vm - EL) - wad)/tauw
                double _wad = (state[i] == 1 ? b : (a*(-EL + vm[i]) - wad[i])/tauw);

                // dz/dt = if state==1:-z+Isp-10 else:-z/tauz
                double _z = (state[i] == 1 ? Isp - z[i] - 10 : (-z[i])/tauz);

                // dVT/dt =if state==2: 0 else: if state==1: -VT +VTMax else:(VTrest - VT)/tauVT
                double _VT = (state[i] == 2 ? 0 : (state[i] == 1 ? -VT[i] + VTMax : (-VT[i] + VTrest)/tauVT));

                // dg_Exc/dt = -g_Exc/tau_gExc
                double _g_Exc = -g_Exc[i]/tau_gExc;

                // dvm/dt = if state>=2:+3.462 else: if state==1:-vm + inter_vm +1/C*( - wad+b + z)+g_Exc else:1/C * ( -gL * (vm - EL) + gL * DeltaT * exp((vm - VT) / DeltaT) - wad + z ) + g_Exc
                vm[i] += dt*_vm ;


                // dvmean/dt = (pos(vm - EL)**2 - vmean)/taumean
                vmean[i] += dt*_vmean ;


                // dumeanLTD/dt = (vm - umeanLTD)/tauLTD
                umeanLTD[i] += dt*_umeanLTD ;


                // dumeanLTP/dt = (vm - umeanLTP)/tauLTP
                umeanLTP[i] += dt*_umeanLTP ;


                // dxtrace /dt = (- xtrace )/taux
                xtrace[i] += dt*_xtrace ;


                // dwad/dt = if state==1:+b else: (a * (vm - EL) - wad)/tauw
                wad[i] += dt*_wad ;


                // dz/dt = if state==1:-z+Isp-10 else:-z/tauz
                z[i] += dt*_z ;


                // dVT/dt =if state==2: 0 else: if state==1: -VT +VTMax else:(VTrest - VT)/tauVT
                VT[i] += dt*_VT ;


                // dg_Exc/dt = -g_Exc/tau_gExc
                g_Exc[i] += dt*_g_Exc ;


                // state = if state > 0: state-1 else:0
                state[i] = (state[i] > 0 ? state[i] - 1 : 0);


                // Spike = 0.0
                Spike[i] = 0.0;


                // dresetvar / dt = 1/(1.0) * (-resetvar)
                double _resetvar = -1.0*resetvar[i];

                // dresetvar / dt = 1/(1.0) * (-resetvar)
                resetvar[i] += dt*_resetvar ;


            }
        } // active



        if( _active ) {
            for (int i = 0; i < size; i++) {
                // Spike emission
                if(state[i] == 0 && vm[i] > VT[i]){ // Condition is met
                    // Reset variables

                    vm[i] = 29.399999999999999;

                    state[i] = 2.0;

                    Spike[i] = 1.0;

                    xtrace[i] += 1.0/taux;

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
        size_in_bytes += sizeof(double);	// gL
        size_in_bytes += sizeof(double);	// DeltaT
        size_in_bytes += sizeof(double);	// tauw
        size_in_bytes += sizeof(double);	// a
        size_in_bytes += sizeof(double);	// b
        size_in_bytes += sizeof(double);	// EL
        size_in_bytes += sizeof(double);	// C
        size_in_bytes += sizeof(double);	// tauz
        size_in_bytes += sizeof(double);	// tauVT
        size_in_bytes += sizeof(double);	// Isp
        size_in_bytes += sizeof(double);	// VTMax
        size_in_bytes += sizeof(double);	// VTrest
        size_in_bytes += sizeof(double);	// taux
        size_in_bytes += sizeof(double);	// tauLTD
        size_in_bytes += sizeof(double);	// tauLTP
        size_in_bytes += sizeof(double);	// taumean
        size_in_bytes += sizeof(double);	// tau_gExc
        size_in_bytes += sizeof(double) * inter_vm.capacity();	// inter_vm
        // Variables
        size_in_bytes += sizeof(double) * vm.capacity();	// vm
        size_in_bytes += sizeof(double) * vmean.capacity();	// vmean
        size_in_bytes += sizeof(double) * umeanLTD.capacity();	// umeanLTD
        size_in_bytes += sizeof(double) * umeanLTP.capacity();	// umeanLTP
        size_in_bytes += sizeof(double) * xtrace.capacity();	// xtrace
        size_in_bytes += sizeof(double) * wad.capacity();	// wad
        size_in_bytes += sizeof(double) * z.capacity();	// z
        size_in_bytes += sizeof(double) * VT.capacity();	// VT
        size_in_bytes += sizeof(double) * g_Exc.capacity();	// g_Exc
        size_in_bytes += sizeof(double) * state.capacity();	// state
        size_in_bytes += sizeof(double) * Spike.capacity();	// Spike
        size_in_bytes += sizeof(double) * resetvar.capacity();	// resetvar
        size_in_bytes += sizeof(double) * r.capacity();	// r

        return size_in_bytes;
    }

    // Memory management: track the memory consumption
    void clear() {
        // Variables
        vm.clear();
        vm.shrink_to_fit();
        vmean.clear();
        vmean.shrink_to_fit();
        umeanLTD.clear();
        umeanLTD.shrink_to_fit();
        umeanLTP.clear();
        umeanLTP.shrink_to_fit();
        xtrace.clear();
        xtrace.shrink_to_fit();
        wad.clear();
        wad.shrink_to_fit();
        z.clear();
        z.shrink_to_fit();
        VT.clear();
        VT.shrink_to_fit();
        g_Exc.clear();
        g_Exc.shrink_to_fit();
        state.clear();
        state.shrink_to_fit();
        Spike.clear();
        Spike.shrink_to_fit();
        resetvar.clear();
        resetvar.shrink_to_fit();
        r.clear();
        r.shrink_to_fit();

    }
};

