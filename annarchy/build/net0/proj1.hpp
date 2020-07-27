#pragma once

#include "pop1.hpp"
#include "pop3.hpp"



extern PopStruct1 pop1;
extern PopStruct3 pop3;


/////////////////////////////////////////////////////////////////////////////
// proj1: pop1 -> N2 with target Exc
/////////////////////////////////////////////////////////////////////////////
struct ProjStruct1{
    // Number of dendrites
    int size;

    // Transmission and plasticity flags
    bool _transmission, _plasticity, _update;
    int _update_period;
    long int _update_offset;


    // Connectivity
    std::vector<int> post_rank;
    std::vector< std::vector< int > > pre_rank;

    // Single weight in the projection
    double w;


    std::map< int, std::vector< std::pair<int, int> > > inv_pre_rank ;
    std::vector< int > inv_post_rank ;








    // Method called to initialize the projection
    void init_projection() {
        _transmission = true;
        _update = true;
        _plasticity = true;
        _update_period = 1;
        _update_offset = 0L;




        // Inverse the connectivity matrix if spiking neurons
        inverse_connectivity_matrix();







    }

    // Spiking networks: inverse the connectivity matrix
    void inverse_connectivity_matrix() {

        inv_pre_rank =  std::map< int, std::vector< std::pair<int, int> > > ();
        for(int i=0; i<pre_rank.size(); i++){
            for(int j=0; j<pre_rank[i].size(); j++){
                inv_pre_rank[pre_rank[i][j]].push_back(std::pair<int, int>(i,j));
            }
        }
        inv_post_rank =  std::vector< int > (pop3.size, -1);
        for(int i=0; i<post_rank.size(); i++){
            inv_post_rank[post_rank[i]] = i;
        }

    }

    // Spiking networks: update maximum delay when non-uniform
    void update_max_delay(int d){

    }

    // Computes the weighted sum of inputs or updates the conductances
    void compute_psp() {

        int nb_post;
        double sum;
        
        // Event-based summation
        if (_transmission && pop3._active){
            // Iterate over all incoming spikes (possibly delayed constantly)
            
            for(int _idx_j = 0; _idx_j < pop1.spiked.size(); _idx_j++){
                // Rank of the presynaptic neuron
                int rk_j = pop1.spiked[_idx_j];
                // Find the presynaptic neuron in the inverse connectivity matrix
                auto inv_post_ptr = inv_pre_rank.find(rk_j);
                if (inv_post_ptr == inv_pre_rank.end())
                    continue;
                // List of postsynaptic neurons receiving spikes from that neuron
                std::vector< std::pair<int, int> >& inv_post = inv_post_ptr->second;
                // Number of post neurons
                int nb_post = inv_post.size();
        
                
                // Iterate over connected post neurons
                for(int _idx_i = 0; _idx_i < nb_post; _idx_i++){
                    // Retrieve the correct indices
                    int i = inv_post[_idx_i].first;
                    int j = inv_post[_idx_i].second;
        
                    // Event-driven integration
                    
                    // Update conductance
                    
                    pop3.g_Exc[post_rank[i]] +=  w;
        
                    // Synaptic plasticity: pre-events
                    
                }
            }
        
            
        } // active
        
    }

    // Draws random numbers
    void update_rng() {

    }

    // Updates synaptic variables
    void update_synapse() {


    }

    // Post-synaptic events
    void post_event() {


    }

    // Accessors for default attributes
    int get_size() { return size; }
    void set_size(int new_size) { size = new_size; }

    // Additional access methods

    // Accessor to connectivity data
    std::vector<int> get_post_rank() { return post_rank; }
    void set_post_rank(std::vector<int> ranks) { post_rank = ranks; }
    std::vector< std::vector<int> > get_pre_rank() { return pre_rank; }
    void set_pre_rank(std::vector< std::vector<int> > ranks) { pre_rank = ranks; }
    int nb_synapses(int n) { return pre_rank[n].size(); }




    // Memory management
    long int size_in_bytes() {
        long int size_in_bytes = 0;
        // local parameter w
        size_in_bytes += sizeof(double);	// w
        
        return size_in_bytes;
    }

    void clear() {
    #ifdef _DEBUG
        std::cout << "PopStruct1::clear()" << std::endl;
    #endif
        // Variables
        
    }
};
