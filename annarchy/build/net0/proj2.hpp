#pragma once

#include "pop2.hpp"
#include "pop3.hpp"



extern PopStruct2 pop2;
extern PopStruct3 pop3;


/////////////////////////////////////////////////////////////////////////////
// proj2: N1 -> N2 with target Exc
/////////////////////////////////////////////////////////////////////////////
struct ProjStruct2{
    // Number of dendrites
    int size;

    // Transmission and plasticity flags
    bool _transmission, _plasticity, _update;
    int _update_period;
    long int _update_offset;


    // Connectivity
    std::vector<int> post_rank;
    std::vector< std::vector< int > > pre_rank;

    // LIL weights
    std::vector< std::vector< double > > w;


    std::map< int, std::vector< std::pair<int, int> > > inv_pre_rank ;
    std::vector< int > inv_post_rank ;





    // Global parameter vmean_fix
    double  vmean_fix ;

    // Global parameter urefsquare
    double  urefsquare ;

    // Global parameter thetaLTD
    double  thetaLTD ;

    // Global parameter thetaLTP
    double  thetaLTP ;

    // Global parameter aLTD
    double  aLTD ;

    // Global parameter aLTP
    double  aLTP ;

    // Global parameter wMin
    double  wMin ;

    // Global parameter wMax
    double  wMax ;

    // Global parameter transmit
    double  transmit ;

    // Global parameter set_fix
    double  set_fix ;

    // Local variable ltdTerm_fix
    std::vector< std::vector<double > > ltdTerm_fix;

    // Local variable ltdTerm
    std::vector< std::vector<double > > ltdTerm;

    // Local variable ltpTerm
    std::vector< std::vector<double > > ltpTerm;

    // Local variable deltaW
    std::vector< std::vector<double > > deltaW;




    // Method called to initialize the projection
    void init_projection() {
        _transmission = true;
        _update = true;
        _plasticity = true;
        _update_period = 1;
        _update_offset = 0L;





        // Inverse the connectivity matrix if spiking neurons
        inverse_connectivity_matrix();



        // Global parameter vmean_fix
        vmean_fix = 0.0;

        // Global parameter urefsquare
        urefsquare = 0.0;

        // Global parameter thetaLTD
        thetaLTD = 0.0;

        // Global parameter thetaLTP
        thetaLTP = 0.0;

        // Global parameter aLTD
        aLTD = 0.0;

        // Global parameter aLTP
        aLTP = 0.0;

        // Global parameter wMin
        wMin = 0.0;

        // Global parameter wMax
        wMax = 0.0;

        // Global parameter transmit
        transmit = 0.0;

        // Global parameter set_fix
        set_fix = 0.0;

        // Local variable ltdTerm_fix
        ltdTerm_fix = std::vector< std::vector<double> >(post_rank.size(), std::vector<double>());

        // Local variable ltdTerm
        ltdTerm = std::vector< std::vector<double> >(post_rank.size(), std::vector<double>());

        // Local variable ltpTerm
        ltpTerm = std::vector< std::vector<double> >(post_rank.size(), std::vector<double>());

        // Local variable deltaW
        deltaW = std::vector< std::vector<double> >(post_rank.size(), std::vector<double>());





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
            
            for(int _idx_j = 0; _idx_j < pop2.spiked.size(); _idx_j++){
                // Rank of the presynaptic neuron
                int rk_j = pop2.spiked[_idx_j];
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
                    
                    pop3.g_Exc[post_rank[i]] +=  transmit*w[i][j];
        
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

        int rk_post, rk_pre;
        double _dt = dt * _update_period;
        
        // Check periodicity
        if(_transmission && _update && pop3._active && ( (t - _update_offset)%_update_period == 0L) ){
            // Global variables
            
            // Local variables
            
            for(int i = 0; i < post_rank.size(); i++){
                rk_post = post_rank[i]; // Get postsynaptic rank
                // Semi-global variables
                
                // Local variables
                for(int j = 0; j < pre_rank[i].size(); j++){
                    rk_pre = pre_rank[i][j]; // Get presynaptic rank
                        
                    // ltdTerm_fix = if w>wMin : (aLTD*(vmean_fix/urefsquare)*pre.Spike * pos(post.umeanLTD - thetaLTD)) else : 0.0
                    ltdTerm_fix[i][j] = (w[i][j] > wMin ? pop2.Spike[rk_pre]*aLTD*vmean_fix*positive(pop3.umeanLTD[rk_post] - thetaLTD)/urefsquare : 0.0);
                    
                    
                    // ltdTerm = if w>wMin : (aLTD*(post.vmean/urefsquare)*pre.Spike * pos(post.umeanLTD - thetaLTD)) else : 0.0
                    ltdTerm[i][j] = (w[i][j] > wMin ? pop3.vmean[rk_post]*pop2.Spike[rk_pre]*aLTD*positive(pop3.umeanLTD[rk_post] - thetaLTD)/urefsquare : 0.0);
                    
                    
                    // ltpTerm = if w<wMax : (aLTP * pos(post.vm - thetaLTP) *(pre.xtrace)* pos(post.umeanLTP - thetaLTD)) else : 0.0
                    ltpTerm[i][j] = (w[i][j] < wMax ? pop2.xtrace[rk_pre]*aLTP*positive(pop3.umeanLTP[rk_post] - thetaLTD)*positive(pop3.vm[rk_post] - thetaLTP) : 0.0);
                    
                    
                    // deltaW = if set_fix==1: ltpTerm - ltdTerm_fix else: ltpTerm - ltdTerm
                    deltaW[i][j] = (set_fix == 1 ? -ltdTerm_fix[i][j] + ltpTerm[i][j] : -ltdTerm[i][j] + ltpTerm[i][j]);
                    
                    
                    // dw/_dt = deltaW
                    double _w = deltaW[i][j];
                    
                    // dw/_dt = deltaW
                    if(_plasticity){
                    w[i][j] += _dt*_w ;
                    if(w[i][j] < 0.0)
                        w[i][j] = 0.0;
                    if(w[i][j] > wMax)
                        w[i][j] = wMax;
                    
                    }
                    
                }
            }
        }
        
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

    // Local parameter w
    std::vector<std::vector< double > > get_w() {
        std::vector< std::vector< double > > w_new(w.size(), std::vector<double>());
        for(int i = 0; i < w.size(); i++) {
            w_new[i] = std::vector<double>(w[i].begin(), w[i].end());
        }
        return w_new;
    }
    std::vector< double > get_dendrite_w(int rk) { return std::vector<double>(w[rk].begin(), w[rk].end()); }
    double get_synapse_w(int rk_post, int rk_pre) { return w[rk_post][rk_pre]; }
    void set_w(std::vector<std::vector< double > >value) {
        w = std::vector< std::vector<double> >( value.size(), std::vector<double>() );
        for(int i = 0; i < value.size(); i++) {
            w[i] = std::vector<double>(value[i].begin(), value[i].end());
        }
    }
    void set_dendrite_w(int rk, std::vector< double > value) { w[rk] = std::vector<double>(value.begin(), value.end()); }
    void set_synapse_w(int rk_post, int rk_pre, double value) { w[rk_post][rk_pre] = value; }


    // Global parameter vmean_fix
    double get_vmean_fix() { return vmean_fix; }
    void set_vmean_fix(double value) { vmean_fix = value; }

    // Global parameter urefsquare
    double get_urefsquare() { return urefsquare; }
    void set_urefsquare(double value) { urefsquare = value; }

    // Global parameter thetaLTD
    double get_thetaLTD() { return thetaLTD; }
    void set_thetaLTD(double value) { thetaLTD = value; }

    // Global parameter thetaLTP
    double get_thetaLTP() { return thetaLTP; }
    void set_thetaLTP(double value) { thetaLTP = value; }

    // Global parameter aLTD
    double get_aLTD() { return aLTD; }
    void set_aLTD(double value) { aLTD = value; }

    // Global parameter aLTP
    double get_aLTP() { return aLTP; }
    void set_aLTP(double value) { aLTP = value; }

    // Global parameter wMin
    double get_wMin() { return wMin; }
    void set_wMin(double value) { wMin = value; }

    // Global parameter wMax
    double get_wMax() { return wMax; }
    void set_wMax(double value) { wMax = value; }

    // Global parameter transmit
    double get_transmit() { return transmit; }
    void set_transmit(double value) { transmit = value; }

    // Global parameter set_fix
    double get_set_fix() { return set_fix; }
    void set_set_fix(double value) { set_fix = value; }

    // Local variable ltdTerm_fix
    std::vector<std::vector< double > > get_ltdTerm_fix() { return ltdTerm_fix; }
    std::vector<double> get_dendrite_ltdTerm_fix(int rk) { return ltdTerm_fix[rk]; }
    double get_synapse_ltdTerm_fix(int rk_post, int rk_pre) { return ltdTerm_fix[rk_post][rk_pre]; }
    void set_ltdTerm_fix(std::vector<std::vector< double > >value) { ltdTerm_fix = value; }
    void set_dendrite_ltdTerm_fix(int rk, std::vector<double> value) { ltdTerm_fix[rk] = value; }
    void set_synapse_ltdTerm_fix(int rk_post, int rk_pre, double value) { ltdTerm_fix[rk_post][rk_pre] = value; }

    // Local variable ltdTerm
    std::vector<std::vector< double > > get_ltdTerm() { return ltdTerm; }
    std::vector<double> get_dendrite_ltdTerm(int rk) { return ltdTerm[rk]; }
    double get_synapse_ltdTerm(int rk_post, int rk_pre) { return ltdTerm[rk_post][rk_pre]; }
    void set_ltdTerm(std::vector<std::vector< double > >value) { ltdTerm = value; }
    void set_dendrite_ltdTerm(int rk, std::vector<double> value) { ltdTerm[rk] = value; }
    void set_synapse_ltdTerm(int rk_post, int rk_pre, double value) { ltdTerm[rk_post][rk_pre] = value; }

    // Local variable ltpTerm
    std::vector<std::vector< double > > get_ltpTerm() { return ltpTerm; }
    std::vector<double> get_dendrite_ltpTerm(int rk) { return ltpTerm[rk]; }
    double get_synapse_ltpTerm(int rk_post, int rk_pre) { return ltpTerm[rk_post][rk_pre]; }
    void set_ltpTerm(std::vector<std::vector< double > >value) { ltpTerm = value; }
    void set_dendrite_ltpTerm(int rk, std::vector<double> value) { ltpTerm[rk] = value; }
    void set_synapse_ltpTerm(int rk_post, int rk_pre, double value) { ltpTerm[rk_post][rk_pre] = value; }

    // Local variable deltaW
    std::vector<std::vector< double > > get_deltaW() { return deltaW; }
    std::vector<double> get_dendrite_deltaW(int rk) { return deltaW[rk]; }
    double get_synapse_deltaW(int rk_post, int rk_pre) { return deltaW[rk_post][rk_pre]; }
    void set_deltaW(std::vector<std::vector< double > >value) { deltaW = value; }
    void set_dendrite_deltaW(int rk, std::vector<double> value) { deltaW[rk] = value; }
    void set_synapse_deltaW(int rk_post, int rk_pre, double value) { deltaW[rk_post][rk_pre] = value; }



    // Memory management
    long int size_in_bytes() {
        long int size_in_bytes = 0;
        // local variable ltdTerm_fix
        size_in_bytes += sizeof(double) * ltdTerm_fix.capacity();
        for(auto it = ltdTerm_fix.begin(); it != ltdTerm_fix.end(); it++)
            size_in_bytes += (it->capacity()) * sizeof(double);
        // local variable ltdTerm
        size_in_bytes += sizeof(double) * ltdTerm.capacity();
        for(auto it = ltdTerm.begin(); it != ltdTerm.end(); it++)
            size_in_bytes += (it->capacity()) * sizeof(double);
        // local variable ltpTerm
        size_in_bytes += sizeof(double) * ltpTerm.capacity();
        for(auto it = ltpTerm.begin(); it != ltpTerm.end(); it++)
            size_in_bytes += (it->capacity()) * sizeof(double);
        // local variable deltaW
        size_in_bytes += sizeof(double) * deltaW.capacity();
        for(auto it = deltaW.begin(); it != deltaW.end(); it++)
            size_in_bytes += (it->capacity()) * sizeof(double);
        // local variable w
        size_in_bytes += sizeof(double) * w.capacity();
        for(auto it = w.begin(); it != w.end(); it++)
            size_in_bytes += (it->capacity()) * sizeof(double);
        // global parameter vmean_fix
        size_in_bytes += sizeof(double);	// vmean_fix
        // global parameter urefsquare
        size_in_bytes += sizeof(double);	// urefsquare
        // global parameter thetaLTD
        size_in_bytes += sizeof(double);	// thetaLTD
        // global parameter thetaLTP
        size_in_bytes += sizeof(double);	// thetaLTP
        // global parameter aLTD
        size_in_bytes += sizeof(double);	// aLTD
        // global parameter aLTP
        size_in_bytes += sizeof(double);	// aLTP
        // global parameter wMin
        size_in_bytes += sizeof(double);	// wMin
        // global parameter wMax
        size_in_bytes += sizeof(double);	// wMax
        // global parameter transmit
        size_in_bytes += sizeof(double);	// transmit
        // global parameter set_fix
        size_in_bytes += sizeof(double);	// set_fix
        
        return size_in_bytes;
    }

    void clear() {
    #ifdef _DEBUG
        std::cout << "PopStruct2::clear()" << std::endl;
    #endif
        // Variables
        ltdTerm_fix.clear();
        ltdTerm_fix.shrink_to_fit();
        ltdTerm.clear();
        ltdTerm.shrink_to_fit();
        ltpTerm.clear();
        ltpTerm.shrink_to_fit();
        deltaW.clear();
        deltaW.shrink_to_fit();
        w.clear();
        w.shrink_to_fit();
        
    }
};
