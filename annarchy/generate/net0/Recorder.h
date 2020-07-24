
/*
 * Recorders
 *
 */
class Monitor
{
public:
    Monitor(std::vector<int> ranks, int period, int period_offset, long int offset) {
        this->ranks = ranks;
        this->period_ = period;
        this->period_offset_ = period_offset;
        this->offset_ = offset;
        if(this->ranks.size() ==1 && this->ranks[0]==-1) // All neurons should be recorded
            this->partial = false;
        else
            this->partial = true;
    };

    ~Monitor() = default;

    virtual void record() = 0;
    virtual void record_targets() = 0;
    virtual long int size_in_bytes() = 0;
    virtual void clear() = 0;

    // Attributes
    bool partial;
    std::vector<int> ranks;
    int period_;
    int period_offset_;
    long int offset_;

};

class PopRecorder0 : public Monitor
{
public:
    PopRecorder0(std::vector<int> ranks, int period, int period_offset, long int offset)
        : Monitor(ranks, period, period_offset, offset) {

        this->rates = std::vector< std::vector< double > >();
        this->record_rates = false; 
        this->p = std::vector< std::vector< double > >();
        this->record_p = false; 
        this->r = std::vector< std::vector< double > >();
        this->record_r = false; 
        this->spike = std::map<int,  std::vector< long int > >();
        if(!this->partial){
            for(int i=0; i<pop0.size; i++) {
                this->spike[i]=std::vector<long int>();
            }
        }
        else{
            for(int i=0; i<this->ranks.size(); i++) {
                this->spike[this->ranks[i]]=std::vector<long int>();
            }
        }
        this->record_spike = false; 

    }

    void record() {

        if(this->record_rates && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->rates.push_back(pop0.rates);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop0.rates[this->ranks[i]]);
                }
                this->rates.push_back(tmp);
            }
        }
        if(this->record_p && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->p.push_back(pop0.p);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop0.p[this->ranks[i]]);
                }
                this->p.push_back(tmp);
            }
        }
        if(this->record_r && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->r.push_back(pop0.r);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop0.r[this->ranks[i]]);
                }
                this->r.push_back(tmp);
            }
        }
        if(this->record_spike){
            for(int i=0; i<pop0.spiked.size(); i++){
                if(!this->partial){
                    this->spike[pop0.spiked[i]].push_back(t);
                }
                else{
                    if( std::find(this->ranks.begin(), this->ranks.end(), pop0.spiked[i])!=this->ranks.end() ){
                        this->spike[pop0.spiked[i]].push_back(t);
                    }
                }
            }
        } 
    }

    void record_targets() {

    }

    long int size_in_bytes() {
        long int size_in_bytes = 0;
        size_in_bytes += sizeof(std::vector<double>) * rates.capacity();	//rates
        
        for(auto it=rates.begin(); it!= rates.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * p.capacity();	//p
        
        for(auto it=p.begin(); it!= p.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * r.capacity();	//r
        
        for(auto it=r.begin(); it!= r.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }
        return size_in_bytes;
    }

    void clear() {
    #ifdef _DEBUG
        std::cout << "PopRecorder0::clear()" << std::endl;
    #endif
        
                for(auto it = this->rates.begin(); it != this->rates.end(); it++)
                    it->clear();
                this->rates.clear();
            
                for(auto it = this->p.begin(); it != this->p.end(); it++)
                    it->clear();
                this->p.clear();
            
                for(auto it = this->r.begin(); it != this->r.end(); it++)
                    it->clear();
                this->r.clear();
            
    }



    // Local variable rates
    std::vector< std::vector< double > > rates ;
    bool record_rates ; 
    // Local variable p
    std::vector< std::vector< double > > p ;
    bool record_p ; 
    // Local variable r
    std::vector< std::vector< double > > r ;
    bool record_r ; 
    // Local variable spike
    std::map<int, std::vector< long int > > spike ;
    bool record_spike ;
    void clear_spike() {
        for ( auto it = spike.begin(); it != spike.end(); it++ ) {
            it->second.clear();
        }
    }

};

class PopRecorder1 : public Monitor
{
public:
    PopRecorder1(std::vector<int> ranks, int period, int period_offset, long int offset)
        : Monitor(ranks, period, period_offset, offset) {

        this->g_vm = std::vector< std::vector< double > >();
        this->record_g_vm = false; 
        this->EL = std::vector< double >();
        this->record_EL = false; 
        this->VTrest = std::vector< double >();
        this->record_VTrest = false; 
        this->taux = std::vector< double >();
        this->record_taux = false; 
        this->Spike = std::vector< std::vector< double > >();
        this->record_Spike = false; 
        this->Reset = std::vector< std::vector< double > >();
        this->record_Reset = false; 
        this->xtrace = std::vector< std::vector< double > >();
        this->record_xtrace = false; 
        this->state = std::vector< std::vector< double > >();
        this->record_state = false; 
        this->r = std::vector< std::vector< double > >();
        this->record_r = false; 
        this->spike = std::map<int,  std::vector< long int > >();
        if(!this->partial){
            for(int i=0; i<pop1.size; i++) {
                this->spike[i]=std::vector<long int>();
            }
        }
        else{
            for(int i=0; i<this->ranks.size(); i++) {
                this->spike[this->ranks[i]]=std::vector<long int>();
            }
        }
        this->record_spike = false; 

    }

    void record() {

        if(this->record_EL && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->EL.push_back(pop1.EL);
        } 
        if(this->record_VTrest && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->VTrest.push_back(pop1.VTrest);
        } 
        if(this->record_taux && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->taux.push_back(pop1.taux);
        } 
        if(this->record_Spike && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->Spike.push_back(pop1.Spike);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop1.Spike[this->ranks[i]]);
                }
                this->Spike.push_back(tmp);
            }
        }
        if(this->record_Reset && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->Reset.push_back(pop1.Reset);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop1.Reset[this->ranks[i]]);
                }
                this->Reset.push_back(tmp);
            }
        }
        if(this->record_xtrace && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->xtrace.push_back(pop1.xtrace);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop1.xtrace[this->ranks[i]]);
                }
                this->xtrace.push_back(tmp);
            }
        }
        if(this->record_state && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->state.push_back(pop1.state);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop1.state[this->ranks[i]]);
                }
                this->state.push_back(tmp);
            }
        }
        if(this->record_r && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->r.push_back(pop1.r);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop1.r[this->ranks[i]]);
                }
                this->r.push_back(tmp);
            }
        }
        if(this->record_spike){
            for(int i=0; i<pop1.spiked.size(); i++){
                if(!this->partial){
                    this->spike[pop1.spiked[i]].push_back(t);
                }
                else{
                    if( std::find(this->ranks.begin(), this->ranks.end(), pop1.spiked[i])!=this->ranks.end() ){
                        this->spike[pop1.spiked[i]].push_back(t);
                    }
                }
            }
        } 
    }

    void record_targets() {

        if(this->record_g_vm && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->g_vm.push_back(pop1.g_vm);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop1.g_vm[this->ranks[i]]);
                }
                this->g_vm.push_back(tmp);
            }
        }
    }

    long int size_in_bytes() {
        long int size_in_bytes = 0;
        size_in_bytes += sizeof(double);	//EL
        size_in_bytes += sizeof(double);	//VTrest
        size_in_bytes += sizeof(double);	//taux
        size_in_bytes += sizeof(std::vector<double>) * Spike.capacity();	//Spike
        
        for(auto it=Spike.begin(); it!= Spike.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * Reset.capacity();	//Reset
        
        for(auto it=Reset.begin(); it!= Reset.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * xtrace.capacity();	//xtrace
        
        for(auto it=xtrace.begin(); it!= xtrace.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * state.capacity();	//state
        
        for(auto it=state.begin(); it!= state.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * r.capacity();	//r
        
        for(auto it=r.begin(); it!= r.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }
        return size_in_bytes;
    }

    void clear() {
    #ifdef _DEBUG
        std::cout << "PopRecorder1::clear()" << std::endl;
    #endif
        
                this->EL.clear();
            
                this->VTrest.clear();
            
                this->taux.clear();
            
                for(auto it = this->Spike.begin(); it != this->Spike.end(); it++)
                    it->clear();
                this->Spike.clear();
            
                for(auto it = this->Reset.begin(); it != this->Reset.end(); it++)
                    it->clear();
                this->Reset.clear();
            
                for(auto it = this->xtrace.begin(); it != this->xtrace.end(); it++)
                    it->clear();
                this->xtrace.clear();
            
                for(auto it = this->state.begin(); it != this->state.end(); it++)
                    it->clear();
                this->state.clear();
            
                for(auto it = this->r.begin(); it != this->r.end(); it++)
                    it->clear();
                this->r.clear();
            
    }



    // Local variable g_vm
    std::vector< std::vector< double > > g_vm ;
    bool record_g_vm ; 
    // Global variable EL
    std::vector< double > EL ;
    bool record_EL ; 
    // Global variable VTrest
    std::vector< double > VTrest ;
    bool record_VTrest ; 
    // Global variable taux
    std::vector< double > taux ;
    bool record_taux ; 
    // Local variable Spike
    std::vector< std::vector< double > > Spike ;
    bool record_Spike ; 
    // Local variable Reset
    std::vector< std::vector< double > > Reset ;
    bool record_Reset ; 
    // Local variable xtrace
    std::vector< std::vector< double > > xtrace ;
    bool record_xtrace ; 
    // Local variable state
    std::vector< std::vector< double > > state ;
    bool record_state ; 
    // Local variable r
    std::vector< std::vector< double > > r ;
    bool record_r ; 
    // Local variable spike
    std::map<int, std::vector< long int > > spike ;
    bool record_spike ;
    void clear_spike() {
        for ( auto it = spike.begin(); it != spike.end(); it++ ) {
            it->second.clear();
        }
    }

};

class PopRecorder2 : public Monitor
{
public:
    PopRecorder2(std::vector<int> ranks, int period, int period_offset, long int offset)
        : Monitor(ranks, period, period_offset, offset) {

        this->g_Exc = std::vector< std::vector< double > >();
        this->record_g_Exc = false; 
        this->gL = std::vector< double >();
        this->record_gL = false; 
        this->DeltaT = std::vector< double >();
        this->record_DeltaT = false; 
        this->tauw = std::vector< double >();
        this->record_tauw = false; 
        this->a = std::vector< double >();
        this->record_a = false; 
        this->b = std::vector< double >();
        this->record_b = false; 
        this->EL = std::vector< double >();
        this->record_EL = false; 
        this->C = std::vector< double >();
        this->record_C = false; 
        this->tauz = std::vector< double >();
        this->record_tauz = false; 
        this->tauVT = std::vector< double >();
        this->record_tauVT = false; 
        this->Isp = std::vector< double >();
        this->record_Isp = false; 
        this->VTMax = std::vector< double >();
        this->record_VTMax = false; 
        this->VTrest = std::vector< double >();
        this->record_VTrest = false; 
        this->taux = std::vector< double >();
        this->record_taux = false; 
        this->tauLTD = std::vector< double >();
        this->record_tauLTD = false; 
        this->tauLTP = std::vector< double >();
        this->record_tauLTP = false; 
        this->taumean = std::vector< double >();
        this->record_taumean = false; 
        this->tau_gExc = std::vector< double >();
        this->record_tau_gExc = false; 
        this->inter_vm = std::vector< std::vector< double > >();
        this->record_inter_vm = false; 
        this->vm = std::vector< std::vector< double > >();
        this->record_vm = false; 
        this->vmean = std::vector< std::vector< double > >();
        this->record_vmean = false; 
        this->umeanLTD = std::vector< std::vector< double > >();
        this->record_umeanLTD = false; 
        this->umeanLTP = std::vector< std::vector< double > >();
        this->record_umeanLTP = false; 
        this->xtrace = std::vector< std::vector< double > >();
        this->record_xtrace = false; 
        this->wad = std::vector< std::vector< double > >();
        this->record_wad = false; 
        this->z = std::vector< std::vector< double > >();
        this->record_z = false; 
        this->VT = std::vector< std::vector< double > >();
        this->record_VT = false; 
        this->state = std::vector< std::vector< double > >();
        this->record_state = false; 
        this->Spike = std::vector< std::vector< double > >();
        this->record_Spike = false; 
        this->resetvar = std::vector< std::vector< double > >();
        this->record_resetvar = false; 
        this->r = std::vector< std::vector< double > >();
        this->record_r = false; 
        this->spike = std::map<int,  std::vector< long int > >();
        if(!this->partial){
            for(int i=0; i<pop2.size; i++) {
                this->spike[i]=std::vector<long int>();
            }
        }
        else{
            for(int i=0; i<this->ranks.size(); i++) {
                this->spike[this->ranks[i]]=std::vector<long int>();
            }
        }
        this->record_spike = false; 

    }

    void record() {

        if(this->record_gL && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->gL.push_back(pop2.gL);
        } 
        if(this->record_DeltaT && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->DeltaT.push_back(pop2.DeltaT);
        } 
        if(this->record_tauw && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->tauw.push_back(pop2.tauw);
        } 
        if(this->record_a && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->a.push_back(pop2.a);
        } 
        if(this->record_b && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->b.push_back(pop2.b);
        } 
        if(this->record_EL && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->EL.push_back(pop2.EL);
        } 
        if(this->record_C && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->C.push_back(pop2.C);
        } 
        if(this->record_tauz && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->tauz.push_back(pop2.tauz);
        } 
        if(this->record_tauVT && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->tauVT.push_back(pop2.tauVT);
        } 
        if(this->record_Isp && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->Isp.push_back(pop2.Isp);
        } 
        if(this->record_VTMax && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->VTMax.push_back(pop2.VTMax);
        } 
        if(this->record_VTrest && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->VTrest.push_back(pop2.VTrest);
        } 
        if(this->record_taux && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->taux.push_back(pop2.taux);
        } 
        if(this->record_tauLTD && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->tauLTD.push_back(pop2.tauLTD);
        } 
        if(this->record_tauLTP && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->tauLTP.push_back(pop2.tauLTP);
        } 
        if(this->record_taumean && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->taumean.push_back(pop2.taumean);
        } 
        if(this->record_tau_gExc && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->tau_gExc.push_back(pop2.tau_gExc);
        } 
        if(this->record_inter_vm && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->inter_vm.push_back(pop2.inter_vm);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.inter_vm[this->ranks[i]]);
                }
                this->inter_vm.push_back(tmp);
            }
        }
        if(this->record_vm && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->vm.push_back(pop2.vm);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.vm[this->ranks[i]]);
                }
                this->vm.push_back(tmp);
            }
        }
        if(this->record_vmean && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->vmean.push_back(pop2.vmean);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.vmean[this->ranks[i]]);
                }
                this->vmean.push_back(tmp);
            }
        }
        if(this->record_umeanLTD && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->umeanLTD.push_back(pop2.umeanLTD);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.umeanLTD[this->ranks[i]]);
                }
                this->umeanLTD.push_back(tmp);
            }
        }
        if(this->record_umeanLTP && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->umeanLTP.push_back(pop2.umeanLTP);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.umeanLTP[this->ranks[i]]);
                }
                this->umeanLTP.push_back(tmp);
            }
        }
        if(this->record_xtrace && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->xtrace.push_back(pop2.xtrace);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.xtrace[this->ranks[i]]);
                }
                this->xtrace.push_back(tmp);
            }
        }
        if(this->record_wad && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->wad.push_back(pop2.wad);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.wad[this->ranks[i]]);
                }
                this->wad.push_back(tmp);
            }
        }
        if(this->record_z && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->z.push_back(pop2.z);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.z[this->ranks[i]]);
                }
                this->z.push_back(tmp);
            }
        }
        if(this->record_VT && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->VT.push_back(pop2.VT);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.VT[this->ranks[i]]);
                }
                this->VT.push_back(tmp);
            }
        }
        if(this->record_state && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->state.push_back(pop2.state);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.state[this->ranks[i]]);
                }
                this->state.push_back(tmp);
            }
        }
        if(this->record_Spike && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->Spike.push_back(pop2.Spike);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.Spike[this->ranks[i]]);
                }
                this->Spike.push_back(tmp);
            }
        }
        if(this->record_resetvar && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->resetvar.push_back(pop2.resetvar);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.resetvar[this->ranks[i]]);
                }
                this->resetvar.push_back(tmp);
            }
        }
        if(this->record_r && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->r.push_back(pop2.r);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.r[this->ranks[i]]);
                }
                this->r.push_back(tmp);
            }
        }
        if(this->record_spike){
            for(int i=0; i<pop2.spiked.size(); i++){
                if(!this->partial){
                    this->spike[pop2.spiked[i]].push_back(t);
                }
                else{
                    if( std::find(this->ranks.begin(), this->ranks.end(), pop2.spiked[i])!=this->ranks.end() ){
                        this->spike[pop2.spiked[i]].push_back(t);
                    }
                }
            }
        } 
    }

    void record_targets() {

        if(this->record_g_Exc && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            if(!this->partial)
                this->g_Exc.push_back(pop2.g_Exc);
            else{
                std::vector<double> tmp = std::vector<double>();
                for (unsigned int i=0; i<this->ranks.size(); i++){
                    tmp.push_back(pop2.g_Exc[this->ranks[i]]);
                }
                this->g_Exc.push_back(tmp);
            }
        }
    }

    long int size_in_bytes() {
        long int size_in_bytes = 0;
        size_in_bytes += sizeof(double);	//gL
        size_in_bytes += sizeof(double);	//DeltaT
        size_in_bytes += sizeof(double);	//tauw
        size_in_bytes += sizeof(double);	//a
        size_in_bytes += sizeof(double);	//b
        size_in_bytes += sizeof(double);	//EL
        size_in_bytes += sizeof(double);	//C
        size_in_bytes += sizeof(double);	//tauz
        size_in_bytes += sizeof(double);	//tauVT
        size_in_bytes += sizeof(double);	//Isp
        size_in_bytes += sizeof(double);	//VTMax
        size_in_bytes += sizeof(double);	//VTrest
        size_in_bytes += sizeof(double);	//taux
        size_in_bytes += sizeof(double);	//tauLTD
        size_in_bytes += sizeof(double);	//tauLTP
        size_in_bytes += sizeof(double);	//taumean
        size_in_bytes += sizeof(double);	//tau_gExc
        size_in_bytes += sizeof(std::vector<double>) * inter_vm.capacity();	//inter_vm
        
        for(auto it=inter_vm.begin(); it!= inter_vm.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * vm.capacity();	//vm
        
        for(auto it=vm.begin(); it!= vm.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * vmean.capacity();	//vmean
        
        for(auto it=vmean.begin(); it!= vmean.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * umeanLTD.capacity();	//umeanLTD
        
        for(auto it=umeanLTD.begin(); it!= umeanLTD.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * umeanLTP.capacity();	//umeanLTP
        
        for(auto it=umeanLTP.begin(); it!= umeanLTP.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * xtrace.capacity();	//xtrace
        
        for(auto it=xtrace.begin(); it!= xtrace.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * wad.capacity();	//wad
        
        for(auto it=wad.begin(); it!= wad.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * z.capacity();	//z
        
        for(auto it=z.begin(); it!= z.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * VT.capacity();	//VT
        
        for(auto it=VT.begin(); it!= VT.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * state.capacity();	//state
        
        for(auto it=state.begin(); it!= state.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * Spike.capacity();	//Spike
        
        for(auto it=Spike.begin(); it!= Spike.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * resetvar.capacity();	//resetvar
        
        for(auto it=resetvar.begin(); it!= resetvar.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }size_in_bytes += sizeof(std::vector<double>) * r.capacity();	//r
        
        for(auto it=r.begin(); it!= r.end(); it++) {
            size_in_bytes += it->capacity() * sizeof(double);
        }
        return size_in_bytes;
    }

    void clear() {
    #ifdef _DEBUG
        std::cout << "PopRecorder2::clear()" << std::endl;
    #endif
        
                this->gL.clear();
            
                this->DeltaT.clear();
            
                this->tauw.clear();
            
                this->a.clear();
            
                this->b.clear();
            
                this->EL.clear();
            
                this->C.clear();
            
                this->tauz.clear();
            
                this->tauVT.clear();
            
                this->Isp.clear();
            
                this->VTMax.clear();
            
                this->VTrest.clear();
            
                this->taux.clear();
            
                this->tauLTD.clear();
            
                this->tauLTP.clear();
            
                this->taumean.clear();
            
                this->tau_gExc.clear();
            
                for(auto it = this->inter_vm.begin(); it != this->inter_vm.end(); it++)
                    it->clear();
                this->inter_vm.clear();
            
                for(auto it = this->vm.begin(); it != this->vm.end(); it++)
                    it->clear();
                this->vm.clear();
            
                for(auto it = this->vmean.begin(); it != this->vmean.end(); it++)
                    it->clear();
                this->vmean.clear();
            
                for(auto it = this->umeanLTD.begin(); it != this->umeanLTD.end(); it++)
                    it->clear();
                this->umeanLTD.clear();
            
                for(auto it = this->umeanLTP.begin(); it != this->umeanLTP.end(); it++)
                    it->clear();
                this->umeanLTP.clear();
            
                for(auto it = this->xtrace.begin(); it != this->xtrace.end(); it++)
                    it->clear();
                this->xtrace.clear();
            
                for(auto it = this->wad.begin(); it != this->wad.end(); it++)
                    it->clear();
                this->wad.clear();
            
                for(auto it = this->z.begin(); it != this->z.end(); it++)
                    it->clear();
                this->z.clear();
            
                for(auto it = this->VT.begin(); it != this->VT.end(); it++)
                    it->clear();
                this->VT.clear();
            
                for(auto it = this->state.begin(); it != this->state.end(); it++)
                    it->clear();
                this->state.clear();
            
                for(auto it = this->Spike.begin(); it != this->Spike.end(); it++)
                    it->clear();
                this->Spike.clear();
            
                for(auto it = this->resetvar.begin(); it != this->resetvar.end(); it++)
                    it->clear();
                this->resetvar.clear();
            
                for(auto it = this->r.begin(); it != this->r.end(); it++)
                    it->clear();
                this->r.clear();
            
    }



    // Local variable g_Exc
    std::vector< std::vector< double > > g_Exc ;
    bool record_g_Exc ; 
    // Global variable gL
    std::vector< double > gL ;
    bool record_gL ; 
    // Global variable DeltaT
    std::vector< double > DeltaT ;
    bool record_DeltaT ; 
    // Global variable tauw
    std::vector< double > tauw ;
    bool record_tauw ; 
    // Global variable a
    std::vector< double > a ;
    bool record_a ; 
    // Global variable b
    std::vector< double > b ;
    bool record_b ; 
    // Global variable EL
    std::vector< double > EL ;
    bool record_EL ; 
    // Global variable C
    std::vector< double > C ;
    bool record_C ; 
    // Global variable tauz
    std::vector< double > tauz ;
    bool record_tauz ; 
    // Global variable tauVT
    std::vector< double > tauVT ;
    bool record_tauVT ; 
    // Global variable Isp
    std::vector< double > Isp ;
    bool record_Isp ; 
    // Global variable VTMax
    std::vector< double > VTMax ;
    bool record_VTMax ; 
    // Global variable VTrest
    std::vector< double > VTrest ;
    bool record_VTrest ; 
    // Global variable taux
    std::vector< double > taux ;
    bool record_taux ; 
    // Global variable tauLTD
    std::vector< double > tauLTD ;
    bool record_tauLTD ; 
    // Global variable tauLTP
    std::vector< double > tauLTP ;
    bool record_tauLTP ; 
    // Global variable taumean
    std::vector< double > taumean ;
    bool record_taumean ; 
    // Global variable tau_gExc
    std::vector< double > tau_gExc ;
    bool record_tau_gExc ; 
    // Local variable inter_vm
    std::vector< std::vector< double > > inter_vm ;
    bool record_inter_vm ; 
    // Local variable vm
    std::vector< std::vector< double > > vm ;
    bool record_vm ; 
    // Local variable vmean
    std::vector< std::vector< double > > vmean ;
    bool record_vmean ; 
    // Local variable umeanLTD
    std::vector< std::vector< double > > umeanLTD ;
    bool record_umeanLTD ; 
    // Local variable umeanLTP
    std::vector< std::vector< double > > umeanLTP ;
    bool record_umeanLTP ; 
    // Local variable xtrace
    std::vector< std::vector< double > > xtrace ;
    bool record_xtrace ; 
    // Local variable wad
    std::vector< std::vector< double > > wad ;
    bool record_wad ; 
    // Local variable z
    std::vector< std::vector< double > > z ;
    bool record_z ; 
    // Local variable VT
    std::vector< std::vector< double > > VT ;
    bool record_VT ; 
    // Local variable state
    std::vector< std::vector< double > > state ;
    bool record_state ; 
    // Local variable Spike
    std::vector< std::vector< double > > Spike ;
    bool record_Spike ; 
    // Local variable resetvar
    std::vector< std::vector< double > > resetvar ;
    bool record_resetvar ; 
    // Local variable r
    std::vector< std::vector< double > > r ;
    bool record_r ; 
    // Local variable spike
    std::map<int, std::vector< long int > > spike ;
    bool record_spike ;
    void clear_spike() {
        for ( auto it = spike.begin(); it != spike.end(); it++ ) {
            it->second.clear();
        }
    }

};

class ProjRecorder0 : public Monitor
{
public:
    ProjRecorder0(std::vector<int> ranks, int period, int period_offset, long int offset)
        : Monitor(ranks, period, period_offset, offset)
    {

        this->w = std::vector< double >();
        this->record_w = false;

    };

    void record() {

        if(this->record_w && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->w.push_back(proj0.w);
        }

    };

    void record_targets() { /* nothing to do here */ }
    long int size_in_bytes() {
        std::cout << "ProjMonitor::size_in_bytes(): not implemented for openMP paradigm." << std::endl;
        return 0;
    }

    void clear() {
        std::cout << "PopMonitor0::clear(): not implemented for openMP paradigm." << std::endl;
    }


    // Global variable w
    std::vector< double > w ;
    bool record_w ;

};

class ProjRecorder1 : public Monitor
{
public:
    ProjRecorder1(std::vector<int> ranks, int period, int period_offset, long int offset)
        : Monitor(ranks, period, period_offset, offset)
    {

        this->vmean_fix = std::vector< double >();
        this->record_vmean_fix = false;

        this->urefsquare = std::vector< double >();
        this->record_urefsquare = false;

        this->thetaLTD = std::vector< double >();
        this->record_thetaLTD = false;

        this->thetaLTP = std::vector< double >();
        this->record_thetaLTP = false;

        this->aLTD = std::vector< double >();
        this->record_aLTD = false;

        this->aLTP = std::vector< double >();
        this->record_aLTP = false;

        this->wMin = std::vector< double >();
        this->record_wMin = false;

        this->wMax = std::vector< double >();
        this->record_wMax = false;

        this->transmit = std::vector< double >();
        this->record_transmit = false;

        this->set_fix = std::vector< double >();
        this->record_set_fix = false;

        this->ltdTerm_fix = std::vector< std::vector< std::vector< double > > >();
        this->record_ltdTerm_fix = false;

        this->ltdTerm = std::vector< std::vector< std::vector< double > > >();
        this->record_ltdTerm = false;

        this->ltpTerm = std::vector< std::vector< std::vector< double > > >();
        this->record_ltpTerm = false;

        this->deltaW = std::vector< std::vector< std::vector< double > > >();
        this->record_deltaW = false;

        this->w = std::vector< std::vector< std::vector< double > > >();
        this->record_w = false;

    };

    void record() {

        if(this->record_vmean_fix && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->vmean_fix.push_back(proj1.vmean_fix);
        }

        if(this->record_urefsquare && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->urefsquare.push_back(proj1.urefsquare);
        }

        if(this->record_thetaLTD && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->thetaLTD.push_back(proj1.thetaLTD);
        }

        if(this->record_thetaLTP && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->thetaLTP.push_back(proj1.thetaLTP);
        }

        if(this->record_aLTD && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->aLTD.push_back(proj1.aLTD);
        }

        if(this->record_aLTP && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->aLTP.push_back(proj1.aLTP);
        }

        if(this->record_wMin && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->wMin.push_back(proj1.wMin);
        }

        if(this->record_wMax && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->wMax.push_back(proj1.wMax);
        }

        if(this->record_transmit && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->transmit.push_back(proj1.transmit);
        }

        if(this->record_set_fix && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            this->set_fix.push_back(proj1.set_fix);
        }

        if(this->record_ltdTerm_fix && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            std::vector< std::vector< double > > tmp;
            for(int i=0; i<this->ranks.size(); i++){
                tmp.push_back(proj1.ltdTerm_fix[this->ranks[i]]);
            }
            this->ltdTerm_fix.push_back(tmp);
            tmp.clear();
        }

        if(this->record_ltdTerm && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            std::vector< std::vector< double > > tmp;
            for(int i=0; i<this->ranks.size(); i++){
                tmp.push_back(proj1.ltdTerm[this->ranks[i]]);
            }
            this->ltdTerm.push_back(tmp);
            tmp.clear();
        }

        if(this->record_ltpTerm && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            std::vector< std::vector< double > > tmp;
            for(int i=0; i<this->ranks.size(); i++){
                tmp.push_back(proj1.ltpTerm[this->ranks[i]]);
            }
            this->ltpTerm.push_back(tmp);
            tmp.clear();
        }

        if(this->record_deltaW && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            std::vector< std::vector< double > > tmp;
            for(int i=0; i<this->ranks.size(); i++){
                tmp.push_back(proj1.deltaW[this->ranks[i]]);
            }
            this->deltaW.push_back(tmp);
            tmp.clear();
        }

        if(this->record_w && ( (t - this->offset_) % this->period_ == this->period_offset_ )){
            std::vector< std::vector< double > > tmp;
            for(int i=0; i<this->ranks.size(); i++){
                tmp.push_back(proj1.w[this->ranks[i]]);
            }
            this->w.push_back(tmp);
            tmp.clear();
        }

    };

    void record_targets() { /* nothing to do here */ }
    long int size_in_bytes() {
        std::cout << "ProjMonitor::size_in_bytes(): not implemented for openMP paradigm." << std::endl;
        return 0;
    }

    void clear() {
        std::cout << "PopMonitor1::clear(): not implemented for openMP paradigm." << std::endl;
    }


    // Global variable vmean_fix
    std::vector< double > vmean_fix ;
    bool record_vmean_fix ;

    // Global variable urefsquare
    std::vector< double > urefsquare ;
    bool record_urefsquare ;

    // Global variable thetaLTD
    std::vector< double > thetaLTD ;
    bool record_thetaLTD ;

    // Global variable thetaLTP
    std::vector< double > thetaLTP ;
    bool record_thetaLTP ;

    // Global variable aLTD
    std::vector< double > aLTD ;
    bool record_aLTD ;

    // Global variable aLTP
    std::vector< double > aLTP ;
    bool record_aLTP ;

    // Global variable wMin
    std::vector< double > wMin ;
    bool record_wMin ;

    // Global variable wMax
    std::vector< double > wMax ;
    bool record_wMax ;

    // Global variable transmit
    std::vector< double > transmit ;
    bool record_transmit ;

    // Global variable set_fix
    std::vector< double > set_fix ;
    bool record_set_fix ;

    // Local variable ltdTerm_fix
    std::vector< std::vector< std::vector< double > > > ltdTerm_fix ;
    bool record_ltdTerm_fix ;

    // Local variable ltdTerm
    std::vector< std::vector< std::vector< double > > > ltdTerm ;
    bool record_ltdTerm ;

    // Local variable ltpTerm
    std::vector< std::vector< std::vector< double > > > ltpTerm ;
    bool record_ltpTerm ;

    // Local variable deltaW
    std::vector< std::vector< std::vector< double > > > deltaW ;
    bool record_deltaW ;

    // Local variable w
    std::vector< std::vector< std::vector< double > > > w ;
    bool record_w ;

};

