//
//  Kitas_ensemble.cpp
//  
//
//  Created by Roberto Moran Tovar on 03.03.21.
//  Simulate an ensemble of kita.

#include "./lib/Kitas.hpp"

#include "./lib/random.cpp"

using namespace std;

int main(int argc, char* argv[]){
    
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
    
    // argv-> 1:how many days ; 2:which days ; 3:tau ; 4:beta ; 5:p_in_inf
    std::string beta_s (argv[4]);
    std::string p_in_inf_s (argv[5]);
    
    
    std::cout<<">Running simulation of Kitas ..."<< std::endl;
    clock_t t1,t2;
    t1=clock();
    
    std::string Text_files_path = "../../../../../Dropbox/Research/Epidemiology_2020/Text_files/Kitas_Schools/Statistics/";
    
    // Parameters
    int n_ensemble = 100000; //Ensemble size
    int N = 20; //Number of kids
    int T = 3*7; //Total number of days for the simulation
    int tau = argv[3][0]-'0';
    int t_inc = 3; //Incubation period in days
    int t_inf = 5; //Infectious period in days
    double beta = std::stod(beta_s); //Infection rate days^-1
    double p_det = 0.98; //Probability of detection.
    double p_in_inf = std::stod(p_in_inf_s); //Probability of influx of a new infection to the kita. Should be proportional to the prevalence in the city.
    
    
    //Vector with the kids of a kita: 0 = Healthy, 1 = Incubation, 2 = Infectious, 3 = Recovered, 4 = Detected.
    std::vector < int > kids;
    kids.resize(N);
    for(int n = 0 ; n<N ; n++){
        kids[n] = 0;
    }
    //--------------------------------
    
    //days of the week
    std::string days_array[7] = {"mon", "tue", "wed", "thu", "fri", "sat", "sun"};
    std::vector< std::string > days;
    days.resize(7);
    
    for(int d = 0 ; d<7 ; d++){
        days[d] = days_array[d];
        //std::cout << d << ":" << days[d] << "\t";
    }
    std::cout << std::endl;
    
    std::vector<int> testing_days;
    int n_testing_days = argv[1][0] - '0';
    std::cout << "We are testing " << n_testing_days << " days: ";
    //std::cout << "How many days do you want to test per week? 1, 2, or 3?" << std::endl;
    //std::cin >> n_testing_days;
    if (n_testing_days!=0) {
        for (int d = 0; d<n_testing_days; d++){
            int m = argv[2][d] - '0';
            testing_days.push_back(m);
            std::cout << m << "(" << days[m] << ")" << ", ";
        }
    }
    std::cout << std::endl;
    
    std::cout << "beta is " << beta << " ,p_in is " << p_in_inf << " and tau is " << tau << std::endl << std::endl;
    
    //for (int d = 1; d<=n_testing_days; d++){
    //    int m;
    //    std::cout << "What is the day " << d<<"?"<< std::endl;
    //    std::cin >> m;
    //    testing_days.push_back(m);
    //}
    //--------------------------------
    
    // initial number of kids in each compartment
    int H = N;
    int Inc = 0;
    int Inf = 0;
    int Rec = 0;
    int Det = 0;
    
    int total_inf = 0;
    int total_det = 0;
    int total_trans = 0;
    int total_det_trans = 0;
    
    double inf;
    double r_inf;
    double r_inf2;
    double r_det;
    std::vector<int>::iterator a;
    //kids[10] = 1;
    
    //Output file
    std::ofstream fout (Text_files_path+"statistics_days-"+std::to_string(n_testing_days)+"-"+argv[2]+"_beta-"+std::to_string(beta)+"_pin-"+std::to_string(p_in_inf)+"_tau-"+std::to_string(tau)+"_1.txt");
    
    //std::ofstream fout_inc (Text_files_path+"tau_inc_days-"+std::to_string(n_testing_days)+"-"+argv[2]+"_beta-"+std::to_string(beta)+"_pin-"+std::to_string(p_in_inf)+"_tau-"+std::to_string(tau)+".txt");
    
    //std::ofstream fout_inf (Text_files_path+"tau_inf_days-"+std::to_string(n_testing_days)+"-"+argv[2]+"_beta-"+std::to_string(beta)+"_pin-"+std::to_string(p_in_inf)+"_tau-"+std::to_string(tau)+".txt");
    
    for(int j = 0 ; j < n_ensemble ; j++){
        
        gsl_rng_set(r, time(NULL));
        //------------ initialize vectors-------------
        
        //Distribution for times
        std::default_random_engine generator;
        double t_inc_n = 0;
        double t_inf_n = 0;
        
        std::vector < int > incubation;
        incubation.resize(N);
        std::vector < int > infectious;
        infectious.resize(N);
        std::vector < int > infected_days;
        infected_days.resize(N);
        std::vector < int > internal_transmission;
        internal_transmission.resize(N);
        
        // Try filling the arrays with constant values and exponentially distributed random values
        for(int n = 0 ; n<N ; n++){
            kids[n] = 0;
            incubation[n] = t_inc;
            infectious[n] = t_inf;
            infected_days[n] = 0;
            internal_transmission[n] = 0;
            
            /*
            while( t_inc_n < 1.0){
                t_inc_n =  gsl_ran_exponential (r, t_inc);
            }
            while( t_inf_n < 1.0){
                t_inf_n =  gsl_ran_exponential (r, t_inf);
            }
            incubation[n] = t_inc_n;
            fout_inc << t_inc_n << endl;
            infectious[n] = t_inf_n;
            fout_inf << t_inf_n << endl;
            infected_days[n] = 0;
            
            t_inc_n = 0;
            t_inf_n = 0;
            */
        }
        //--------------------------------
        
        H = N;
        Inc = 0;
        Inf = 0;
        Rec = 0;
        Det = 0;
        total_inf = 0;
        total_det = 0;
        total_trans = 0;
        total_det_trans = 0;
        //--------------------------------------------
        
        // for-loop running over T days
        for(int t = 0; t<=T ; t++){
            int d = t%7;
            update_state_kids(N, &H, &Inc, &Inf, &Rec, &Det, &kids);
            inf = Inf/double(N);
            std::string day = days[d];
            
            if (d<5) { //Week days
                a = std::find(testing_days.begin(), testing_days.end(), d);
                // for-loop running over N kids
                for(int n = 0; n<N ; n++){
                    if(kids[n]==0){ //Kid is healthy
                        r_inf = randX(0,1);
                        if(r_inf<(inf*beta)){
                            kids[n] = 1;
                            total_inf ++;
                            total_trans ++;
                            internal_transmission[n] = 1;
                        }else{
                            r_inf2 = randX(0,1);
                            if(r_inf2<(p_in_inf*0.5*beta)){
                                kids[n] = 1;
                                total_inf ++;
                            }
                        }
                    } else{ //Kid isn't healthy
                        
                        if(kids[n]==1){ //Kid is incubation period
                            if (infected_days[n]>incubation[n]) { //Incubation period is over?
                                kids[n]=2;
                            }
                            infected_days[n] = infected_days[n]+1;
                        } else if(kids[n]==2){ //Kid is infectious period
                            if (infected_days[n]<(incubation[n]+infectious[n])) { //Incubation period isn't over?
                                infected_days[n] = infected_days[n]+1;
                            }else{
                                kids[n]=3;
                            }
                        }
                        
                        if((a != testing_days.end())){ // Testing day?
                            if(kids[n] < 3){ // Kid isn't yet recovered or tested
                                r_det = randX(0,1);
                                if(r_det < (p_det*(det_rate(infected_days[n], t_inc, tau)))){ //Kid is detected
                                    kids[n] = 4;
                                    total_det++;
                                    if (internal_transmission[n]==1) {
                                        total_det_trans++;
                                    }
                                /*
                                if(infected_days[n] > (incubation[n]-tau)){ //Kid is detectable
                                    r_det = randX(0,1);
                                    if(r_det < p_det){ //Kid is detected
                                        kids[n] = 4;
                                        total_det++;
                                 
                                    }
                                 */
                                }
                            }
                        }
                    }
                }
            }else { //Weekends
                // for-loop running over N kids
                for(int n = 0; n<N ; n++){
                    if(kids[n]==0){ //Kid is healthy
                        r_inf2 = randX(0,1);
                        if(r_inf2<(p_in_inf*beta)){
                            kids[n] = 1;
                            total_inf ++;
                        }
                    } else{ //Kid isn't healthy
                        
                        if(kids[n]==1){ //Kid is incubation period
                            if (infected_days[n]>incubation[n]) { //Incubation period is over?
                                kids[n]=2;
                            }
                            infected_days[n] = infected_days[n]+1;
                        } else if(kids[n]==2){ //Kid is in infectious period
                            if (infected_days[n]<(incubation[n]+infectious[n])) { //Incubation period isn't over?
                                infected_days[n] = infected_days[n]+1;
                            }else{
                                kids[n]=3;
                            }
                        }
                    }
                }
            }
          
        }
        
        fout << total_inf << "\t" << total_det << "\t" << total_trans << "\t" << total_det_trans << std::endl;
    }
    fout.close();
    //fout_inc.close();
    //fout_inf.close();
    //------------------------------------------------------------------------------
    std::cout<< ">Simulation completed…"<< std::endl;
    t2= clock();
    std::cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< std::endl;
    return 0;
}
