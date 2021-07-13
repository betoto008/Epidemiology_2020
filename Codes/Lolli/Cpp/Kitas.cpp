//
//  Kitas.cpp
//  
//
//  Created by Roberto Moran Tovar on 27.02.21.
//  Simulate the dynamics of infections in a kita.

#include "./lib/Kitas.hpp"

#include "./lib/random.cpp"


using namespace std;

int main(int argc, char* argv[]){
    
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
    
    gsl_rng_set(r, time(NULL));
    
    // argv-> 1:how many days ; 2:which days ; 3:tau ; 4:beta ; 5:p_in_inf
    std::string beta_s (argv[4]);
    std::string p_in_inf_s (argv[5]);
    
    
    std::cout<<">Running simulation of Kitas ..."<< std::endl;
    clock_t t1,t2;
    t1=clock();
    
    std::string Text_files_path = "../../../../../Dropbox/Research/Epidemiology_2020/Text_files/Kitas_Schools/";
    
    // Parameters
    int N = 200; //Number of kids
    int T = 12*7; //Total number of days for the simulation
    int t_inc = 4; //Incubation period in days
    int t_inf = 6; //Infectious period in days
    int tau = argv[3][0]-'0';
    double beta = std::stod(beta_s); //Infection rate days^-1
    double p_det = 0.9; //Probability of detection.
    double p_in_inf = std::stod(p_in_inf_s); //Probability of influx of a new infection to the kita. Should be proportional to the prevalence in the city.
    
    
    //Vector with the kids of a kita: 0 = Healthy, 1 = Incubation, 2 = Infectious, 3 = Recovered, 4 = Detected.
    std::vector < int > kids;
    kids.resize(N);
    for(int n = 0 ; n<N ; n++){
        kids[n] = 0;
    }
    //--------------------------------
    
    //Vectors with the counter of a incubations days, infectious days and actual infected days of each kid.
    
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
    for(int n = 0 ; n<N ; n++){
        kids[n] = 0;
        incubation[n] = t_inc;
        infectious[n] = t_inf;
        
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
    cout << endl;
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
    
    // initial number of kids in each compartment
    int H = N;
    int Inc = 0;
    int Inf = 0;
    int Rec = 0;
    int Det = 0;
    
    double inf;
    double r_inf;
    double r_det;
    
    std::vector<int>::iterator a;
    //kids[10] = 1;
    
    //Output file
    std::ofstream fout (Text_files_path+"output_testing_days-"+std::to_string(n_testing_days)+"-"+argv[2]+"_beta-"+std::to_string(beta)+"_pin-"+std::to_string(p_in_inf)+"_tau-"+std::to_string(tau)+".txt");

    // for-loop running over T days
    for(int t = 0; t<=T ; t++){
        int d = t%7;
        update_state_kids(N, &H, &Inc, &Inf, &Rec, &Det, &kids);
        fout << H << " " << Inc << " " << Inf << " " << Rec << " " << Det << std::endl;
        inf = Inf/double(N);
        //std::string day = days[d];
        
        if (d<5) { //Week days
            a = std::find(testing_days.begin(), testing_days.end(), d);
            // for-loop running over N kids
            for(int n = 0; n<N ; n++){
                if(kids[n]==0){ //Kid is healthy
                    r_inf = randX(0,1);
                    if(r_inf<(inf*beta)){
                        kids[n] = 1;
                    }
                } else{ //Kid isn't healthy
                    
                    if((a != testing_days.end())){ // Testing day?
                        if(kids[n] < 3){ // Kid isn't yet recovered or tested
                            r_det = randX(0,1);
                            if(r_det < (p_det*((1/(1+pow((log10(exp(3.0*infected_days[n]))/(2.5)),(-4))))))){ //Kid is detected
                                kids[n] = 4;
                            /*
                            if(infected_days[n] > (incubation[n]-tau)){ //Kid is detectable
                                r_det = randX(0,1);
                                if(r_det < p_det){ //Kid is detected
                                    kids[n] = 4;
                             
                                }
                             */
                            }
                        }
                        
                    }
                    
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
                }
                
            }
        }else { //Weekends
            // for-loop running over N kids
            for(int n = 0; n<N ; n++){
                if(kids[n]==0){ //Kid is healthy
                    r_inf = randX(0,1);
                    if(r_inf<(p_in_inf)){
                        kids[n] = 1;
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
    fout.close();
    //------------------------------------------------------------------------------
    std::cout<< ">Simulation completed…"<< std::endl;
    t2= clock();
    std::cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< std::endl;
    return 0;
}
