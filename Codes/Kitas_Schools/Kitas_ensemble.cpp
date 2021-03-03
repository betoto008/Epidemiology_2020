//
//  Kitas_ensemble.cpp
//  
//
//  Created by Roberto Moran Tovar on 03.03.21.
//  Simulate an ensemble of kita.

#include "./lib/Kitas.hpp"

#include "./lib/random.cpp"

//using namespace std;

int main(int argc, char* argv[]){
    
    // argv-> 1:how many days ; 2:which days ; 3:t_inc ; 4:beta ; 5:p_in_inf
    std::string beta_s (argv[4]);
    std::string p_in_inf_s (argv[5]);
    
    
    std::cout<<">Running simulation of Kitas ..."<< std::endl;
    clock_t t1,t2;
    t1=clock();
    
    std::string Text_files_path = "../../../../../Dropbox/Research/Epidemiology_2020/Text_files/Kitas_Schools/";
    
    // Parameters
    int n_ensemble = 2000; //Ensemble size
    int N = 20; //Number of kids
    int T = 4*7; //Total number of days for the simulation
    int t_inc = argv[3][0]-'0'; //Incubation period in days
    int t_inf = 6; //Infectious period in days
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
    
    //Vectors with the counter of a incubations days and infectious days.
    std::vector < int > incubation;
    incubation.resize(N);
    std::vector < int > infectious;
    infectious.resize(N);
    for(int n = 0 ; n<N ; n++){
        incubation[n] = t_inc;
        infectious[n] = t_inf;
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
    
    std::cout << "beta is " << beta << " ,p_in is " << p_in_inf << " and t_inc is " << t_inc << std::endl << std::endl;
    
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
    
    double inf;
    double r1;
    std::vector<int>::iterator a;
    //kids[10] = 1;
    
    //Output file
    std::ofstream fout (Text_files_path+"statistics_days-"+std::to_string(n_testing_days)+"-"+argv[2]+"_beta-"+std::to_string(beta)+"_pin-"+std::to_string(p_in_inf)+"_tau-"+std::to_string(t_inc)+".txt");
    
    for(int j = 0 ; j < n_ensemble ; j++){
        
        //------------ initialize vectors-------------
        for(int n = 0 ; n<N ; n++){
            kids[n] = 0;
            incubation[n] = t_inc;
            infectious[n] = t_inf;
        }
        H = N;
        Inc = 0;
        Inf = 0;
        Rec = 0;
        Det = 0;
        total_inf = 0;
        //--------------------------------------------
        
        // for-loop running over T days
        for(int t = 0; t<=T ; t++){
            int d = t%7;
            update_state_kids(N, &H, &Inc, &Inf, &Rec, &Det, &kids);
            inf = Inf/double(N);
            std::string day = days[d];
            
            if (d<5) { //Week days
                // for-loop running over N kids
                for(int n = 0; n<N ; n++){
                    if(kids[n]==0){ //Kid is healthy
                        r1 = randX(0,1);
                        if(r1<(inf*beta)){
                            kids[n] = 1;
                        }
                    }else if(kids[n]==1){ //Kid is in incubation period
                        if (incubation[n]==0) { //Incubation period is over?
                            kids[n]=2;
                            total_inf++;
                        }else{
                            incubation[n]--;
                        }
                    }else if(kids[n]==2){ //Kid is in infectious period
                        if (infectious[n]==0) { //Infectious period is over?
                            kids[n]=3;
                        } else {
                            infectious[n]--;
                            a = std::find(testing_days.begin(), testing_days.end(), d);
                            if (a != testing_days.end()) { //Is it a testing day?
                                kids[n]=4;
                            } else {
                                
                            }
                        }

                    }
                    else{
                        
                    }
                }
            } else { //Weekends
                // for-loop running over N kids
                for(int n = 0; n<N ; n++){
                    if(kids[n]==0){ //Kid is healthy
                        r1 = randX(0,1);
                        if(r1<(p_in_inf)){
                            kids[n] = 1;
                        }
                    }else if(kids[n]==1){ //Kid is in incubation period
                        if (incubation[n]==0) { //Incubation period is over?
                            kids[n]=2;
                            total_inf++;
                        }else{
                            incubation[n]--;
                        }
                    }else if(kids[n]==2){ //Kid is in infectious period
                        if (infectious[n]==0) { //Infectious period is over?
                            kids[n]=3;
                        } else {
                            infectious[n]--;
                        }
                    }
                    else{
                        
                    }
                }
            }
          
        }
        
        fout << total_inf << std::endl;
    }
    fout.close();
    //------------------------------------------------------------------------------
    std::cout<< ">Simulation completed…"<< std::endl;
    t2= clock();
    std::cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< std::endl;
    return 0;
}
