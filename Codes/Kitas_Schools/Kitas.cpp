//
//  Kitas.cpp
//  
//
//  Created by Roberto Moran Tovar on 27.02.21.
//  Simulate the dynamics of infections in a kita.

#include "./lib/Kitas.hpp"

#include "./lib/random.cpp"

//using namespace std;

int main(int argc, char* argv[]){
    
    std::cout<<">Running simulation of Kitas ..."<< std::endl;
    clock_t t1,t2;
    t1=clock();
    
    std::string Text_files_path = "../../../../../Dropbox/Research/Epidemiology_2020/Text_files/Kitas_Schools/";
    
    // Parameters
    int N = 200; //Number of kids
    int T = 8*7; //Total number of days for the simulation
    int t_inc = 4; //Incubation period in days
    int t_inf = 6; //Infectious period in days
    double beta = .5; //Infection rate days^-1
    double p_det = 0.8; //Probability of detection.
    double p_in_inf = 0.1; //Probability of influx of a new infection to the kita. Should be proportional to the prevalence in the city.
    
    
    //Vector with the kids of a kita: 0 = Healthy, 1 = Incubation, 2 = Infected, 3 = Recovered, 4 = Detected.
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
        std::cout << d << ":" << days[d] << "\t";
    }
    std::cout << std::endl;
    
    std::vector<int> testing_days;
    int n_testing_days;
    
    std::cout << "How many days do you want to test per week? 1, 2, or 3?" << std::endl;
    std::cin >> n_testing_days;
    
    for (int d = 1; d<=n_testing_days; d++){
        int m;
        std::cout << "What is the day " << d<<"?"<< std::endl;
        std::cin >> m;
        testing_days.push_back(m);
    }
    //--------------------------------
    
    // initial number of kids in each compartment
    int H = N;
    int Inc = 0;
    int Inf = 0;
    int Rec = 0;
    int Det = 0;
    
    double inf;
    double r1;
    std::vector<int>::iterator a;
    kids[10] = 1;
    
    //Output file
    std::ofstream fout (Text_files_path+"output_"+std::to_string(n_testing_days)+"-testing_days.txt");
    
    // for-loop running over T days
    for(int t = 0; t<=T ; t++){
        int d = t%7;
        update_state_kids(N, &H, &Inc, &Inf, &Rec, &Det, &kids);
        fout << H << " " << Inc << " " << Inf << " " << Rec << " " << Det << std::endl;
        inf = Inf/double(N);
        std::string day = days[d];
        
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
    }
    fout.close();
    //------------------------------------------------------------------------------
    std::cout<< ">Simulation completed…"<< std::endl;
    t2= clock();
    std::cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< std::endl;
    return 0;
}
