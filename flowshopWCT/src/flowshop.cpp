/***************************************************************************
 HEURISTIC OPTIMIZATION : Implementation excercise 1
 
 Surname : Wallemacq
 First name: Nicolas
 Student number: 000445394
 ***************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include "pfspinstance.h"
#include <string>
#include <cstring>
#include <fstream>
#include <math.h>

using namespace std;
using namespace std::chrono;

/* #################### PART USEFUL FUNCTIONS TO DEBUG: BEGIN #################### */

void printVector(vector< int > & solution, PfspInstance &instance) {
  int i;
  for (i = 1; i <= solution.size()-1; ++i) {
  //for (i = 1; i <= instance.getNbJob(); ++i) {
    cout << solution[i] << " " ;
  }
  cout << endl;
}

void countCheckSol(vector< int > & solution, PfspInstance &instance) {
  /* This function was useful to assess that there were no duplicates in the solution */
  int i;
  long int sum=0;
  long int sumtrue=0;
  for (i = 1; i <= instance.getNbJob(); ++i) {
    sum+=solution[i];
  }
  for (i = 1; i <= instance.getNbJob(); ++i) {
    sumtrue+=i;
  }
  if (sum==sumtrue) {
    cout << "NO DUPLICATES" << endl;
  }
  else {
    cout << "ERROR --> PRESENCE OF DUPLICATES" << endl;
  }
}

void getBestSolutionsVector(vector <double> &BestSolutionVector){
  ifstream fileInBestSolutions;
  fileInBestSolutions.open("../instances/bestSolutions.txt");
  std::string str; //will be used to read unwanted values (for ex the ',')
  double valueBestSol;
  /* We read the values stored in the bestSolutions txt file and we put it in a vector */
    if ( fileInBestSolutions.is_open() ) {
        //cout << "File  is now open, start to read..." << std::endl;
        for (int i=1; i<=3; i++) { fileInBestSolutions >> str;} // names of the rows
        for (int i=1;i<=60;i++){
            fileInBestSolutions >> str;
            fileInBestSolutions >> str;
            fileInBestSolutions >> valueBestSol;
            BestSolutionVector.push_back(valueBestSol);
        }
    }
    else {cout    << "error while opening file " << std::endl;}
    fileInBestSolutions.close();
}


int verifyArguments(int argc, string mode, string algorithm, string solInit, string neighborhood, string pivotingRule){
  /*  This function verifies that the arguments passed in the command line were correct 
  Possible values for mode: {individual, all}, for algorithm: {ii, vnd}, for solInit: {randinit, srz}, for neighborhood: {transpose, exchange, insert},
  for pivotingRule: {first, best} */
  cout << "RENTRE DANS VERIFY ARGUMENTS" << endl;

  /* Mode "individual" --> we take into account all the arguments passed in the command line*/ 
  if (mode =="individual"){
    if (algorithm=="ii" || algorithm=="ILS"|| algorithm=="RII"|| algorithm=="SA"){
      std::cout << "rentreeeeee" << std::endl;
      if (argc == 1) {cout << "Usage: ./flowshopWCT <instance_file>" << endl; return 0;}
      if (!((solInit=="randinit")==1 || (solInit=="srz")==1)) { 
        cout << "Argument given for initial solution not recognized" << endl;
        return 0;
      }
      if (!((neighborhood=="transpose")==1 || (neighborhood=="exchange")==1 || (neighborhood=="insert")==1)) { 
        cout << "Argument given for neighborhood not recognized" << endl;
        return 0;
      }
      if (!((pivotingRule=="first")==1 || (pivotingRule=="best")==1)) { 
        cout << "Argument given for pivoting rule not recognized" << endl;
        return 0;
      }
    }
    else if (algorithm=="vnd"){
      if (argc == 1) {cout << "Usage: ./flowshopWCT <instance_file>" << endl; return 0;}
      /*if (!((solInit=="srz")==1)) { 
        cout << "Argument given for initial solution must be srz" << endl;
        return 0;
      }
      */
      if (!((neighborhood=="tr_ex_in")==1 || (neighborhood=="tr_in_ex")==1)) {
        cout << "Argument given for neighborhood not recognized" << endl;
        return 0;
      }
      
      if (!((pivotingRule=="first")==1)) { 
        cout << "Argument given for initial solution must be first" << endl;
        return 0;
      }
    }
    else {cout << "Argument given for algorithm not recognized" << endl; return 0;}
  }
  /* Mode "all" --> we neglect all the other arguments passed in the command line since we're going to compute all the algorithms on all the instances*/ 
  else if (mode == "all"){
    cout << "All arguments are accepted since they won't be taken into account except 'all' " << endl;
  }
  else { cout << "Argument given for mode not recognized" << endl; return 0;}

  cout << "Arguments given: " << endl;
  cout << "Mode: " << mode << endl;
  cout << "Algorithm: " << algorithm << endl;
  cout << "Initial solution: " << solInit << endl;
  cout << "Neighborhood:     " << neighborhood << endl;
  cout << "Pivoting rule:    " << pivotingRule << endl;
  return 1;
  }
void PrintResults(long int initialWCT, long int finalWCT,vector< int > &initialSol, vector< int > &solution, PfspInstance instance) {
  /* This function was useful to print the results (initial solution, solution) and verifies that there are no duplicates. 
  This was useful to debug but is not used in the final implementation */
  cout << "Total weighted completion time (initial): " << initialWCT <<endl;
  cout << "Total weighted completion time (final): " << finalWCT << endl;
  cout << "Gain in WCT: " << initialWCT-finalWCT << endl;
  printVector(initialSol, instance);
  printVector(solution, instance);
  countCheckSol(solution, instance);
}

/* #################### PART USEFUL FUNCTIONS TO VISUALIZE: END #################### */


/* #################### PART INITIAL SOLUTION: BEGIN #################### */

int generateRndPosition(int min, int max)
/* The following functions help generating a random initial solution */ 
{
  return ( rand() % max + min );
}

void randomPermutation(int nbJobs, vector< int > & sol)
/* Fill the solution with numbers between 1 and nbJobs, shuffled in order to get a random initial solution */
{
  vector<bool> alreadyTaken(nbJobs+1, false); // nbJobs elements with value false
  vector<int > choosenNumber(nbJobs+1, 0);

  int nbj;
  int rnd, i, j, nbFalse;
  nbj = 0;
  for (i = nbJobs; i >= 1; --i)
  {
    rnd = generateRndPosition(1, i);
    nbFalse = 0;
    /* find the rndth cell with value = false : */

    for (j = 1; nbFalse < rnd; ++j)
      if ( ! alreadyTaken[j] )
        ++nbFalse;
    --j;

    sol[j] = i;

    ++nbj;
    choosenNumber[nbj] = j;
    alreadyTaken[j] = true;
  }
}

void randinit(vector< int > & solution, PfspInstance &instance) {
  /* Initializes the solution with a random solution*/
  randomPermutation(instance.getNbJob(), solution);
}

void srzinit(vector< int > & solution, PfspInstance &instance){
  /* Generates an initial solution using the Simplified RZ heuristic */
  long int WCTseq;
  long int WCTcandidate;
  int temp;
  /* Compute the starting sequence */
  //<int>std::vector start(instance.getNbJob()+1);
  std::vector<int> start(instance.getNbJob()+1);
  start = instance.computeInitVectorHeuristics();

  /* Compute the sequence of jobs*/
  std::vector<int> sequence(1); //first value of this vector is of no importance, I want to make it match the other vectors size used for the solutions
  sequence.push_back(start[1]);
  WCTseq = instance.computeWCT(sequence);
  /* for all the remaining jobs: add one job in sequence and find its position so that it minimizes the WCT */ 
  for (int j=2; j<=instance.getNbJob(); j++) {
    sequence.push_back(start[j]);
    temp = start[j];
    std::vector<int> sequenceCandidate(sequence);
    std::vector<int> bestCandidate(sequence);
    WCTseq = instance.computeWCT(sequence);
    for (int i=1; i<sequenceCandidate.size()-1; i++) { 
      sequenceCandidate.erase(sequenceCandidate.begin()+j);
      sequenceCandidate.insert(sequenceCandidate.begin()+i, temp);
      WCTcandidate = instance.computeWCT(sequenceCandidate);
      if ((WCTcandidate < WCTseq)==1){
        WCTseq = WCTcandidate;
        copy(sequenceCandidate.begin(), sequenceCandidate.end(), bestCandidate.begin());
      }
      copy(sequence.begin(), sequence.end(), sequenceCandidate.begin());
    }
    copy(bestCandidate.begin(), bestCandidate.end(), sequence.begin());
  }
  /* The last value of sequence is the heuristic initial solution. We inject it to solution */
  copy(sequence.begin(), sequence.end(), solution.begin()); 
}

void GenerateInitialSol(string solInit, vector< int > & solution, vector< int > & initialSol, long int &initialWCT, PfspInstance &instance) {
  /* Depending on the solInit chosen, calls the appropriate function to give an initial solution */
  if ((solInit=="randinit")==1) {
    randinit(solution, instance);
  }
  else if ((solInit=="srz")==1){
    srzinit(solution, instance);
  }
  initialWCT = instance.computeWCT(solution);
  copy(solution.begin(), solution.end(), initialSol.begin());
}

/* #################### PART INITIAL SOLUTION: END #################### */


/* #################### PART NEIGHBORHOOD: BEGIN #################### */
int Transpose(vector< int > & solution, vector< int > & solutionNeighbor,string pivotingRule, PfspInstance &instance) {
  int i;
  copy(solution.begin(), solution.end(), solutionNeighbor.begin());
  std::vector< int > best(solution);
  /* --------- First improvement -------*/
  if ((pivotingRule=="first")==1) {
    for (i = 2; i <= instance.getNbJob(); ++i) {
      solutionNeighbor[i]=solution[i-1]; 
      solutionNeighbor[i-1]=solution[i]; 
      if (instance.computeWCT(solutionNeighbor) < instance.computeWCT(solution)) {return 0;} //a first improvement was found, exits the function
      copy(solution.begin(), solution.end(), solutionNeighbor.begin()); //no improvement, we suppress the modifications made to solutionNeighbor
    }
  }
  /* --------- Best improvement -------*/
  else if ((pivotingRule=="best")==1) {
    for (i = 2; i <= instance.getNbJob(); ++i) {
      solutionNeighbor[i]=solution[i-1];
      solutionNeighbor[i-1]=solution[i];
      if (instance.computeWCT(solutionNeighbor) < instance.computeWCT(best)) {
        copy(solutionNeighbor.begin(), solutionNeighbor.end(), best.begin()); //an improvement was found, we keep it in memory in other vector 'best'
      }
      copy(solution.begin(), solution.end(), solutionNeighbor.begin()); //we suppress the modifications made to solutionNeighbor
    }
    /* After the for loop, the vector best contains the best solution. We inject it to solutionNeighbor, not to solution !
    This trick is important for the function FindLocalOptimum() to behave well */
    copy(best.begin(), best.end(), solutionNeighbor.begin()); 
  }
  return 0;
}

int Exchange(vector< int > & solution, vector< int > & solutionNeighbor,string pivotingRule, PfspInstance &instance) {
  int i;
  int j;
  copy(solution.begin(), solution.end(), solutionNeighbor.begin());
  std::vector< int > best(solution);
  /* --------- First improvement -------*/
  if ((pivotingRule=="first")==1) {
    for (i = 1; i <= instance.getNbJob(); ++i) {
      for (j=i+1; j <= instance.getNbJob(); ++j) {
          swap(solutionNeighbor[i], solutionNeighbor[j]);
          if (instance.computeWCT(solutionNeighbor) < instance.computeWCT(solution)) {return 0;} //a first improvement was found, exits the function
          copy(solution.begin(), solution.end(), solutionNeighbor.begin()); // there was no improvement
      }
    }
  } 

  /* --------- Best improvement -------*/
  else if ((pivotingRule=="best")==1) {
    for (i = 1; i <= instance.getNbJob(); ++i) {
      for (j=i+1; j <= instance.getNbJob(); ++j) {
        swap(solutionNeighbor[i], solutionNeighbor[j]);
        if (instance.computeWCT(solutionNeighbor) < instance.computeWCT(best)) {
          copy(solutionNeighbor.begin(), solutionNeighbor.end(), best.begin()); //an improvement was found, we keep it in memory in other vector 'best'
        } 
        copy(solution.begin(), solution.end(), solutionNeighbor.begin()); //we suppress the modifications made to solutionNeighbor
      }
    }
    /* After the for loop, the vector best contains the best solution. We inject it to solutionNeighbor, not to solution !
    This trick is important for the function FindLocalOptimum() to behave well */
    copy(best.begin(), best.end(), solutionNeighbor.begin());
  }
  return 0;
}

int Insert(vector< int > & solution, vector< int > & solutionNeighbor,string pivotingRule, PfspInstance &instance) {
  int i;
  int j;
  int temp;
  copy(solution.begin(), solution.end(), solutionNeighbor.begin());
  std::vector< int > best(solution);
  /* --------- First improvement -------*/
  if ((pivotingRule=="first")==1) {
    for (i = 1; i <= instance.getNbJob(); ++i) {
      for (j=1; j <= instance.getNbJob(); ++j) {
        temp = solutionNeighbor[i];
        solutionNeighbor.erase(solutionNeighbor.begin()+i);
        solutionNeighbor.insert(solutionNeighbor.begin()+j, temp);
        if (instance.computeWCT(solutionNeighbor) < instance.computeWCT(solution)) {return 0;} //a first improvement was found, exits the function
        copy(solution.begin(), solution.end(), solutionNeighbor.begin()); // there was no improvement
      }
    }
  } 
  /* --------- Best improvement -------*/
  else if ((pivotingRule=="best")==1) {
    for (i = 1; i <= instance.getNbJob(); ++i) {
      for (j=1; j <= instance.getNbJob(); ++j) {
        temp = solutionNeighbor[i];
        solutionNeighbor.erase(solutionNeighbor.begin()+i);
        solutionNeighbor.insert(solutionNeighbor.begin()+j, temp);
        if (instance.computeWCT(solutionNeighbor) < instance.computeWCT(best)) {
          copy(solutionNeighbor.begin(), solutionNeighbor.end(), best.begin()); //an improvement was found, we keep it in memory in other vector 'best'
        } 
        copy(solution.begin(), solution.end(), solutionNeighbor.begin()); //we suppress the modifications made to solutionNeighbor
      }
    }
    /* After the for loop, the vector best contains the best solution. We inject it to solutionNeighbor, not to solution !
    This trick is important for the function FindLocalOptimum() to behave well */
    copy(best.begin(), best.end(), solutionNeighbor.begin());
  }
  return 0;
}

/* Performs a Random Walk on solution. It randomly chooses (with equal probability) in which neighborhood to perform the walk */
std::vector <int> randomWalk(vector< int > solution, PfspInstance &instance){
  std::vector <int> solutionNeighbor(solution);
  int neighborhood = rand()%3 + 1; // Randomly chooses a number between 1 and 3
  copy(solution.begin(), solution.end(), solutionNeighbor.begin());
  /* neighborhood=1 --> randomly applies a transposition */
  if (neighborhood==1) {
    int ind1 = (rand()%(instance.getNbJob()-1) +2); // math trick to be sure to get ind1 comprised between [2;nbJobs]
    solutionNeighbor[ind1]=solution[ind1-1]; 
    solutionNeighbor[ind1-1]=solution[ind1];
  }
  /* neighborhood=2 --> randomly applies an exchange */
  if (neighborhood==2) {
    //cout << "random exchange" << endl;
    int ind1 = (rand()%instance.getNbJob() ) +1;
    int ind2 = (rand()%instance.getNbJob() ) +1;
    swap(solutionNeighbor[ind1], solutionNeighbor[ind2]);
  }
  /* neighborhood=3 --> randomly applies an insertion */
  if (neighborhood==3) {
    int ind1 = (rand()%instance.getNbJob() ) +1;
    int ind2 = (rand()%instance.getNbJob() ) +1;
    int temp = solutionNeighbor[ind1];
    solutionNeighbor.erase(solutionNeighbor.begin()+ind1);
    solutionNeighbor.insert(solutionNeighbor.begin()+ind2, temp);
  }
  return solutionNeighbor;
}

/* #################### PART NEIGHBORHOOD: END #################### */

/* #################### PART ALGORITHMS: BEGIN #################### */

void CreateNameFile(string &nameFileResults, string algorithm, string solInit, string neighborhood, string pivotingRule){
  nameFileResults.append("./resultsfiles/");
  nameFileResults.append(algorithm);
  nameFileResults.append("_");
  nameFileResults.append(solInit);
  nameFileResults.append("_");
  nameFileResults.append(neighborhood);
  nameFileResults.append("_");
  nameFileResults.append(pivotingRule);
  nameFileResults.append(".txt");
}

void IfThenElseVND(int &i, long int &WCT, long int &WCTNeighbor, PfspInstance &instance, vector< int > & solution, vector< int > & solutionNeighbor){
  WCT = instance.computeWCT(solution);
  WCTNeighbor = instance.computeWCT(solutionNeighbor);
  copy(solutionNeighbor.begin(), solutionNeighbor.end(), solution.begin());
  if (WCTNeighbor < WCT) {
    copy(solutionNeighbor.begin(), solutionNeighbor.end(), solution.begin());
    i=1;
  }
  else {i=i+1;}
}

/* Metropolis condition. Returns of probability of acceptance */
float ProbaAcceptMetropolis(float T, vector< int > &solution, vector< int > & solutionNeighbor, PfspInstance &instance){
  float fs = instance.computeWCT(solution);
  float fsprime = instance.computeWCT(solutionNeighbor);
  if (fsprime <= fs){
    return 1;
  }
  else {
    return ( exp((fs-fsprime)/T) );
  }
}

/* Returns 1 if the new proposed solution is accepted. Returns 0 if it is rejected*/
int Accept(vector< int > &solution, vector< int > &solutionNeighbor, PfspInstance &instance, float T){
  float pAccept = ProbaAcceptMetropolis(T, solution, solutionNeighbor, instance); //Metropolis condition for the probability of accepting
  int proba = rand()%100 + 1;
  if (proba <= pAccept*100){ //we accept the worsening solution
    return 1;
  }
  return 0;
}

/* The Perturbation implemented performs nbPerturbations times a Random Walk. This Random Walk randomly chooses a neighbor
with either the neighborhood insert, transpose or exchange */
std::vector <int> Perturbation(vector< int > solution, PfspInstance &instance, int nbPerturbations){
  std::vector <int> solutionNeighbor(solution);
  for (int i=1; i<=nbPerturbations; i++){solutionNeighbor=randomWalk(solution, instance);}
  return solutionNeighbor;
}

std::vector <int> localSearch(vector< int > solution, string neighborhood, string pivotingRule, PfspInstance &instance){
  std::vector <int> solutionNeighbor(solution);
  long int WCTNeighbor = 0; // to enter in the while the first time
  long int WCT = instance.computeWCT(solution);
  while (WCTNeighbor < WCT) {
    if ((neighborhood=="transpose")==1){ Transpose(solution, solutionNeighbor, pivotingRule, instance);}
    else if ((neighborhood=="exchange")==1){ Exchange(solution, solutionNeighbor, pivotingRule, instance);}
    else if ((neighborhood=="insert")==1){ Insert(solution, solutionNeighbor, pivotingRule, instance);}
    WCT = instance.computeWCT(solution); 
    WCTNeighbor = instance.computeWCT(solutionNeighbor); 
    copy(solutionNeighbor.begin(), solutionNeighbor.end(), solution.begin()); 
  }
  return(solution);
}


void FindLocalOptimum(string algorithm, string neighborhood, string pivotingRule, PfspInstance &instance, vector< int > &solution, vector< int > &solutionNeighbor, long int &finalWCT){
    /* Finds a local optimum, starting from an initial solution.
    The algorithm chosen depends on the variables algorithm, neighborhood, pivotingRule that were passed as arguments in the command line */
    long int WCTNeighbor = 0;
    long int WCT = instance.computeWCT(solution);
    long int WCTprev;

    /* Local Search in a single neighborhood */
    if (algorithm=="ii") {
      solution=localSearch(solution,neighborhood, pivotingRule, instance);
    }
    else if (algorithm=="ILS"){
      // the vector solution already contains the heuristic initial solution
      std::vector <int> Sol_init(instance.getNbJob()+1); 
      copy(solution.begin(), solution.end(),Sol_init.begin()); // Sol_init = SRZheuristic(instance)
      
      std::vector <int> Sol(instance.getNbJob()+1); 
      Sol = localSearch(Sol_init,neighborhood, pivotingRule, instance); // Sol = LocalSearch(Sol_init)

      std::vector <int> Sol_best(instance.getNbJob()+1); 
      copy(Sol.begin(), Sol.end(),Sol_best.begin()); //Sol_best = Sol;

      std::vector <int> Sol_prev(instance.getNbJob()+1);
      std::vector <int> Sol_pert(instance.getNbJob()+1);
      std::vector <int> Sol_current(instance.getNbJob()+1);

      /* Initialization of the ILS parameters */
      float T = 700; // Temperature (for the Metropolis condition)
      int nbPerturbationStart = (instance.getNbJob())/50;
      int adaptativeNbPerturbations= nbPerturbationStart;
      int countStuckLocalMin = 0;
      double stop_criterion;

      if (instance.getNbJob()==50){ stop_criterion = 500*0.07726; } //taken from my results of imp1 
      else if (instance.getNbJob()==100){ stop_criterion = 350;} // max (cf tasks)
      high_resolution_clock::time_point t_start = high_resolution_clock::now(); // start the chrono
      bool continue_search = true;

      /* while stop criterion is not met */
      while(continue_search){
        copy(Sol.begin(), Sol.end(), Sol_prev.begin()); // Sol_prev = Sol

        Sol_pert = Perturbation(Sol, instance, adaptativeNbPerturbations); // Sol_pert = perturbation(Sol)

        Sol_current = localSearch(Sol_pert, neighborhood, pivotingRule, instance); // Sol_current = LocalSearch(Sol_pert)

        if (instance.computeWCT(Sol_current) == instance.computeWCT(Sol_prev)){
          countStuckLocalMin+=1; 
          if (countStuckLocalMin>=4){ // We need t augment the number of perturbations
            adaptativeNbPerturbations+=1;
            countStuckLocalMin=0;
          }
        }
        else{
            if (instance.computeWCT(Sol_current) < instance.computeWCT(Sol_best)){
              copy(Sol_current.begin(), Sol_current.end(), Sol_best.begin()); // Sol_best = Sol_current
            }
        }
        /* Metropolis acceptance criterion */
        if (Accept(Sol, Sol_current, instance,T)==1){
          copy(Sol_current.begin(), Sol_current.end(), Sol.begin()); // solution accepted
          //std::cout << "Metropis ACCEPTED << std::endl;
        }
        else{
          //std::cout << "Metropolis REJECTED" << std::endl;
        }

        if (instance.computeWCT(Sol) != instance.computeWCT(Sol_prev)){
          /* We got a different solution than the previous round so we seem to have escaped the local minima */
          countStuckLocalMin = 0;
          adaptativeNbPerturbations = nbPerturbationStart;
        }
        /* Verifying stop criterion*/
        high_resolution_clock::time_point t_now = high_resolution_clock::now(); // end the chrono
        duration<double> elapsed_time = duration_cast<duration<double> >(t_now - t_start);
        if (elapsed_time.count() >= stop_criterion) {continue_search=false;}
      }
      copy(Sol_best.begin(), Sol_best.end(), solution.begin()); // final solution
    }

    else if (algorithm=="RII"){
      std::vector <int> Sol_prev(instance.getNbJob()+1);
      std::vector <int> Sol_best(instance.getNbJob()+1); 
      copy(solution.begin(), solution.end(),Sol_best.begin()); //Sol_best = Sol;

      /* Initializng the parameters of the RII */
      int wp=2;
      int p;
      int nbPerturbations = (instance.getNbJob())/25;
      double stop_criterion;
      if (instance.getNbJob()==50){ stop_criterion = 500*0.07726; } //taken from my results of imp1 
      else if (instance.getNbJob()==100){ stop_criterion = 350;} // max (cf tasks)
      high_resolution_clock::time_point t_start = high_resolution_clock::now(); // start the chrono
      bool continue_search = true;

      while(continue_search){
        copy(solution.begin(), solution.end(),Sol_prev.begin());
        int p = rand()%100 + 1;

        if (p<=wp){ // We perturb the solution
          solution = Perturbation(solution, instance, nbPerturbations);
        }
        else{
          Insert(solution, solutionNeighbor, pivotingRule, instance);
          copy(solutionNeighbor.begin(), solutionNeighbor.end(), solution.begin());
        }
        if (instance.computeWCT(solution) == instance.computeWCT(Sol_prev)){
          solution = Perturbation(solution, instance, nbPerturbations);
        }
        if (instance.computeWCT(solution) < instance.computeWCT(Sol_best)){
          copy(solution.begin(), solution.end(), Sol_best.begin()); // Sol_best = Sol_current
        }

        /* Verifying stop criterion*/
        high_resolution_clock::time_point t_now = high_resolution_clock::now(); // end the chrono
        duration<double> elapsed_time = duration_cast<duration<double> >(t_now - t_start);
        if (elapsed_time.count() >= stop_criterion) {continue_search=false;}
      }
      copy(Sol_best.begin(), Sol_best.end(), solution.begin()); // final solution
    }

    /* Variable Neighborhood Descent */
    else if (algorithm=="vnd"){
      int i=1;
      while (i<=3){
        if (neighborhood=="tr_ex_in"){
          if (i==1) {
            Transpose(solution, solutionNeighbor, pivotingRule, instance);
            IfThenElseVND(i, WCT, WCTNeighbor, instance, solution, solutionNeighbor);
          }
          if (i==2) {
            Exchange(solution, solutionNeighbor, pivotingRule, instance);
            IfThenElseVND(i, WCT, WCTNeighbor, instance, solution, solutionNeighbor);
          }
          if (i==3) {
            Insert(solution, solutionNeighbor, pivotingRule, instance);
            IfThenElseVND(i, WCT, WCTNeighbor, instance, solution, solutionNeighbor);
          }
        }
        else if (neighborhood=="tr_in_ex"){
          if (i==1) {
            Transpose(solution, solutionNeighbor, pivotingRule, instance);
            IfThenElseVND(i, WCT, WCTNeighbor, instance, solution, solutionNeighbor);
          }
          if (i==2) {
            Insert(solution, solutionNeighbor, pivotingRule, instance);
            IfThenElseVND(i, WCT, WCTNeighbor, instance, solution, solutionNeighbor);
          }
          if (i==3) {
            Exchange(solution, solutionNeighbor, pivotingRule, instance);
            IfThenElseVND(i, WCT, WCTNeighbor, instance, solution, solutionNeighbor);
          }
        }
      }
    }

  /* We end up with a final soluton, it is possible to compute the final Weighted sum Completion Time*/
  finalWCT = instance.computeWCT(solution);  
  }

int PFSP(string algorithm, string solInit, string neighborhood, string pivotingRule){
  cout << " -------------------------------------------------------- PFSP begins on one instance" << endl; //to spot the beggining of an instance in the terminal
  std::string nameFileResults;
  std::string onetothirty[30] = { "01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20", "21", "22","23","24","25","26","27","28","29","30"};
  std::string str("../instances/100_20_");
  int seeds[] = {1997,1998,1999,2000,2001,2002};
  int positionToInsert;
  /* Opening/Creating a txt file to later print the solutions inside it*/
  CreateNameFile(nameFileResults, algorithm, solInit, neighborhood, pivotingRule);
  ofstream myfile;
  myfile.open (nameFileResults);
  myfile << "Instance" << " , " << "WCT_1" << " , "<< "WCT_2" << " , " << "WCT_3" << " , "<< "WCT_4" << " , "<< "WCT_5" << " , "<< "Avg_RPD" <<'\n'; //names of the columns

  for (int f=1; f<=2; f++){ // f=1 is for the 100 jobs files and f=2 is for the 50 jobs files
    if (f==1){positionToInsert=20;}
    else if (f==2){ str="../instances/50_20_"; positionToInsert=19;}
    for (int i=1; i<=30; i++) { //TO DO: PUT IT AGAIN
      /* Naming the instance accordingly */
      str.insert(positionToInsert, onetothirty[i-1]);
      char* fileName = strcpy(new char[str.length() + 1], str.c_str());
      // /* Initializing local variables */ 
      long int initialWCT;
      long int finalWCT;
      float avgRPD;
      high_resolution_clock::time_point t1;
      high_resolution_clock::time_point t2;
      duration<double> time_span;
      std::vector <int> WCT_vec(5);
      std::vector <float> RPD_vec(5);
      float sumWCT=0;
      duration<double>  sumTime;
      std::vector <double> Time_vec(5);
      std::vector <double> BestSolutionVector;
      getBestSolutionsVector(BestSolutionVector);

      /* Creating instance object */
      PfspInstance instance;
      /* Read data from file */
      if (! instance.readDataFromFile(fileName))
        return 1;

      /* Create vectors of int to represent the initial solution, the solution given by the neighbor, and the (final) solution
        WARNING: By convention, we store the jobs starting from index 1, thus the size nbJob + 1. */
      std::vector <int> solution ( instance.getNbJob()+ 1 );
      std::vector <int> solutionNeighbor ( instance.getNbJob()+ 1 );
      std::vector <int> initialSol( instance.getNbJob()+ 1 );

      if (algorithm=="SA" || algorithm == "ILS" || algorithm == "RII"){
        avgRPD=0;
        for (int iter=0; iter<5; iter++){ //We run each instance 5 times and then we'll compute the average RPD
          /* Initializing local variables */ 
          initialWCT=0;
          finalWCT=0;
          std::fill(solution.begin(), solution.end(), 0);
          std::fill(solutionNeighbor.begin(), solutionNeighbor.end(), 0);
          std::fill(initialSol.begin(), initialSol.end(), 0);
          srand(seeds[iter]);
          t1 = high_resolution_clock::now(); // start the chrono
          /* Generation of an initial solution */
          GenerateInitialSol(solInit, solution, initialSol, initialWCT, instance);
          /* Finding a local optimum */ 
          FindLocalOptimum(algorithm, neighborhood, pivotingRule, instance, solution, solutionNeighbor, finalWCT); 
          t2 = high_resolution_clock::now(); // end the chrono
          time_span = duration_cast<duration<double> >(t2 - t1);
          /* Printing the results in the terminal */
          cout << "Initial Solution: " ; printVector(initialSol, instance);
          cout << "Final Solution:   " ; printVector(solution, instance);
          cout << "Initial WCT: " << initialWCT << endl;
          cout << "Final WCT:   " << finalWCT << endl;

          WCT_vec[iter]= finalWCT;
          RPD_vec[iter] = (finalWCT-BestSolutionVector[f*30 -30 +i-1])/BestSolutionVector[f*30 -30 +i-1];
          cout << " ---> BestSol=" << BestSolutionVector[f*30 -30 +i -1] << " ---->  RPD[iter] = " << RPD_vec[iter] << endl;
          avgRPD += RPD_vec[iter];
        }
      avgRPD = (avgRPD/5)*100;
      myfile << str << " , " << WCT_vec[0] << " , " << WCT_vec[1] << " , " << WCT_vec[2] << " , " << WCT_vec[3] << " , " << WCT_vec[4] << " , " << avgRPD << '\n';

      }
      else {
          /* Initializing local variables */ 
        initialWCT=0;
        finalWCT=0;
        std::fill(solution.begin(), solution.end(), 0);
        std::fill(solutionNeighbor.begin(), solutionNeighbor.end(), 0);
        std::fill(initialSol.begin(), initialSol.end(), 0);
        srand (1997);
        t1 = high_resolution_clock::now(); // start the chrono
        /* Generation of an initial solution */
        GenerateInitialSol(solInit, solution, initialSol, initialWCT, instance);
        /* Finding a local optimum */ 
        FindLocalOptimum(algorithm, neighborhood, pivotingRule, instance, solution, solutionNeighbor, finalWCT); 
        t2 = high_resolution_clock::now(); // end the chrono
        time_span = duration_cast<duration<double> >(t2 - t1);
        /* Printing the results in the results txt file */
        myfile << str << " , " << finalWCT << " , " << time_span.count() <<'\n';
        /* Printing the results in the terminal */
        cout << "Initial Solution: " ; printVector(initialSol, instance);
        cout << "Final Solution:   " ; printVector(solution, instance);
        cout << "Initial WCT: " << initialWCT << endl;
        cout << "Final WCT:   " << finalWCT << endl;
      }
      
      str.erase(str.begin()+positionToInsert, str.end());
      
    }
  }
  myfile.close();
  cout << " -------------------------------------------------------- PFSP ends on one instance (ran 5 times if algorithm = RII or ILS" << endl; //to spot the end of an instance in the terminal
  return 0;
}

/* #################### PART ALGORITHMS: END #################### */


/* #################### PART PERFORMANCE: BEGIN #################### */
/* The following functions help assess the performances of the algorithms */ 
void computeAvgTime(vector< float > ComputationTimeVector, double &avgTime100jobs, double &avgTime50jobs) {
  /* Computes the average time it takes to compute an algorithm for instances of 1O0 jobs and 50 jobs */
    for (int i=0; i<=29; i++) {
        avgTime100jobs += ComputationTimeVector[i];
        avgTime50jobs += ComputationTimeVector[i+30];
    }
    avgTime100jobs = avgTime100jobs/30;
    avgTime50jobs = avgTime50jobs/30;
    cout << "Average time for 100 jobs:  " << avgTime100jobs << endl;
    cout << "Average time for 50 jobs:  " << avgTime50jobs << endl;
}


void computeAvgPercentageDeviation (vector< double > SolutionVector, vector< double > BestSolutionVector, double &relativePercentageDeviation) {
    /* Computes the average percentage deviation from the best known solution */
    relativePercentageDeviation =0; 
    for (int i=0; i<=59; i++) {
        relativePercentageDeviation = relativePercentageDeviation + (SolutionVector[i]-BestSolutionVector[i])/BestSolutionVector[i] ;
    }
    relativePercentageDeviation= relativePercentageDeviation*100/60;
    cout << "Average percentage deviation from the best known solution:   " << relativePercentageDeviation << endl;
}

void computeAvgRPDInstances(vector< double >RPDVector, double &AvgRPD_100, double &AvgRPD_50){
  AvgRPD_100=0;
  AvgRPD_50=0;
  for (int i=0; i<=29; i++) {
        AvgRPD_100 = AvgRPD_100 + RPDVector[i] ;
  }
  AvgRPD_100 = AvgRPD_100/30;
  for (int i=30; i<=59; i++) {
        AvgRPD_50 = AvgRPD_50 + RPDVector[i] ;
  }
  AvgRPD_50 = AvgRPD_50/30;
}

void performance(string algorithm, string solInit, string neighborhood, string pivotingRule){
  /* Assesses the performance of an algorithm and prints the results in a txt file with the CSV format */
    std::string performanceFile;
    performanceFile.append("./resultsfiles/"); performanceFile.append(algorithm); performanceFile.append("_");
    performanceFile.append(solInit); performanceFile.append("_"); performanceFile.append(neighborhood); performanceFile.append("_");
    performanceFile.append(pivotingRule); performanceFile.append("_performance"); performanceFile.append(".txt");
    ifstream fileInSolutions, fileInBestSolutions; //input files
    ofstream myOutputFile; //output file
    std::string str; //will be used to read unwanted values (for ex the ',')

    /* Creating the name of the output file where we'll print the performance results of the algorithm */
    std::string resultFile;
    CreateNameFile(resultFile,algorithm, solInit, neighborhood, pivotingRule);
    double valueSol, valueBestSol;
    double avgTime50jobs, avgTime100jobs, relativePercentageDeviation;
    float time;
    std::vector <double> SolutionVector, BestSolutionVector;
    std::vector <float> ComputationTimeVector;

    /* Opening the different txt files */
    fileInBestSolutions.open("../instances/bestSolutions.txt");
    fileInSolutions.open(resultFile);
    myOutputFile.open(performanceFile);

  /* We read the values stored in the solution txt file and put we it in vectors */
	if ( fileInSolutions.is_open() ) {
        //cout << "File  is now open, start to read..." << std::endl;
        for (int i=1; i<=5; i++) { fileInSolutions >> str;} // names of the rows
        for (int i=1;i<=60;i++){
            fileInSolutions >> str; //name of the instance, we don't want it
            fileInSolutions >> str; // ',' , we don't want it
            /* We recover the WCT for the algorithm on one instance*/
            fileInSolutions >> valueSol; 
            SolutionVector.push_back(valueSol);
            fileInSolutions >> str;  // ',' , we don't want it
            /* We recover the time it took to compute the solution for the algorithm on one instance */
            fileInSolutions >> time;
            ComputationTimeVector.push_back(time);
        }
    }
    else {cout    << "error while opening file " << std::endl;}
    fileInSolutions.close();
    /* We read the values stored in the bestSolutions txt file and we put it in a vector */
    if ( fileInBestSolutions.is_open() ) {
        //cout << "File  is now open, start to read..." << std::endl;
        for (int i=1; i<=3; i++) { fileInBestSolutions >> str;} // names of the rows
        for (int i=1;i<=60;i++){
            fileInBestSolutions >> str;
            fileInBestSolutions >> str;
            fileInBestSolutions >> valueBestSol;
            BestSolutionVector.push_back(valueBestSol);
        }
    }
    else {cout    << "error while opening file " << std::endl;}
    fileInBestSolutions.close();
    /* Now that we read the values in the files and that we put them in vectors, we can compute the average values */
    computeAvgTime(ComputationTimeVector, avgTime100jobs, avgTime50jobs);
    computeAvgPercentageDeviation(SolutionVector, BestSolutionVector, relativePercentageDeviation);
    /* We write those values in the output file */
    myOutputFile << "Avg_time_100_jobs" << " , " << "Avg_time_50_jobs" << " , " << "RPD" <<'\n';
    myOutputFile << avgTime100jobs << " , " << avgTime50jobs << " , " << relativePercentageDeviation <<'\n';
    
    myOutputFile.close();
}

void performanceSLS(string algorithm, string solInit, string neighborhood, string pivotingRule){
  /* Assesses the performance of an algorithm and prints the results in a txt file with the CSV format.
   Relatively the same function than "performance() but adapted for the results of the SLS methods" */
    std::string performanceFile;
    performanceFile.append("./resultsfiles/"); performanceFile.append(algorithm); performanceFile.append("_");
    performanceFile.append(solInit); performanceFile.append("_"); performanceFile.append(neighborhood); performanceFile.append("_");
    performanceFile.append(pivotingRule); performanceFile.append("_performance"); performanceFile.append(".txt");
    ifstream fileInSolutions, fileInBestSolutions; //input files
    ofstream myOutputFile; //output file
    std::string str; //will be used to read unwanted values (for ex the ',')

    /* Creating the name of the output file where we'll print the performance results of the algorithm */
    std::string resultFile;
    CreateNameFile(resultFile,algorithm, solInit, neighborhood, pivotingRule);
    double valueRPD;

    std::vector <double> RPDVector; 
    fileInSolutions.open(resultFile);
    myOutputFile.open(performanceFile);

  /* We read the avgRPD values stored in the solution txt file and put we it in a vector */
	if ( fileInSolutions.is_open() ) {
        //cout << "File  is now open, start to read..." << std::endl;
        for (int i=1; i<=13; i++) { fileInSolutions >> str;} // names of the rows
        for (int i=1;i<=60;i++){
          for (int i=1; i<=12; i++) { fileInSolutions >> str;}
          /* We recover the (avg) RPD for the algorithm on one instance*/
          fileInSolutions >> valueRPD; 
          RPDVector.push_back(valueRPD);
        }
    }
    else {cout    << "error while opening file " << std::endl;}
    fileInSolutions.close();
    double AvgRPD_100, AvgRPD_50, AvgRPD_Global;
    computeAvgRPDInstances(RPDVector, AvgRPD_100, AvgRPD_50);
    AvgRPD_Global = (AvgRPD_100+AvgRPD_50)/2;
    /* We write those values in the output file */
    myOutputFile << "Avg_RPD_100_jobs" << " , " << "Avg_RPD_50_jobs" << " , " << "Avg_RPD_Global" <<'\n';
    myOutputFile << AvgRPD_100 << " , " << AvgRPD_50 << " , " << AvgRPD_Global <<'\n';
    myOutputFile.close();
}



/* #################### PART PERFORMANCE: END #################### */

/* ####################*/
/* Main function of the program */ 
int main(int argc, char *argv[])
{
  /* Depending the arguments passed in the commandline, it either:
  - mode='individual' --> Computes the results of the PFSP with an individual algorithm ON EVERY INSTANCES. 
  This individual algorithm depends on the other arguments passed in the command line by the user: 
  he chooses the type of algorithm to use (ii or vnd), the type of solInit (randinit or srz), 
  the type of neighborhood (transpose, exchange or insert) and the type of pivotingRule (first or best)

  - mode='all' --> Mode "all" --> It neglects all the other arguments passed in the command line 
  since we're going to compute ALL THE ALGORITHMS ON EVERY INSTANCES
  
  It outputs the results in the command line and creates 2 txt files per algorithm (one for to store the final WCT 
  and one to compute the performances). Those results are stored in the folder "resultsfiles"
  */ 
  std::string mode(argv[1]); //all algorithms or one individual individual algorithms (both on every instances)
  std::string algorithm;
  std::string solInit;
  std::string neighborhood;
  std::string pivotingRule;
  std::string algorithmsVector[2] = {"ii", "vnd"};
  std::string solInitVector[2] = {"randinit", "srz"};
  std::string neigborhoodsVector[3] = {"transpose", "exchange", "insert"};
  std::string neigborhoodsVNDVector[2] = {"tr_ex_in", "tr_in_ex"};
  std::string pivotingRulesVector[2] = {"first", "best"};

  
  if (mode=="individual"){ 
    /* The user specified he want to compute the results on all the instances for ONE specific algorithm 
    so we recovere those specifications in the variables algorithm, solInit, neighborhood, pivotingRule */
    std::string algorithm(argv[2]);
    std::string solInit(argv[3]);
    std::string neighborhood(argv[4]);
    std::string pivotingRule(argv[5]);
    if (algorithm=="ii"){
      if (verifyArguments(argc, mode, algorithm, solInit, neighborhood, pivotingRule)==0){ return 0 ;}
      PFSP(algorithm, solInit, neighborhood, pivotingRule);
      performance(algorithm, solInit, neighborhood, pivotingRule);
    }
    else if (algorithm=="vnd"){
      if (verifyArguments(argc, mode, algorithm, solInit, neighborhood, pivotingRule)==0){ return 0 ;}
      PFSP(algorithm, solInit, neighborhood, pivotingRule);
      performance(algorithm, solInit, neighborhood, pivotingRule);
    }
    else if (algorithm=="ILS"){
      if (verifyArguments(argc, mode, algorithm, solInit, neighborhood, pivotingRule)==0){ return 0 ;}
      PFSP(algorithm, solInit, neighborhood, pivotingRule);
      performanceSLS(algorithm, solInit, neighborhood, pivotingRule);
      
    }
    else if (algorithm=="RII"){
      if (verifyArguments(argc, mode, algorithm, solInit, neighborhood, pivotingRule)==0){ return 0 ;}
      PFSP(algorithm, solInit, neighborhood, pivotingRule);
      performanceSLS(algorithm, solInit, neighborhood, pivotingRule);
    }
    else if (algorithm=="SA"){
      if (verifyArguments(argc, mode, algorithm, solInit, neighborhood, pivotingRule)==0){ return 0 ;}
      PFSP(algorithm, solInit, neighborhood, pivotingRule);
      performanceSLS(algorithm, solInit, neighborhood, pivotingRule);
    }
    
  }
  else if (mode=="all"){
    /* The user specified he wants to compute all the algorithms on all the instances so beside "mode" 
     we do not listen to the other arguments passed in the command line */
    for (int a=0; a<sizeof(algorithmsVector)/sizeof(algorithmsVector[0]); a++){
      algorithm = algorithmsVector[a];
      if (algorithm=="ii"){
        for (int s=0; s<sizeof(solInitVector)/sizeof(solInitVector[0]); s++){
          solInit = solInitVector[s];
          /* We initialize the seed. We work with the same seed in order to be able to compare all the algorithms between each other */
          srand (1997);
          for (int n=0; n<sizeof(neigborhoodsVector)/sizeof(neigborhoodsVector[0]); n++){
            neighborhood = neigborhoodsVector[n];
            for (int p=0; p<sizeof(pivotingRulesVector)/sizeof(pivotingRulesVector[0]); p++){
              pivotingRule = pivotingRulesVector[p];
              PFSP(algorithm, solInit, neighborhood, pivotingRule);
              performance(algorithm, solInit, neighborhood, pivotingRule);
            }
          }
        }
      }

      else if (algorithm=="vnd"){
        solInit = "srz";
        pivotingRule = "first";
        for (int n=0; n<sizeof(neigborhoodsVNDVector)/sizeof(neigborhoodsVNDVector[0]); n++){
            neighborhood = neigborhoodsVNDVector[n];
            PFSP(algorithm, solInit, neighborhood, pivotingRule);
            performance(algorithm, solInit, neighborhood, pivotingRule);
        }
      }
    }
  }
  else {cout << "Argument given for mode not recognized" << endl;} 
  return 0;
}