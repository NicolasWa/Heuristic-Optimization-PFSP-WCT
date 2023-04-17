--------- README FOR THE SECOND IMPLEMENTATION EXERCISE ----------

  ###### PART ALGORITHM
  To compile, go to the flowshopWCT directory and type 'make' in the terminal to compile everything.

  Then, the user can execute the the algorithm of his choice:
      - Either a Random Iterative Improvement: command ---> ./flowshopWCT individual RII srz insert first
      - Or an Iterative Local Search: 	     command --> ./flowshopWCT individual RII srz insert first
  ./flowshopWCT MODE ALGORITHM INITSOL NEIGHBORHOOD PIVOTINGRULE

  where the different values for those variables are (choose one for each):
  - RII stands for Random Iterative Improvement
  - srz stands for Simplified RZ Heuristic
  - insert stands for the neighborhood "insert"
  - first stands for first improvement

The program will create results txt files written in the CSV format in the folder 'resultsfiles'
/!\ If MODE=individual :



/!\ The txt file 'all_performances.txt' inside the folder 'resultsfiles' is not generated automatically.
   I manually gathered all the results of the PREFIX_performance.txt (where PREFIX is the name 
   of the algorithm, e.g. 'RII_srz_insert_first') to get a better overview of the overall performances.

###### PART STATISTICS

Compile and execute the R file. "statistical_tests_imp2.R" present in resultsfiles. It is NOT in the 'src' folder in order to avoid path problems.

/!\ Be sure to set the working directory to the source file location ! 

  