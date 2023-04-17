------ README FOR THE FIRST IMPLEMENTATION EXERCISE ----------

  ###### PART ALGORITHM
  To compile, go to the flowshopWCT directory and type 'make' in the terminal to compile everything.

  Then, the user can execute the program in different modes, at his convenience. 
  He simply has to write in the command line
  ./flowshopWCT MODE ALGORITHM INITSOL NEIGHBORHOOD PIVOTINGRULE

  where the different values for those variables are (choose one for each):
  - MODE: {individual, all}
  - ALGORITHM: {ii, vnd}
  - INITSOL: {randinit, srz}
  - NEIGHBORHOOD: {transpose, exchange, insert}
  - PIVOTINGRULE: {first, best}

The program will create results txt files written in the CSV format in the folder 'resultsfiles'
/!\ If MODE=individual :
  -  The program computes the results of the PFSP with an individual algorithm (defined by 
  the different values of ALGORITHM, INITSOL, NEIGHBORHOOD, PIVOTINGRULE passed as arguments).
  It does so ON EVERY INSTANCES.

/!\ If MODE=all :
  -  Then it is important to insert bogus values for ALGORITHM, INITSOL, NEIGHBORHOOD, PIVOTINGRULE
  even though those values won't be taken into account since the program will compute ALL the algorithms 
  on ALL the instances
  - Expect more or less 30min of computation time

An exemple of correct input in the command line is:
./flowshopWCT individual ii srz transpose first

/!\ The txt file 'all_performances.txt' inside the folder 'resultsfiles' is not generated automatically.
   I manually gathered all the results of the PREFIX_performance.txt (where PREFIX is the name 
   of the algorithm, e.g. 'ii_srz_transpose_first') to get a better overview of the overall performances.

###### PART STATISTICS

Compile and execute the R file present in resultsfiles. It is NOT in the 'src' folder in order to avoid path problems.
The R file itself does the Wilcoxon test given 2 input sets.
It outputs the p value and whether the two input sets are statistically different or not.

/!\ The input files must be changed manually INSIDE the R file. To see the name of the input files, see the name
    of the files inside the 'resultsfiles' folder.

/!\ Be sure to set the working directory to the source file location ! 

  