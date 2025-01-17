# Interdependence-between-individual-and-social-learning-shape-cooperation-with-asymmetric-exploration
Our study analyzes a theoretical model using both mathematical and numerical methods. 
Here, we provide a computer simulation program that can be used to obtain all the data involved in our article ("readme_data_analyze" provides how to  obtain the desired data from the program). 
Readers can download the code here, and run it in MATLAB.

main_b can calculate the evolution data of the cooperator ratio and the results of Fig.4e, main_a can calculate all other results.

The estimated running time of main_a is more than 2*10^5 seconds, which is much longer than the calculation of fixed probability. 
This is because the evolutionary process we are focusing on needs to calculate the average frequency of cooperators in the steady state. Each run needs to simulate 2000*N steps as transient state and 2000*N steps as steady state. 
The program can be run more than ten times faster by converting program_c to a mex file.
