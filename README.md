# Explaining human cooperation through a dual mechanism of individual and social learning
Our study analyzes a theoretical model using both mathematical and numerical methods. 
Here, we provide a base code with parameter settings (λ_1,…,λ_(μ-1),λ_μ )=(0,…,0,1).
All the data presented in our study can be obtained by appropriately extending this code.
Readers can download the code here, and run it in Julia.

The estimated running time of main_a is more than 2*10^5 seconds, which is much longer than the calculation of fixed probability. 
This is because the evolutionary process we are focusing on needs to calculate the average frequency of cooperators in the steady state. Each run needs to simulate 2000*N steps as transient state and 2000*N steps as steady state. 

