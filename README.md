# Explaining human cooperation through a dual mechanism of individual and social learning
Our study analyzes a theoretical model using both mathematical and numerical methods.
We provide the code capable of generating all the data presented in this article.
Readers can download the code and run it in Julia.

main last step: a base script with parameter settings (λ_1,…,λ_(μ-1),λ_μ)=(0,…,0,1)
main any lambda: a script that allows arbitrary combinations of (λ_1,…,λ_μ). The runtime is longer than "main last step".
main with calculation of q: a script that calculates the stationary average IBD probability 𝑞 using parameter settings (λ_1,…,λ_(μ-1),λ_μ)=(0,…,0,1). The runtime is substantially longer than "main last step", as it requires solving the system of equations defined in Eq. (10) of the paper.
If the memory required for large 𝑁 exceeds your computer’s capacity, you can use a series of smaller 𝑁 values to fit the stationary average IBD probability 𝑞 for larger 𝑁.


The estimated running time of "main last step" is over 10^5 seconds, which is significantly longer than the calculation for fixed probability. 
This is because the evolutionary process we focus on requires computing the average frequency of cooperators at the steady state. 
Each data point is obtained by averaging over 200 independent runs, and and each run simulates 2000*N steps as transient state and 2000*N steps as steady state. 

