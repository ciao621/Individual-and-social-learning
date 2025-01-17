

tic;

%population size
N=500;
%intensity of selection, w<<1/N
w_all=0.01;
%the duration of transient and steady states.
time_all=2000*N;

  
%degree
degree=4;
%payoff parameter;
b_all=[3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13];
c=1;

%parameter mu can be any distribution. In practice, we only need to generate mu_all according to a specific distribution.
mu_all=500*ones(1,N);

%probability of individual learning
probIL_all=tanh(w_all*[c*degree]);

%the time interval used to calculate neighborhood changes
scale=200;

%number of network structures
number_structural=1;


[neighbor_difference_CD,neighbor_difference_DC,C_to_C,C_to_D,D_to_C,D_to_D,trade_CC,trade_CD,trade_DC,trade_DD,clustering_C,clustering_D,steady_SL_cooperation,steady_SL_defection,steady_IL_cooperation,steady_IL_defection,steady_outcome]=program(b_all,mu_all,probIL_all,c,N,w_all,degree,scale,time_all,number_structural);






%output the total runtime of the program (in seconds)
elapsed_time = toc;
disp(['The total runtime of the program: ', num2str(elapsed_time), ' seconds']);