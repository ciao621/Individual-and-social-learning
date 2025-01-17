# Interdependence-between-individual-and-social-learning-shape-cooperation-with-asymmetric-exploration
This study analyzes a theoretical model using both mathematical and numerical methods. Here, we provide a computer simulation program.


%main_b can calculate the evolution data of the cooperator ratio and the results of Fig.4e, main_a can calculate all other results.

%The estimated running time of main_a is more than 2*10^5 seconds, which is much longer than the calculation of fixed probability. 
%This is because the evolutionary process we are focusing on needs to calculate the average frequency of cooperators in the steady state. Each run needs to simulate 2000*N steps as transient state and 2000*N steps as steady state. 
%The program can be run more than ten times faster by converting program_c to a mex file.

%Our program focuses on the random regular network
%For any heterogeneous network, just delete the part of the program that generates the random regular network and load the specific adjacency matrix adjacency_matrix

%All image parameters are either fixed values ​​as described in Supplementary Table 1 or follow the ranges described below.
%Fig.4c and Supplementary Fig.S3b, b_all=[0.5,1,1.5,2,2.5,3 3.5 4 4.5 5 5.5 6 6.5 7 7.5]
b_all=[3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13];
probIL_all=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95];
w_all=0.01:0.01:0.20;
%Note that this program only supports a single range parameter, i.e. if b_all takes the above range, probIL_all and w_all are fixed values.

%Total number of updates
number_all=sum(steady_SL_cooperation)+sum(steady_SL_defection)+sum(steady_IL_cooperation)+sum(steady_IL_defection);

%Exploration direction
number_trade_C=sum(steady_IL_cooperation(3:4));
number_trade_D=sum(steady_IL_defection(3:4));
number_trial_C=steady_IL_cooperation(2);
number_trial_D=steady_IL_defect(2);

radio_trade_C=number_trade_C/number_all;
radio_trade_D=number_trade_D/number_all;
radio_trial_C=number_trial_C/number_all;
radio_trial_D=number_trial_D/number_all;

%Assortment
assortment=clustering_C-clustering_D;

%mutation frequency,f_{C\toD}=mutation_C(5),f_{D\toC}=mutation_D(1)
mutation_C=C_to_D./(C_to_C+C_to_D);
mutation_D=D_to_C./(D_to_C+D_to_D);

%probability of playing C; the i-th column indicates that there are i-1 cooperators among the neighbors
probability_C=C_to_C./(C_to_C+C_to_D);
probability_D=D_to_C./(D_to_C+D_to_D);

%probability of playing C in the trade-off phase of individual learning; the i-th column indicates that there are i-1 cooperators among the neighbors
probability_IL_C=trade_CC./(trade_CC+trade_CD);
probability_IL_D=trade_DC./(trade_DC+trade_DD);

%coopeartion frequency in social learning
coopeartion_SL=sum(steady_SL_cooperation)./(sum(steady_SL_cooperation)+sum(steady_SL_defection));

%coopeartion frequency in individual learning, i.e.,the colors of phase diagram of individual learning
coopeartion_IL=sum(steady_IL_cooperation)./(sum(steady_IL_cooperation)+sum(steady_IL_defection));

%proportion of defectors in trial-and-error
propdefector_IL=steady_IL_cooperation(2)./(steady_IL_cooperation(2)+steady_IL_defection(2));

%trade-off frequency,f_{C\toD}=fretrade_CD(5),f_{D\toC}=fretrade_DC(1)
fretrade_CD=trade_CD./(C_to_C+C_to_D);
fretrade_DC=trade_DC./(D_to_C+D_to_D);


%the colors of phase diagram of trade-off choice
coopeartion_trade=sum(steady_IL_cooperation(3:4))./(sum(steady_IL_cooperation(3:4))+sum(steady_IL_defection(3:4)));
