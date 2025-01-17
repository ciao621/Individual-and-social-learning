
tic;

%population size
N=500;
%intensity of selection, w<<1/N
w=0.01;
%duration of transient and steady states.
time_all=50*N;

%degree
degree=4;
%payoff parameter
b=5;
c=1;

%parameter mu can be any distribution. In practice, we only need to generate mu_all according to a specific distribution.
mu_ave=20;
mu_all=randi([mu_ave-min(1000-mu_ave,mu_ave-1), mu_ave+min(1000-mu_ave,mu_ave-1)], 1, N);

%initial ratio of cooperators
initial=0.5;

%time interval for recording the evolving ratio of cooperators
scale=N/2;

%probability of individual learning
probIL=tanh(w*[c*degree]);

%number of network structures
number_structural=5;



[steady_evolution,strategy_change_before,strategy_change_after,payoff_own_C,payoff_other_C,payoff_own_D,payoff_other_D]=program_b(b,mu_all,probIL,c,N,w,degree,initial,scale,time_all,number_structural);




%Fig.4e
%{
number_all=1;
condition_all=zeros(2,4);
%The first column is the number of times β^own>β^other>0, and the second column is the number of times the opposite result occurs.
outcome=zeros(1,2);

%number of network structures
number_structural=1;


while number_all <= 20

    [steady_evolution,strategy_change_before,strategy_change_after,payoff_own_C,payoff_other_C,payoff_own_D,payoff_other_D]=program_b(b,mu_all,probIL,c,N,w,degree,initial,scale,time_all,number_structural);




    %own流
    data_C_own=payoff_own_C;
    data_D_own=payoff_own_D;

    data_C_other=payoff_other_C;
    data_D_other=payoff_other_D;

    data_befor1=strategy_change_before;
    data_after1=strategy_change_after;

    delt_own1=data_C_own;
    delt_own1(:,:)=99;

    delt_other1=data_C_other;
    delt_other1(:,:)=99;

    c_before1=data_befor1;
    c_before1(:,:)=99;
    c_after1=data_after1;
    c_after1(:,:)=99;

    for i=1:size(data_C_own,1)
        idx_C_own=find(data_C_own(i,:)~=99);
        idx_D_own=find(data_D_own(i,:)~=99);
        idx_own=intersect(idx_C_own,idx_D_own);

        idx_C_other=find(data_C_other(i,:)~=99);
        idx_D_other=find(data_D_other(i,:)~=99);
        idx_other=intersect(idx_C_other,idx_D_other);

        idx=intersect(idx_own,idx_other);

        data_colu_own=data_C_own(i,idx)-data_D_own(i,idx);
        delt_own1(i,1:length(data_colu_own))=data_colu_own;

        data_colu_other=data_C_other(i,idx)-data_D_other(i,idx);
        delt_other1(i,1:length(data_colu_other))=data_colu_other;

        data_colu_before=data_befor1(i,idx);
        c_before1(i,1:length(data_colu_before))=data_colu_before;

        data_colu_after=data_after1(i,idx);
        c_after1(i,1:length(data_colu_after))=data_colu_after;

    end

    idx=zeros(1,size(data_C_own,1));
    for i=1:size(data_C_own,1)

        idx_own=find(delt_own1(i,:)==99,1)-1;
        idx_other=find(delt_other1(i,:)==99,1)-1;
        idx_i=min(idx_other,idx_own);
        idx(i)=idx_i;
    end
    idx=min(idx);


    delt_own1=delt_own1(:,1:idx);
    delt_other1=delt_other1(:,1:idx);
    c_before1=c_before1(:,1:idx);
    c_after1=c_after1(:,1:idx);

    delt_own1 = reshape(delt_own1, [], 1);
    delt_other1 = reshape(delt_other1, [], 1);
    c_before1 = reshape(c_before1, [], 1);
    c_after1 = reshape(c_after1, [], 1);






    % Create a table (similar to DataFrame in pandas)
    T = table(c_after1, c_before1, delt_own1, delt_other1);
    % Fit a Generalized Linear Model (GLM) with a binomial distribution
    % Assuming the response variable is binary, fit a logistic regression model
    mdl = fitglm(T, 'c_after1 ~ c_before1 + delt_own1 + delt_other1', 'Distribution', 'binomial');

    condition1_1=mdl.Coefficients{3,1}>0;
    condition1_2=mdl.Coefficients{4,1}>0;
    condition1_3=mdl.Coefficients{3,1}>mdl.Coefficients{4,1};
    condition1=condition1_1*condition1_2*condition1_3;

    condition2_1=mdl.Coefficients{3,4}<0.05;
    condition2_2=mdl.Coefficients{4,4}<0.05;
    condition2=condition2_1*condition2_2;
    condition=condition1*condition2;

    condition_all=zeros(2,4);
    condition2_3=condition2*condition1_1*condition1_2*(1-condition1_3);
    condition_all(1,1)=condition_all(1,1)+1-condition1_1;
    condition_all(1,2)=condition_all(1,2)+1-condition1_2;
    condition_all(1,3)=condition_all(1,3)+1-condition1_3;
    condition_all(1,4)=condition_all(1,4)+1;
    condition_all(2,1)=condition_all(2,1)+1-condition2_1;
    condition_all(2,2)=condition_all(2,2)+1-condition2_2;
    condition_all(2,3)=condition_all(2,3)+condition2_3;


    while condition ~= 1 %|| condition_value ~= 1


        [steady_evolution,strategy_change_before,strategy_change_after,payoff_own_C,payoff_other_C,payoff_own_D,payoff_other_D]=program_b(b,mu_all,probIL,c,N,w,degree,initial,scale,time_all,number_structural);




        %own流
        data_C_own=payoff_own_C;
        data_D_own=payoff_own_D;

        data_C_other=payoff_other_C;
        data_D_other=payoff_other_D;

        data_befor1=strategy_change_before;
        data_after1=strategy_change_after;

        delt_own1=data_C_own;
        delt_own1(:,:)=99;

        delt_other1=data_C_other;
        delt_other1(:,:)=99;

        c_before1=data_befor1;
        c_before1(:,:)=99;
        c_after1=data_after1;
        c_after1(:,:)=99;

        for i=1:size(data_C_own,1)
            idx_C_own=find(data_C_own(i,:)~=99);
            idx_D_own=find(data_D_own(i,:)~=99);
            idx_own=intersect(idx_C_own,idx_D_own);

            idx_C_other=find(data_C_other(i,:)~=99);
            idx_D_other=find(data_D_other(i,:)~=99);
            idx_other=intersect(idx_C_other,idx_D_other);

            idx=intersect(idx_own,idx_other);

            data_colu_own=data_C_own(i,idx)-data_D_own(i,idx);
            delt_own1(i,1:length(data_colu_own))=data_colu_own;

            data_colu_other=data_C_other(i,idx)-data_D_other(i,idx);
            delt_other1(i,1:length(data_colu_other))=data_colu_other;

            data_colu_before=data_befor1(i,idx);
            c_before1(i,1:length(data_colu_before))=data_colu_before;

            data_colu_after=data_after1(i,idx);
            c_after1(i,1:length(data_colu_after))=data_colu_after;

        end

        idx=zeros(1,size(data_C_own,1));
        for i=1:size(data_C_own,1)

            idx_own=find(delt_own1(i,:)==99,1)-1;
            idx_other=find(delt_other1(i,:)==99,1)-1;
            idx_i=min(idx_other,idx_own);
            idx(i)=idx_i;
        end
        idx=min(idx);

        delt_own1=delt_own1(:,1:idx);
        delt_other1=delt_other1(:,1:idx);
        c_before1=c_before1(:,1:idx);
        c_after1=c_after1(:,1:idx);

        delt_own1 = reshape(delt_own1, [], 1);
        delt_other1 = reshape(delt_other1, [], 1);
        c_before1 = reshape(c_before1, [], 1);
        c_after1 = reshape(c_after1, [], 1);






        % Create a table (similar to DataFrame in pandas)
        T = table(c_after1, c_before1, delt_own1, delt_other1);
        % Fit a Generalized Linear Model (GLM) with a binomial distribution
        % Assuming the response variable is binary, fit a logistic regression model
        mdl = fitglm(T, 'c_after1 ~ c_before1 + delt_own1 + delt_other1', 'Distribution', 'binomial');

        condition1_1=mdl.Coefficients{3,1}>0;
        condition1_2=mdl.Coefficients{4,1}>0;
        condition1_3=mdl.Coefficients{3,1}>mdl.Coefficients{4,1};
        condition1=condition1_1*condition1_2*condition1_3;

        condition2_1=mdl.Coefficients{3,4}<0.05;
        condition2_2=mdl.Coefficients{4,4}<0.05;
        condition2=condition2_1*condition2_2;
        condition=condition1*condition2;


        condition2_3=condition2*condition1_1*condition1_2*(1-condition1_3);
        condition_all(1,1)=condition_all(1,1)+1-condition1_1;
        condition_all(1,2)=condition_all(1,2)+1-condition1_2;
        condition_all(1,3)=condition_all(1,3)+1-condition1_3;
        condition_all(1,4)=condition_all(1,4)+1;
        condition_all(2,1)=condition_all(2,1)+1-condition2_1;
        condition_all(2,2)=condition_all(2,2)+1-condition2_2;
        condition_all(2,3)=condition_all(2,3)+condition2_3;

    end

    if condition_all(2,3) == 0
        outcome(1)=outcome(1)+1;
    else
        outcome(2)=outcome(2)+1;
    end
    number_all=number_all+1;
end

%}




%}
%}

% 在完成程序后调用 toc 函数
elapsed_time = toc;

% 显示花费的时间（以秒为单位）
disp(['程序运行时间: ', num2str(elapsed_time), ' 秒']);