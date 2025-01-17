function [steady_evolution,strategy_change_before,strategy_change_after,payoff_own_C,payoff_other_C,payoff_own_D,payoff_other_D]=program_b(b,mu_all,probIL,c,N,w,degree,initial,scale,time_all,number_structural)


%输出结果





%evolving ratio of cooperators
steady_evolution=zeros(number_structural,time_all/scale+1);


strategy_change_before=zeros(number_structural,time_all+1);
strategy_change_after=zeros(number_structural,time_all+1);
payoff_own_C=zeros(number_structural,time_all+1);
payoff_other_C=zeros(number_structural,time_all+1);
payoff_own_D=zeros(number_structural,time_all+1);
payoff_other_D=zeros(number_structural,time_all+1);





% Initialization
% Strategy initial value,1 represents a cooperator, 0 represents a defector
strategy_initial=zeros(1,N);
indices=randperm(N,floor(N*initial));
strategy_initial(indices)=1;

steady_evolution(:,1)=mean(strategy_initial);

%99 means that the data does not exist at this time step
strategy_change_before(:,:)=99;
strategy_change_after(:,:)=99;
payoff_own_C(:,:)=99;
payoff_other_C(:,:)=99;
payoff_own_D(:,:)=99;
payoff_other_D(:,:)=99;

payoff_own_C_initial=zeros(1,N);
payoff_own_C_initial(:)=99;
payoff_own_D_initial=zeros(1,N);
payoff_own_D_initial(:)=99;

payoff_other_C_initial=zeros(1,N);
payoff_other_C_initial(:)=99;
payoff_other_D_initial=zeros(1,N);
payoff_other_D_initial(:)=99;





for number1=1:number_structural

    judgment_regular=0;

    

    % createRegularGraph - creates a simple d-regular undirected graph


    % simple = without loops or double edges


    % d-reglar = each vertex is adjecent to d edges


    %


    % input arguments :


    %   vertNum - number of vertices


    %   deg - the degree of each vertex




    % output arguments :


    %   A - A sparse matrix representation of the graph


    % algorithm :


    % "The pairing model" : create n*d 'half edges'.


    % repeat as long as possible: pick a pair of half edges


    %   and if it's legal (doesn't creat a loop nor a double edge)


    %   add it to the graph



    % reference: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.67.7957&rep=rep1&type=pdf


    n = N;


    d = degree;


    matIter = 10;


    %a list of open half-edges


    U = repmat(1:n,1,d);


    %the graphs adajency matrix


    adjacency_matrix=zeros(n,n);


    edgesTested=0;


    repetition=1;


    %continue until a proper graph is formed


    while ~isempty(U) && repetition < matIter



        edgesTested = edgesTested + 1;


        %chose at random 2 half edges


        i1 = ceil(rand*length(U));


        i2 = ceil(rand*length(U));


        v1 = U(i1);


        v2 = U(i2);


        %check that there are no loops nor parallel edges


        if (v1 == v2) || (adjacency_matrix(v1,v2) == 1)



            %restart process if needed


            if (edgesTested == n*d)


                repetition=repetition+1;


                edgesTested = 0;


                U = repmat(1:n,1,d);


                adjacency_matrix = zeros(n,n);


            end


        else


            %add edge to graph


            adjacency_matrix(v1, v2)=1;


            adjacency_matrix(v2, v1)=1;



            %remove used half-edges


            v = sort([i1,i2]);


            U = [U(1:v(1)-1), U(v(1)+1:v(2)-1), U(v(2)+1:end)];


        end


    end


    degree_vector=sum(adjacency_matrix);

    for i=1:N

        if degree_vector(i)==degree 

            judgment_regular=1;

        end

    end

    while judgment_regular~=1

        % createRegularGraph - creates a simple d-regular undirected graph

        % simple = without loops or double edges

        % d-reglar = each vertex is adjecent to d edges

        %

        % input arguments :

        %   vertNum - number of vertices

        %   deg - the degree of each vertex



        % output arguments :

        %   A - A sparse matrix representation of the graph

        % algorithm :

        % "The pairing model" : create n*d 'half edges'.

        % repeat as long as possible: pick a pair of half edges

        %   and if it's legal (doesn't creat a loop nor a double edge)

        %   add it to the graph

        % reference: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.67.7957&rep=rep1&type=pdf

        n = N;

        d = degree;

        matIter = 10;

        %a list of open half-edges

        U = repmat(1:n,1,d);

        %the graphs adajency matrix

        adjacency_matrix=zeros(n,n);

        edgesTested=0;

        repetition=1;

        %continue until a proper graph is formed

        while ~isempty(U) && repetition < matIter


            edgesTested = edgesTested + 1;

            %chose at random 2 half edges

            i1 = ceil(rand*length(U));

            i2 = ceil(rand*length(U));

            v1 = U(i1);

            v2 = U(i2);

            %check that there are no loops nor parallel edges

            if (v1 == v2) || (adjacency_matrix(v1,v2) == 1)


                %restart process if needed

                if (edgesTested == n*d)

                    repetition=repetition+1;

                    edgesTested = 0;

                    U = repmat(1:n,1,d);

                    adjacency_matrix = zeros(n,n);

                end

            else

                %add edge to graph

                adjacency_matrix(v1, v2)=1;


                adjacency_matrix(v2, v1)=1;


                %remove used half-edges

                v = sort([i1,i2]);

                U = [U(1:v(1)-1), U(v(1)+1:v(2)-1), U(v(2)+1:end)];

            end

        end

        degree_vector=sum(adjacency_matrix);
        for i=1:N
            if degree_vector(i)==degree 
                judgment_regular=1;
            end
        end
    end







    % Initialization
    strategy_vector=strategy_initial;
    study_vector=zeros(1,N);
    payoff_study=zeros(1,N);
    payoff_vector=zeros(1,N);
    payoff_study_old=zeros(1,N);

    payoff_ownvector_C=payoff_own_C_initial;
    payoff_ownvector_D=payoff_own_D_initial;
    payoff_othervector_C=payoff_other_C_initial;
    payoff_othervector_D=payoff_other_D_initial;

    for i=1:N
        neighbors_focal=adjacency_matrix(:,i) == 1;
        neighbors_focal_cooperation=neighbors_focal(strategy_vector(neighbors_focal)==1);
        payoff_vector(i)=sum(strategy_vector(neighbors_focal_cooperation))*b-strategy_vector(i)*c*degree;

        payoff_neighbors_cooperation=payoff_vector(neighbors_focal_cooperation);
        neighbors_focal_defect=neighbors_focal(strategy_vector(neighbors_focal)==0);
        payoff_neighbors_defect=payoff_vector(neighbors_focal_defect);
        if isempty(neighbors_focal_cooperation) ~= 1 %
            payoff_othervector_C(i)=sum(payoff_neighbors_cooperation)/length(neighbors_focal_cooperation);
        end
        if isempty(neighbors_focal_defect) ~= 1 
            payoff_othervector_D(i)=sum(payoff_neighbors_defect)/length(neighbors_focal_defect);
        end
    end

    payoff_ownvector_C(strategy_vector==1)=payoff_vector(strategy_vector==1);
    payoff_ownvector_D(strategy_vector==0)=payoff_vector(strategy_vector==0);

    time_ing=1;


    while time_ing <= time_all

        % The trade-off phase of individual learning
        study_node=study_vector~=0;
        payoff_study((study_node))=payoff_study((study_node))+payoff_vector(study_node)/2;
        study_over=find(study_vector == 1);
        for i=1:length(study_over)
            neighbors_focal=find(adjacency_matrix(:,study_over(i)) == 1);
            neighbors_focal_cooperation=neighbors_focal(strategy_vector(neighbors_focal)==1);

            %Data collection
            payoff_neighbors_cooperation=payoff_vector(neighbors_focal_cooperation);
            neighbors_focal_defect=neighbors_focal(strategy_vector(neighbors_focal)==0);
            payoff_neighbors_defect=payoff_vector(neighbors_focal_defect);
            strategy_change_before(number1,time_ing-mu_all(study_over(i))+1)=1-strategy_vector(study_over(i));
            if strategy_vector(study_over(i)) == 1 
                payoff_own_C(number1,time_ing-mu_all(study_over(i))+1)=payoff_vector(study_over(i));
            else 
                payoff_own_D(number1,time_ing-mu_all(study_over(i))+1)=payoff_vector(study_over(i));
            end
            if isempty(neighbors_focal_cooperation) ~= 1 
                payoff_othervector_C(i)=sum(payoff_neighbors_cooperation)/length(neighbors_focal_cooperation);
            end
            if isempty(neighbors_focal_defect) ~= 1 
                payoff_othervector_D(i)=sum(payoff_neighbors_defect)/length(neighbors_focal_defect);
            end
            payoff_other_C(number1,time_ing-mu_all(study_over(i))+1)=payoff_othervector_C(study_over(i));
            payoff_other_D(number1,time_ing-mu_all(study_over(i))+1)=payoff_othervector_D(study_over(i));


            if payoff_study(study_over(i)) < payoff_study_old(study_over(i))
                strategy_vector(study_over(i))=1-strategy_vector(study_over(i));
                %payoff update
                payoff_vector(study_over(i))=payoff_vector(study_over(i))+(1-2*strategy_vector(study_over(i)))*c*degree;
                payoff_vector(neighbors_focal)=payoff_vector(neighbors_focal)+(2*strategy_vector(study_over(i))-1)*b;

            end

            %Data collection
            strategy_change_after(number1,time_ing-mu_all(study_over(i))+1)=strategy_vector(study_over(i));

        end

        payoff_study(study_over)=0;
        study_vector(study_node)=study_vector(study_node)-1;




        %decision making
        decision_game_maker=randi(N); 
        if study_vector(decision_game_maker) ~= 0

            %Data collection
            if mod(time_ing,scale)==0
                steady_evolution(number1,time_all/scale+1)=mean(strategy_vector);
            end

            time_ing=time_ing+1;
            continue;

        else

            game_old=strategy_vector(decision_game_maker);
            neighbors_focal=find(adjacency_matrix(:,decision_game_maker) == 1);
            neighbors_focal_cooperation=neighbors_focal(strategy_vector(neighbors_focal)==1);

            %Data collection
            payoff_neighbors_cooperation=payoff_vector(neighbors_focal_cooperation);
            neighbors_focal_defect=neighbors_focal(strategy_vector(neighbors_focal)==0);
            payoff_neighbors_defect=payoff_vector(neighbors_focal_defect);
            if isempty(neighbors_focal_cooperation) ~= 1 
                payoff_othervector_C(i)=sum(payoff_neighbors_cooperation)/length(neighbors_focal_cooperation);
            end
            if isempty(neighbors_focal_defect) ~= 1
                payoff_othervector_D(i)=sum(payoff_neighbors_defect)/length(neighbors_focal_defect);
            end
            payoff_other_C(number1,time_ing+1)=payoff_othervector_C(decision_game_maker);
            payoff_other_D(number1,time_ing+1)=payoff_othervector_D(decision_game_maker);

            if rand < probIL

                %The trial-and-error phase of individual learning

                %Data collection
                if strategy_vector(decision_game_maker) == 1
                    payoff_ownvector_C(decision_game_maker)=payoff_vector(decision_game_maker);
                    payoff_own_C(number1,time_ing+1)=payoff_ownvector_C(decision_game_maker);
                else
                    payoff_ownvector_D(decision_game_maker)=payoff_vector(decision_game_maker);
                    payoff_own_D(number1,time_ing+1)=payoff_ownvector_D(decision_game_maker);
                end

                strategy_vector(decision_game_maker)=1-strategy_vector(decision_game_maker);
                payoff_study(decision_game_maker)=payoff_vector(decision_game_maker)+(1-2*strategy_vector(decision_game_maker))*c*degree;
                study_vector(decision_game_maker)=mu_all(decision_game_maker);
                payoff_study_old(decision_game_maker)=payoff_vector(decision_game_maker);
               
            else  

                %social learning

                %Data collection
                payoff_own_C(number1,time_ing+1)=0;
                payoff_own_D(number1,time_ing+1)=0;
                strategy_change_before(number1,time_ing+1)=strategy_vector(decision_game_maker);

                fitness_focal_cooperation=1-w+w*payoff_neighbors_cooperation;
                fitness_focal_defect=1-w+w*payoff_neighbors_defect;
                probability_cooperation=sum(fitness_focal_cooperation)/(sum(fitness_focal_cooperation)+sum(fitness_focal_defect));
                probability_cooperation=max(probability_cooperation,0);
                if rand < probability_cooperation
                    strategy_vector(decision_game_maker)=1;
                else
                    strategy_vector(decision_game_maker)=0;
                end
                strategy_change_after(number1,time_ing+1)=strategy_vector(decision_game_maker);

            end

            %Data collection
            if mod(time_ing,scale)==0
                steady_evolution(number1,time_ing/scale+1)=mean(strategy_vector);
            end

            %payoff update
            if game_old ~= strategy_vector(decision_game_maker)
                payoff_vector(decision_game_maker)=payoff_vector(decision_game_maker)+(1-2*strategy_vector(decision_game_maker))*c*degree;
                payoff_vector(neighbors_focal)=payoff_vector(neighbors_focal)+(2*strategy_vector(decision_game_maker)-1)*b;
            end


            time_ing=time_ing+1;
        end

    end
end
end
