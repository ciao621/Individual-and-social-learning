function [neighbor_difference_CD,neighbor_difference_DC,C_to_C,C_to_D,D_to_C,D_to_D,trade_CC,trade_CD,trade_DC,trade_DD,clustering_C,clustering_D,steady_SL_cooperation,steady_SL_defection,steady_IL_cooperation,steady_IL_defection,steady_outcome]=program_a(b_all,mu_all,probIL_all,c,N,w_all,degree,scale,time_all,number_structural)

%Data collection
parameter_length=max([length(b_all),length(probIL_all),length(w_all)]);
if length(b_all) == 1
    b_all = repmat(b_all, 1, parameter_length);
end
if length(probIL_all) == 1
    probIL_all = repmat(probIL_all, 1, parameter_length);
end
if length(w_all) == 1
    w_all = repmat(w_all, 1, parameter_length);
end
parameter_matrix = [b_all; probIL_all; w_all];


%Output
%Strategy update through social learning, odd columns are the number of times the strategy remains unchanged, even columns are the number of times the strategy changes.
steady_SL_cooperation=zeros(parameter_length,2);
steady_SL_defection=zeros(parameter_length,2);

%Strategy updates through individual learning, the first two columns are trial-and-error updates, the last two columns are trade-off updates, the odd columns are the number of times the strategy remains unchanged, even columns are the number of times the strategy changes.
steady_IL_cooperation=zeros(parameter_length,4);
steady_IL_defection=zeros(parameter_length,4);

steady_outcome=zeros(parameter_length,1);
clustering_C=zeros(parameter_length,1);
clustering_D=zeros(parameter_length,1);
neighbor_difference_CD=zeros(parameter_length,1);
neighbor_difference_DC=zeros(parameter_length,1);

%The number of strategy switches; the i-th row indicates that there are i-1 cooperators among the neighbors.
C_to_D=zeros(parameter_length,5);
C_to_C=zeros(parameter_length,5);
D_to_D=zeros(parameter_length,5);
D_to_C=zeros(parameter_length,5);

%The number of strategy switches in the trade-off phase of individual learning; the i-th row indicates that there are i-1 cooperators among the neighbors.
trade_CC=zeros(parameter_length,5);
trade_CD=zeros(parameter_length,5);
trade_DC=zeros(parameter_length,5);
trade_DD=zeros(parameter_length,5);



% Initialization
% Strategy initial value, half cooperation and half defection,1 represents a cooperator, 0 represents a defector
strategy_initial=ones(1,N);
indices=randperm(N,N/2);
strategy_initial(indices)=0;









for number_parameter=1:parameter_length


    b=parameter_matrix(1,number_parameter);
    probIL=parameter_matrix(2,number_parameter);
    w=parameter_matrix(3,number_parameter);
    

    %Data collection
    number_IL_CD=0;
    number_IL_DC=0;





    %Change the network structure
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
        neighbor_cooperation_old=zeros(1,N);

        %Compute the initial payoff vector
        for i=1:N
            neighbors_focal=adjacency_matrix(:,i) == 1;
            payoff_vector(i)=sum(strategy_vector(neighbors_focal))*b-strategy_vector(i)*c*degree;
        end

        time_ing=1;

        %Transient
        while time_ing <= time_all

            % The trade-off phase of individual learning
            study_node=study_vector~=0;
            payoff_study((study_node))=(payoff_study((study_node))+payoff_vector(study_node))/2;
            study_over=find(study_vector == 1);
            for i=1:length(study_over)
                if payoff_study(study_over(i)) <= payoff_study_old(study_over(i))
                    strategy_vector(study_over(i))=1-strategy_vector(study_over(i));
                    neighbors_focal=find(adjacency_matrix(:,study_over(i)) == 1);
                    %payoff update
                    payoff_vector(study_over(i))=payoff_vector(study_over(i))+(1-2*strategy_vector(study_over(i)))*c*degree;
                    payoff_vector(neighbors_focal)=payoff_vector(neighbors_focal)+(2*strategy_vector(study_over(i))-1)*b;
                end
            end
            payoff_study(study_over)=0;
            study_vector(study_node)=study_vector(study_node)-1;


            %decision making
            decision_game_maker=randi(N);

            if study_vector(decision_game_maker) ~= 0

                %If the selected node is in the process of individual learning, there will be no strategy update this turn.

                time_ing=time_ing+1;

                continue;

            else


                game_old=strategy_vector(decision_game_maker);
                neighbors_focal=find(adjacency_matrix(:,decision_game_maker) == 1);
                neighbors_focal_cooperation=neighbors_focal(strategy_vector(neighbors_focal)==1);


                if rand < probIL

                    %The trial-and-error phase of individual learning

                    strategy_vector(decision_game_maker)=1-strategy_vector(decision_game_maker);
                    payoff_study(decision_game_maker)=payoff_vector(decision_game_maker)+(1-2*strategy_vector(decision_game_maker))*c*degree;
                    study_vector(decision_game_maker)=mu_all(decision_game_maker);
                    payoff_study_old(decision_game_maker)=payoff_vector(decision_game_maker);

                else

                    %social learning
                    payoff_neighbors_cooperation=payoff_vector(neighbors_focal_cooperation);
                    fitness_focal_cooperation=1-w+w*payoff_neighbors_cooperation;
                    neighbors_focal_defect=neighbors_focal(strategy_vector(neighbors_focal)==0);
                    payoff_neighbors_defect=payoff_vector(neighbors_focal_defect);
                    fitness_focal_defect=1-w+w*payoff_neighbors_defect;
                    probability_cooperation=sum(fitness_focal_cooperation)/(sum(fitness_focal_cooperation)+sum(fitness_focal_defect));
                    probability_cooperation=max(probability_cooperation,0);

                    if rand < probability_cooperation

                        strategy_vector(decision_game_maker)=1;

                    else

                        strategy_vector(decision_game_maker)=0;

                    end

                end

                %payoff update
                if game_old ~= strategy_vector(decision_game_maker)

                    payoff_vector(decision_game_maker)=payoff_vector(decision_game_maker)+(1-2*strategy_vector(decision_game_maker))*c*degree;

                    payoff_vector(neighbors_focal)=payoff_vector(neighbors_focal)+(2*strategy_vector(decision_game_maker)-1)*b;

                end

                time_ing=time_ing+1;
            end

        end


        %Steady-state
        time_ing=1;
        while time_ing <= time_all

            % The trade-off phase of individual learning
            study_node=study_vector~=0;
            payoff_study((study_node))=(payoff_study((study_node))+payoff_vector(study_node))/2;
            study_over=find(study_vector == 1);


            for i=1:length(study_over)

                neighbors_focal=find(adjacency_matrix(:,study_over(i)) == 1);
                neighbors_focal_cooperation=neighbors_focal(strategy_vector(neighbors_focal)==1);

                %Data collection
                if strategy_vector(study_over(i))==0
                    neighbor_difference_CD(number_parameter)=neighbor_difference_CD(number_parameter)+length(neighbors_focal_cooperation)-neighbor_cooperation_old(study_over(i));
                elseif strategy_vector(study_over(i))==1
                    neighbor_difference_DC(number_parameter)=neighbor_difference_DC(number_parameter)+length(neighbors_focal_cooperation)-neighbor_cooperation_old(study_over(i));

                end

                if payoff_study(study_over(i)) < payoff_study_old(study_over(i))
                    strategy_vector(study_over(i))=1-strategy_vector(study_over(i));
                    neighbors_focal=find(adjacency_matrix(:,study_over(i)) == 1);
                    %payoff update
                    payoff_vector(study_over(i))=payoff_vector(study_over(i))+(1-2*strategy_vector(study_over(i)))*c*degree;
                    payoff_vector(neighbors_focal)=payoff_vector(neighbors_focal)+(2*strategy_vector(study_over(i))-1)*b;

                    %Data collection
                    steady_IL_cooperation(number_parameter,4)=steady_IL_cooperation(number_parameter,4)+strategy_vector(study_over(i));
                    steady_IL_defection(number_parameter,4)=steady_IL_defection(number_parameter,4)+1-strategy_vector(study_over(i));
                    C_to_D(number_parameter,length(neighbors_focal_cooperation)+1)=C_to_D(number_parameter,length(neighbors_focal_cooperation)+1)+1-strategy_vector(study_over(i));
                    D_to_C(number_parameter,length(neighbors_focal_cooperation)+1)=D_to_C(number_parameter,length(neighbors_focal_cooperation)+1)+strategy_vector(study_over(i));
                    trade_CD(number_parameter,length(neighbors_focal_cooperation)+1)=trade_CD(number_parameter,length(neighbors_focal_cooperation)+1)+1-strategy_vector(study_over(i));
                    trade_DC(number_parameter,length(neighbors_focal_cooperation)+1)=trade_DC(number_parameter,length(neighbors_focal_cooperation)+1)+strategy_vector(study_over(i));
                else
                    %Data collection
                    steady_IL_cooperation(number_parameter,3)=steady_IL_cooperation(number_parameter,3)+strategy_vector(study_over(i));
                    steady_IL_defection(number_parameter,3)=steady_IL_defection(number_parameter,3)+1-strategy_vector(study_over(i));
                    C_to_C(number_parameter,length(neighbors_focal_cooperation)+1)=C_to_C(number_parameter,length(neighbors_focal_cooperation)+1)+strategy_vector(study_over(i));
                    D_to_D(number_parameter,length(neighbors_focal_cooperation)+1)=D_to_D(number_parameter,length(neighbors_focal_cooperation)+1)+1-strategy_vector(study_over(i));
                    trade_CC(number_parameter,length(neighbors_focal_cooperation)+1)=trade_CC(number_parameter,length(neighbors_focal_cooperation)+1)+strategy_vector(study_over(i));
                    trade_DD(number_parameter,length(neighbors_focal_cooperation)+1)=trade_DD(number_parameter,length(neighbors_focal_cooperation)+1)+1-strategy_vector(study_over(i));
                end
            end
            payoff_study(study_over)=0;
            study_vector(study_node)=study_vector(study_node)-1;




            %decision making
            decision_game_maker=randi(N);

            if study_vector(decision_game_maker) ~= 0

                %If the selected node is in the process of individual learning, there will be no strategy update this turn.
                steady_outcome(number_parameter)=steady_outcome(number_parameter)+mean(strategy_vector);

                %Data collection
                if mod(time_ing,scale) == 0
                    clustering_C_ing=0;
                    cooperator_node=find(strategy_vector==1);

                    for i=1:length(cooperator_node)
                        neighbors_focal=find(adjacency_matrix(:,cooperator_node(i)) == 1);
                        neighbors_focal_cooperation=neighbors_focal(strategy_vector(neighbors_focal)==1);
                        clustering_C_ing=clustering_C_ing+length(neighbors_focal_cooperation)/degree;
                    end
                    clustering_C(number_parameter)=clustering_C(number_parameter)+clustering_C_ing/max(length(cooperator_node),1);
                    clustering_D_ing=0;
                    defector_node=find(strategy_vector==0);

                    for i=1:length(defector_node)
                        neighbors_focal=find(adjacency_matrix(:,defector_node(i)) == 1);
                        neighbors_focal_cooperation=neighbors_focal(strategy_vector(neighbors_focal)==1);
                        clustering_D_ing=clustering_D_ing+length(neighbors_focal_cooperation)/degree;
                    end

                    clustering_D(number_parameter)=clustering_D(number_parameter)+clustering_D_ing/max(length(defector_node),1);
                end

                time_ing=time_ing+1;

                continue;

            else


                game_old=strategy_vector(decision_game_maker);
                neighbors_focal=find(adjacency_matrix(:,decision_game_maker) == 1);
                neighbors_focal_cooperation=neighbors_focal(strategy_vector(neighbors_focal)==1);

                if rand < probIL

                    %The trial-and-error phase of individual learning
                    strategy_vector(decision_game_maker)=1-strategy_vector(decision_game_maker);
                    payoff_study(decision_game_maker)=payoff_vector(decision_game_maker)+(1-2*strategy_vector(decision_game_maker))*c*degree;
                    neighbor_cooperation_old(decision_game_maker)=length(neighbors_focal_cooperation);

                    %Data collection
                    steady_IL_cooperation(number_parameter,2)=steady_IL_cooperation(number_parameter,2)+strategy_vector(decision_game_maker);
                    steady_IL_defection(number_parameter,2)=steady_IL_defection(number_parameter,2)+1-strategy_vector(decision_game_maker);
                    if game_old==1
                        number_IL_CD=number_IL_CD+1;
                    elseif game_old==0
                        number_IL_DC=number_IL_DC+1;
                    end

                    study_vector(decision_game_maker)=mu_all(decision_game_maker);
                    payoff_study_old(decision_game_maker)=payoff_vector(decision_game_maker);

                else
                    %social learning
                    payoff_neighbors_cooperation=payoff_vector(neighbors_focal_cooperation);
                    fitness_focal_cooperation=1-w+w*payoff_neighbors_cooperation;
                    neighbors_focal_defect=neighbors_focal(strategy_vector(neighbors_focal)==0);
                    payoff_neighbors_defect=payoff_vector(neighbors_focal_defect);
                    fitness_focal_defect=1-w+w*payoff_neighbors_defect;
                    probability_cooperation=sum(fitness_focal_cooperation)/(sum(fitness_focal_cooperation)+sum(fitness_focal_defect));
                    probability_cooperation=max(probability_cooperation,0);

                    if rand < probability_cooperation

                        strategy_vector(decision_game_maker)=1;

                    else

                        strategy_vector(decision_game_maker)=0;

                    end

                    %Data collection
                    if game_old ~= strategy_vector(decision_game_maker)

                        steady_SL_cooperation(number_parameter,2)=steady_SL_cooperation(number_parameter,2)+strategy_vector(decision_game_maker);
                        steady_SL_defection(number_parameter,2)=steady_SL_defection(number_parameter,2)+1-strategy_vector(decision_game_maker);

                    else

                        steady_SL_cooperation(number_parameter,1)=steady_SL_cooperation(number_parameter,1)+strategy_vector(decision_game_maker);
                        steady_SL_defection(number_parameter,1)=steady_SL_defection(number_parameter,1)+1-strategy_vector(decision_game_maker);

                    end
                end

                %payoff update
                if game_old ~= strategy_vector(decision_game_maker)

                    payoff_vector(decision_game_maker)=payoff_vector(decision_game_maker)+(1-2*strategy_vector(decision_game_maker))*c*degree;
                    payoff_vector(neighbors_focal)=payoff_vector(neighbors_focal)+(2*strategy_vector(decision_game_maker)-1)*b;

                end

                %Data collection
                if game_old ~= strategy_vector(decision_game_maker)

                    C_to_D(number_parameter,length(neighbors_focal_cooperation)+1)=C_to_D(number_parameter,length(neighbors_focal_cooperation)+1)+1-strategy_vector(decision_game_maker);
                    D_to_C(number_parameter,length(neighbors_focal_cooperation)+1)=D_to_C(number_parameter,length(neighbors_focal_cooperation)+1)+strategy_vector(decision_game_maker);

                else

                    C_to_C(number_parameter,length(neighbors_focal_cooperation)+1)=C_to_C(number_parameter,length(neighbors_focal_cooperation)+1)+strategy_vector(decision_game_maker);
                    D_to_D(number_parameter,length(neighbors_focal_cooperation)+1)=D_to_D(number_parameter,length(neighbors_focal_cooperation)+1)+1-strategy_vector(decision_game_maker);

                end


                steady_outcome(number_parameter)=steady_outcome(number_parameter)+mean(strategy_vector);

                if mod(time_ing,scale) == 0
                    clustering_C_ing=0;
                    cooperator_node=find(strategy_vector==1);
                    for i=1:length(cooperator_node)
                        neighbors_focal=find(adjacency_matrix(:,cooperator_node(i)) == 1);
                        neighbors_focal_cooperation=neighbors_focal(strategy_vector(neighbors_focal)==1);
                        clustering_C_ing=clustering_C_ing+length(neighbors_focal_cooperation)/degree;
                    end
                    clustering_C(number_parameter)=clustering_C(number_parameter)+clustering_C_ing/max(length(cooperator_node),1);
                    clustering_D_ing=0;
                    defector_node=find(strategy_vector==0);
                    for i=1:length(defector_node)
                        neighbors_focal=find(adjacency_matrix(:,defector_node(i)) == 1);
                        neighbors_focal_cooperation=neighbors_focal(strategy_vector(neighbors_focal)==1);
                        clustering_D_ing=clustering_D_ing+length(neighbors_focal_cooperation)/degree;
                    end
                    clustering_D(number_parameter)=clustering_D(number_parameter)+clustering_D_ing/max(length(defector_node),1);
                end

                time_ing=time_ing+1;

            end

        end






    end

    neighbor_difference_CD(number_parameter)=neighbor_difference_CD(number_parameter)/max(number_IL_CD,1);
    neighbor_difference_DC(number_parameter)=neighbor_difference_DC(number_parameter)/max(number_IL_DC,1);

end

%Data collection
steady_outcome(:)=steady_outcome(:)/(time_all*number_structural);
clustering_C=clustering_C/(time_all/scale*number_structural);
clustering_D=clustering_D/(time_all/scale*number_structural);

end
