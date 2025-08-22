using Random
using LinearAlgebra
using Statistics
using StatsBase
using JLD2
using Dates
using Base.Threads


function create_regular_graph(n, d, max_iter=10)
    U = repeat(1:n, d)
    adjacency_matrix = zeros(Int, n, n)
    edges_tested = 0
    repetition = 1

    while !isempty(U) && repetition < max_iter
        edges_tested += 1
        i1 = rand(1:length(U))
        i2 = rand(1:length(U))
        v1 = U[i1]
        v2 = U[i2]

        if v1 != v2 && adjacency_matrix[v1, v2] == 0
            adjacency_matrix[v1, v2] = 1
            adjacency_matrix[v2, v1] = 1
            U = deleteat!(U, sort([i1, i2]))
        elseif edges_tested == n * d
            repetition += 1
            edges_tested = 0
            U = repeat(1:n, d)
            adjacency_matrix = zeros(Int, n, n)
        end
    end
    return adjacency_matrix
end


function program_c(b_all, mu, PIL, c, N, w, degree, time_all, number_structural)
  
    xiangtu_x_axis = length(b_all)


    steady_SL_cooperation = zeros(xiangtu_x_axis, 4)
    steady_SL_defect = zeros(xiangtu_x_axis, 4)


    steady_IL_cooperation = zeros(xiangtu_x_axis, 8)
    steady_IL_defect = zeros(xiangtu_x_axis, 8)

    steady_social = zeros(xiangtu_x_axis, 2)
    steady_outcome = zeros(xiangtu_x_axis)



    strategy_initial = ones(Int, N)
    indices = sample(1:N, N ÷ 2; replace=false)  
    strategy_initial[indices] .= 0  


    study_vector_initial = zeros(Int, N)

    
    payoff_vector_initial = zeros(N)

    for number_b in 1:xiangtu_x_axis
        b = b_all[number_b]

  
        for number1 in 1:number_structural
            judgment_regular = 0

            adjacency_matrix = create_regular_graph(N, degree)
            degree_vector = sum(adjacency_matrix, dims=2)

            while !all(degree_vector .== degree)
                adjacency_matrix = create_regular_graph(N, degree)
                degree_vector = sum(adjacency_matrix, dims=2)
            end




            neighbors_list = [findall(adjacency_matrix[i, :] .== 1) for i in 1:N]

            strategy_vector = copy(strategy_initial)
            study_vector = copy(study_vector_initial)
            payoff_study = copy(payoff_vector_initial)
            payoff_vector = copy(payoff_vector_initial)
            payoff_study_old = copy(payoff_vector_initial)

            adjacency_mat = adjacency_matrix .* 1 
            payoff_vector = (adjacency_mat * strategy_vector) .* b .- strategy_vector .* c * degree

   
            time_ing = 1
            while time_ing <= time_all

        
                study_node = study_vector .!= 0
                payoff_study[study_node] .= payoff_vector[study_node]
                study_over = findall(study_vector .== 1)

                for i in study_over
                    if payoff_study[i] <= payoff_study_old[i]
                        strategy_vector[i] = 1 - strategy_vector[i]
                        neighbors_focal = neighbors_list[i]
                        payoff_vector[i] += (1 - 2 * strategy_vector[i]) * c * degree
                        payoff_vector[neighbors_focal] .+= (2 * strategy_vector[i] - 1) * b

           
                        if strategy_vector[i] == 1
                            steady_IL_cooperation[number_b, 4] += 1
                        else
                            steady_IL_defect[number_b, 4] += 1
                        end
                    else
  
                        if strategy_vector[i] == 1
                            steady_IL_cooperation[number_b, 3] += 1
                        else
                            steady_IL_defect[number_b, 3] += 1
                        end
                    end
                end

                study_vector[study_node] .-= 1

 
                decision_maker = rand(1:N)
                study_vector[decision_maker] = 0
                game_old = strategy_vector[decision_maker]
                neighbors_focal = neighbors_list[decision_maker]


                if rand() < PIL
                    strategy_vector[decision_maker] = 1 - strategy_vector[decision_maker]
                    payoff_study[decision_maker] = 0


                    if strategy_vector[decision_maker] == 1
                        steady_IL_cooperation[number_b, 2] += 1
                    else
                        steady_IL_defect[number_b, 2] += 1
                    end

                    study_vector[decision_maker] = mu
                    payoff_study_old[decision_maker] = payoff_vector[decision_maker]
                else

                    Coneibor_idx = neighbors_focal[strategy_vector[neighbors_focal].==1]
                    if isempty(Coneibor_idx)
                        fitness_focal_cooperation = 0
                    else
                        payoff_neighbors_cooperation = payoff_vector[Coneibor_idx]
                        fitness_focal_cooperation = 1 .- w .+ w .* payoff_neighbors_cooperation
                    end

                    Deneibor_idx = neighbors_focal[strategy_vector[neighbors_focal].==0]
                    if isempty(Deneibor_idx)
                        fitness_focal_defect = 0
                    else
                        payoff_neighbors_defect = payoff_vector[Deneibor_idx]
                        fitness_focal_defect = 1 .- w .+ w .* payoff_neighbors_defect
                    end

                    Sumfitness_focal_cooperation = sum(fitness_focal_cooperation)
                    Sumfitness_focal_defect = sum(fitness_focal_defect)

                    probability_cooperation = Sumfitness_focal_cooperation / (Sumfitness_focal_cooperation + Sumfitness_focal_defect)
                    probability_cooperation = max(probability_cooperation, 0)

                    if rand() < probability_cooperation
                        strategy_vector[decision_maker] = 1
                    else
                        strategy_vector[decision_maker] = 0
                    end

                    if game_old != strategy_vector[decision_maker]
                        if strategy_vector[decision_maker] == 1
                            steady_SL_cooperation[number_b, 2] += 1
                        else
                            steady_SL_defect[number_b, 2] += 1
                        end
                    else
                        if strategy_vector[decision_maker] == 1
                            steady_SL_cooperation[number_b, 1] += 1
                        else
                            steady_SL_defect[number_b, 1] += 1
                        end
                    end

                    steady_social[number_b, 1] += 1
                end

                if game_old != strategy_vector[decision_maker]
                    payoff_vector[decision_maker] += (1 - 2 * strategy_vector[decision_maker]) * c * degree
                    payoff_vector[neighbors_focal] .+= (2 * strategy_vector[decision_maker] - 1) * b
                end

                time_ing += 1


            end


            time_ing = 1
            while time_ing <= time_all

                study_node = study_vector .!= 0
                payoff_study[study_node] .= payoff_vector[study_node]
                study_over = findall(study_vector .== 1)

                for i in study_over
                    if payoff_study[i] <= payoff_study_old[i]
                        strategy_vector[i] = 1 - strategy_vector[i]
                        neighbors_focal = neighbors_list[i]
                        payoff_vector[i] += (1 - 2 * strategy_vector[i]) * c * degree
                        payoff_vector[neighbors_focal] .+= (2 * strategy_vector[i] - 1) * b

  
                        if strategy_vector[i] == 1
                            steady_IL_cooperation[number_b, 8] += 1
                        else
                            steady_IL_defect[number_b, 8] += 1
                        end
                    else

                        if strategy_vector[i] == 1
                            steady_IL_cooperation[number_b, 7] += 1
                        else
                            steady_IL_defect[number_b, 7] += 1
                        end
                    end
                end

                study_vector[study_node] .-= 1

   
                decision_maker = rand(1:N)
                study_vector[decision_maker] = 0
                game_old = strategy_vector[decision_maker]
                neighbors_focal = neighbors_list[decision_maker]


                if rand() < PIL
                    strategy_vector[decision_maker] = 1 - strategy_vector[decision_maker]
                    payoff_study[decision_maker] = 0

  
                    if strategy_vector[decision_maker] == 1
                        steady_IL_cooperation[number_b, 6] += 1
                    else
                        steady_IL_defect[number_b, 6] += 1
                    end

                    study_vector[decision_maker] = mu
                    payoff_study_old[decision_maker] = payoff_vector[decision_maker]
                else

                    Coneibor_idx = neighbors_focal[strategy_vector[neighbors_focal].==1]
                    if isempty(Coneibor_idx)
                        fitness_focal_cooperation = 0
                    else
                        payoff_neighbors_cooperation = payoff_vector[Coneibor_idx]
                        fitness_focal_cooperation = 1 .- w .+ w .* payoff_neighbors_cooperation
                    end

                    Deneibor_idx = neighbors_focal[strategy_vector[neighbors_focal].==0]
                    if isempty(Deneibor_idx)
                        fitness_focal_defect = 0
                    else
                        payoff_neighbors_defect = payoff_vector[Deneibor_idx]
                        fitness_focal_defect = 1 .- w .+ w .* payoff_neighbors_defect
                    end

                    Sumfitness_focal_cooperation = sum(fitness_focal_cooperation)
                    Sumfitness_focal_defect = sum(fitness_focal_defect)



                    probability_cooperation = Sumfitness_focal_cooperation / (Sumfitness_focal_cooperation + Sumfitness_focal_defect)
                    probability_cooperation = max(probability_cooperation, 0)

                    if rand() < probability_cooperation
                        strategy_vector[decision_maker] = 1
                    else
                        strategy_vector[decision_maker] = 0
                    end

                    if game_old != strategy_vector[decision_maker]
                        if strategy_vector[decision_maker] == 1
                            steady_SL_cooperation[number_b, 4] += 1
                        else
                            steady_SL_defect[number_b, 4] += 1
                        end
                    else
                        if strategy_vector[decision_maker] == 1
                            steady_SL_cooperation[number_b, 3] += 1
                        else
                            steady_SL_defect[number_b, 3] += 1
                        end
                    end

                    steady_social[number_b, 2] += 1
                end

                if game_old != strategy_vector[decision_maker]
                    payoff_vector[decision_maker] += (1 - 2 * strategy_vector[decision_maker]) * c *degree
                    payoff_vector[neighbors_focal] .+= (2 * strategy_vector[decision_maker] - 1) * b
                end

                steady_outcome[number_b] += mean(strategy_vector)




                time_ing += 1


            end
        end
    end


    steady_outcome ./= (time_all * number_structural)
    steady_social ./= (time_all * number_structural)

    return steady_SL_cooperation, steady_SL_defect, steady_IL_cooperation, steady_IL_defect, steady_outcome, steady_social
end



function main()
    # Start timing
    start_time = now()

    # Parameter settings
    N = 500                  # Network size
    w = 0.01                 # Selection strength 
    time_all = 2000 * N      # Total simulation time
    degree = 4               # Network degree
    b_all = [3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13]  # Benefit
    c = 1                    # Cost of cooperation
    number_structural = 200    # Number of network structures to simulate

    # Trial-and-error step lengths μ
    mus = 20:20:600   
    mus = vcat(1, mus)

    # Individual learning probability
    PILs = fill(0.05, length(mus))

    # Parallel computation using Threads.@threads
    Threads.@threads for (i, (mu, PIL)) in collect(enumerate(zip(mus, PILs)))
        println("Computing parameter set $i: mu = $mu, PIL = $PIL...")

        # Call the program_c function
        steady_SL_cooperation, steady_SL_defect, steady_IL_cooperation, steady_IL_defect, steady_outcome, steady_social = program_c(b_all, mu, PIL, c, N, w, degree, time_all, number_structural)
        
        # Save results to a JLD2 file (file name includes parameter identifiers)
        filename = "results_mu$(mu)_prob$(PIL).jld2"
        @save filename steady_SL_cooperation steady_SL_defect steady_IL_cooperation steady_IL_defect steady_outcome steady_social
        
        println("Results saved to $filename")
    end

    # End timing
    elapsed_time = (now() - start_time).value / 1000  # Convert to seconds
    println("Total runtime: $elapsed_time seconds")
end

# Here, we provide a base code with parameter settings (λ_1,…,λ_(μ-1),λ_μ )=(0,…,0,1). 
# All the data presented in our study can be obtained by appropriately extending this code.

main()

