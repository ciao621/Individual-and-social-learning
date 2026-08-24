using Random
using Statistics


@inline function count_C(i, strategy, neighbors_list)
    count = 0
    for nb in neighbors_list[i]
        count += strategy[nb]
    end
    return count
end


@inline function node_payoff(i, strategy, neighbors_list, degree_vector, b, c)
    return b * count_C(i, strategy, neighbors_list) - c * degree_vector[i] * strategy[i]
end


function cluster_probabilities(strategy, neighbors_list)
    q_C_given_C_sum = 0.0
    q_C_given_D_sum = 0.0
    n_C_nodes = 0
    n_D_nodes = 0

    for i in eachindex(strategy)
        degree = length(neighbors_list[i])
        if degree == 0
            continue
        end

        frac_C_neighbors = count_C(i, strategy, neighbors_list) / degree
        if strategy[i] == 1
            q_C_given_C_sum += frac_C_neighbors
            n_C_nodes += 1
        else
            q_C_given_D_sum += frac_C_neighbors
            n_D_nodes += 1
        end
    end

    q_C_given_C = n_C_nodes > 0 ? q_C_given_C_sum / n_C_nodes : NaN
    q_C_given_D = n_D_nodes > 0 ? q_C_given_D_sum / n_D_nodes : NaN
    return q_C_given_C, q_C_given_D
end


function validate_lambda_weights(lambda_weights, mu)
    if mu <= 0
        error("mu must be positive.")
    end

    weights = Float64.(lambda_weights)
    if length(weights) != mu
        error("lambda_weights must have length mu.")
    end
    if any(weight -> weight < 0.0, weights)
        error("lambda_weights must be non-negative.")
    end
    if !isapprox(sum(weights), 1.0; atol=1e-10, rtol=1e-10)
        error("lambda_weights must sum to one.")
    end

    return weights
end


function run_simulation(;
    neighbors_list,
    degree_vector,
    b::Float64,
    c::Float64,
    mu::Int,
    lambda_weights,
    p_trial::Float64,
    w::Float64,
    time_all::Int,
    n_runs::Int,
)
    # Strategy encoding: 0 = defection, 1 = cooperation.
    # p_trial is the probability of individual learning, P^IL, in the paper.
    lambda_weights = validate_lambda_weights(lambda_weights, mu)
    degree_vector = Int.(vec(degree_vector))
    N = length(neighbors_list)
    max_C = maximum(degree_vector)

    # transition_counts[update_type, pre_state + 1, pre_C + 1, post_state + 1]
    # update_type: 1 = reflection, 2 = individual-learning trial, 3 = social learning.
    transition_counts = zeros(Int, 3, 2, max_C + 1, 2)

    # trial_results[old_strategy + 1, C_start + 1, C_end + 1, reverted + 1]
    # old_strategy: 0 = D tried C, 1 = C tried D.
    # reverted: 1 = exploratory action kept, 2 = exploratory action reverted.
    trial_results = zeros(Int, 2, max_C + 1, max_C + 1, 2)

    steady_outcome_sum = 0.0
    q_C_given_C_sum = 0.0
    q_C_given_D_sum = 0.0
    q_C_given_C_count = 0
    q_C_given_D_count = 0

    trial_active = zeros(Bool, N)
    trial_start_C = zeros(Int, N)
    trial_old_strategy = zeros(Int, N)
    trial_original_payoff = zeros(Float64, N)
    trial_cognition = zeros(Float64, N)

    for _ in 1:n_runs
        strategy_vector = ones(Int, N)
        strategy_vector[randperm(N)[1:div(N, 2)]] .= 0

        study_vector = zeros(Int, N)
        fill!(trial_active, false)
        fill!(trial_cognition, 0.0)

        # Transient stage: evolve the system without recording statistics.
        for _ in 1:time_all
            strategy_snapshot = copy(strategy_vector)

            # Select the focal individual first, so reselection interrupts an
            # active trial before that trial collects the current payoff.
            dm = rand(1:N)
            if trial_active[dm]
                trial_active[dm] = false
                study_vector[dm] = 0
                trial_cognition[dm] = 0.0
            end

            # Remaining active individual learners collect f(t+1), ..., f(t+mu).
            for i in findall(trial_active)
                experience_index = mu - study_vector[i] + 1
                payoff_now = node_payoff(i, strategy_snapshot, neighbors_list, degree_vector, b, c)
                trial_cognition[i] += lambda_weights[experience_index] * payoff_now
            end

            # Reflection after mu collected payoff experiences.
            study_over = findall(==(1), study_vector)
            for i in study_over
                if trial_cognition[i] <= trial_original_payoff[i]
                    strategy_vector[i] = trial_old_strategy[i]
                end
                trial_active[i] = false
                study_vector[i] = 0
                trial_cognition[i] = 0.0
            end

            study_node = study_vector .!= 0
            study_vector[study_node] .-= 1

            pre_state = strategy_snapshot[dm]

            if rand() < p_trial
                strategy_vector[dm] = 1 - pre_state
                study_vector[dm] = mu
                trial_active[dm] = true
                trial_start_C[dm] = count_C(dm, strategy_snapshot, neighbors_list)
                trial_old_strategy[dm] = pre_state
                trial_original_payoff[dm] = node_payoff(dm, strategy_snapshot, neighbors_list, degree_vector, b, c)
                trial_cognition[dm] = 0.0
            else
                fit_C = 0.0
                fit_D = 0.0
                for nb in neighbors_list[dm]
                    payoff_nb = node_payoff(nb, strategy_snapshot, neighbors_list, degree_vector, b, c)
                    fitness = 1.0 - w + w * payoff_nb
                    if strategy_snapshot[nb] == 1
                        fit_C += fitness
                    else
                        fit_D += fitness
                    end
                end
                prob_C = fit_C / (fit_C + fit_D)
                strategy_vector[dm] = rand() < prob_C ? 1 : 0
            end
        end

        # Recording stage: use the same update rule while collecting statistics.
        for _ in 1:time_all
            strategy_snapshot = copy(strategy_vector)

            # Select the focal individual first, so reselection interrupts an
            # active trial before that trial collects the current payoff.
            dm = rand(1:N)
            if trial_active[dm]
                trial_active[dm] = false
                study_vector[dm] = 0
                trial_cognition[dm] = 0.0
            end

            # Remaining active individual learners collect payoff experiences.
            for i in findall(trial_active)
                experience_index = mu - study_vector[i] + 1
                payoff_now = node_payoff(i, strategy_snapshot, neighbors_list, degree_vector, b, c)
                trial_cognition[i] += lambda_weights[experience_index] * payoff_now
            end

            # Reflection compares experience-based cognition with the old payoff.
            study_over = findall(==(1), study_vector)
            for i in study_over
                pre_state = strategy_snapshot[i]
                pre_C = count_C(i, strategy_snapshot, neighbors_list)
                reverted = trial_cognition[i] <= trial_original_payoff[i]

                if reverted
                    strategy_vector[i] = trial_old_strategy[i]
                end

                transition_counts[1, pre_state + 1, pre_C + 1, strategy_vector[i] + 1] += 1

                if trial_active[i]
                    cs = min(trial_start_C[i], max_C) + 1
                    ce = min(pre_C, max_C) + 1
                    trial_results[trial_old_strategy[i] + 1, cs, ce, Int(reverted) + 1] += 1
                end

                trial_active[i] = false
                study_vector[i] = 0
                trial_cognition[i] = 0.0
            end

            study_node = study_vector .!= 0
            study_vector[study_node] .-= 1

            pre_state = strategy_snapshot[dm]
            pre_C = count_C(dm, strategy_snapshot, neighbors_list)

            if rand() < p_trial
                strategy_vector[dm] = 1 - pre_state
                transition_counts[2, pre_state + 1, pre_C + 1, strategy_vector[dm] + 1] += 1

                study_vector[dm] = mu
                trial_active[dm] = true
                trial_start_C[dm] = pre_C
                trial_old_strategy[dm] = pre_state
                trial_original_payoff[dm] = node_payoff(dm, strategy_snapshot, neighbors_list, degree_vector, b, c)
                trial_cognition[dm] = 0.0
            else
                fit_C = 0.0
                fit_D = 0.0
                for nb in neighbors_list[dm]
                    payoff_nb = node_payoff(nb, strategy_snapshot, neighbors_list, degree_vector, b, c)
                    fitness = 1.0 - w + w * payoff_nb
                    if strategy_snapshot[nb] == 1
                        fit_C += fitness
                    else
                        fit_D += fitness
                    end
                end
                prob_C = fit_C / (fit_C + fit_D)
                strategy_vector[dm] = rand() < prob_C ? 1 : 0
                transition_counts[3, pre_state + 1, pre_C + 1, strategy_vector[dm] + 1] += 1
            end

            steady_outcome_sum += mean(strategy_vector)
            q_C_given_C_t, q_C_given_D_t = cluster_probabilities(strategy_vector, neighbors_list)
            if !isnan(q_C_given_C_t)
                q_C_given_C_sum += q_C_given_C_t
                q_C_given_C_count += 1
            end
            if !isnan(q_C_given_D_t)
                q_C_given_D_sum += q_C_given_D_t
                q_C_given_D_count += 1
            end
        end
    end

    steady_outcome = steady_outcome_sum / (time_all * n_runs)
    q_C_given_C = q_C_given_C_count > 0 ? q_C_given_C_sum / q_C_given_C_count : NaN
    q_C_given_D = q_C_given_D_count > 0 ? q_C_given_D_sum / q_C_given_D_count : NaN

    return (;
        steady_outcome,
        transition_counts,
        trial_results,
        q_C_given_C,
        q_C_given_D,
        max_C,
    )
end
