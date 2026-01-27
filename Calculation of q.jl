using SparseArrays
using LinearAlgebra
using IterativeSolvers  # idrs
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





"""
    f_value(PIL, N, k, q) -> Float64

Compute 2*(1 - ((1-PIL)/N)*q)^(2k) - 2*(1 - ((1-PIL)/N)*q)^k + 1
in a numerically stable way.
"""
function f_value(PIL::Real, N::Integer, k::Integer)
    @assert N ≥ 1
    @assert k ≥ 0

    # q 
    deg = 4
    q0 = 1 / deg
    # y = ((1-PIL)/N)*q
    y = ((1 - float(PIL)) * float(deg * q0 + q0 - deg * q0^2)) / float(N)

    # need 1-y > 0 for log1p(-y) (i.e., y < 1)
    @assert y < 1 "Require ((1-PIL)/N)*q < 1 so that 1-y is positive."

    loga = log1p(-y)          # log(1 - y)  (stable when y is small)
    powk = exp(k * loga)      # (1 - y)^k
    pow2k = powk * powk       # (1 - y)^(2k)

    return muladd(2.0, pow2k, muladd(-2.0, powk, 1.0))  # 2*pow2k - 2*powk + 1
end


"""
    weights_uniform(mu) -> Vector{Float64}

λ_k = 1/mu
"""
weights_uniform(mu::Integer) = fill(1.0 / mu, mu)

"""
    weights_geometric(mu, r) -> Vector{Float64}

λ_{k+1}/λ_k = r, normalized so sum(λ)=1.
For r=1 it reduces to uniform.
"""
function weights_geometric(mu::Integer, r::Real)
    @assert mu ≥ 1
    r = float(r)
    if isapprox(r, 1.0; atol=0, rtol=1e-12)
        return weights_uniform(mu)
    end
    # λ_k ∝ r^(k-1)
    # normalize using closed-form sum to avoid overflow
    # sum_{k=1..mu} r^(k-1) = (1 - r^mu)/(1 - r)
    denom = (1 - r^mu) / (1 - r)
    return [r^(k - 1) / denom for k in 1:mu]
end

"""
    weighted_sum_f(PTE, N, mu, weights) -> Float64

Compute Σ_{k=1..mu} λ_k f(k).
"""
function weighted_sum_f(PTE::Real, N::Integer, mu::Integer, λ::AbstractVector{<:Real})
    @assert length(λ) == mu
    # optional sanity check:
    # @assert abs(sum(λ) - 1) ≤ 1e-10
    s = 0.0
    for k in 1:mu
        s += float(λ[k]) * f_value(PTE, N, k)
    end
    return s
end






"""
Last step
"""
function tau(adjancey_matrix::AbstractMatrix, PTE::Float64, mu::Int; row_stochastic::Bool=true)

    n1, n2 = size(adjancey_matrix)
    @assert n1 == n2 "Adjacency must be square"
    N = n1

    # work in Float64 sparse
    M = sparse(Float64.(adjancey_matrix))  # if A already sparse, this keeps it sparse
    deg = sum(M, dims=1)
    P = M ./ deg
    Pcol = sparse(transpose(P))  # column i of Pcol == row i of P

    # coefficients in Eq.(10)
    PRE = (mu > 0) ? PTE * (1 - 1 / N)^mu : 0.0
    # g_xing = 2 * (1 - (1 - PTE) / N * (deg[1] / (2 * (deg[1] - 1))))^(2 * mu) - 2 * (1 - (1 - PTE) / N * (deg[1] / (2 * (deg[1] - 1))))^(mu) + 1
    q0 = 1 / deg[1]
    g_xing = 2 * (1 - (1 - PTE) / N * (deg[1] * q0 + q0 - deg[1] * q0^2))^(2 * mu) - 2 * (1 - (1 - PTE) / N * (deg[1] * q0 + q0 - deg[1] * q0^2))^(mu) + 1
    PREg = PRE * g_xing
    #PREg = PRE * (2^(-1/N))^mu
    denom = 1 + PTE + 2 * PREg

    c = (PTE + PREg) / denom
    d = (1 - PTE) / (2 * denom)

    # index in vec(τ) (column-major)
    idx(i, j) = i + (j - 1) * N

    # diagonal mask
    isdiag = falses(N^2)
    @inbounds for i in 1:N
        isdiag[idx(i, i)] = true
    end

    off_idx = findall(!, isdiag)
    n_off = length(off_idx)

    # map full index -> position in off vector
    map = zeros(Int, N^2)
    @inbounds for (pos, fid) in enumerate(off_idx)
        map[fid] = pos
    end

    # RHS
    b = fill(c, n_off)

    # sparse triplets for A_off
    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, n_off * (2d + 1))
    sizehint!(J, n_off * (2d + 1))
    sizehint!(V, n_off * (2d + 1))

    @inline function full_to_ij(fid::Int, N::Int)
        j = (fid - 1) ÷ N + 1
        i = fid - (j - 1) * N
        return i, j
    end

    @inbounds for r in 1:n_off
        fid = off_idx[r]
        i, j = full_to_ij(fid, N)  # equation for τ_ij (i≠j)

        # + τ_ij
        push!(I, r)
        push!(J, r)
        push!(V, 1.0)

        # term: -d * Σ_k p_{ik} τ_{jk}
        ks, vals = findnz(Pcol[:, i])  # row i of P
        for t in eachindex(ks)
            k = ks[t]
            pik = vals[t]
            fid_var = idx(j, k)        # τ_{jk}
            if isdiag[fid_var]
                # τ_{jj}=1 when k=j
                b[r] += d * pik
            else
                push!(I, r)
                push!(J, map[fid_var])
                push!(V, -d * pik)
            end
        end

        # term: -d * Σ_k p_{jk} τ_{ik}
        ks2, vals2 = findnz(Pcol[:, j]) # row j of P
        for t in eachindex(ks2)
            k = ks2[t]
            pjk = vals2[t]
            fid_var = idx(i, k)         # τ_{ik}
            if isdiag[fid_var]
                # τ_{ii}=1 when k=i
                b[r] += d * pjk
            else
                push!(I, r)
                push!(J, map[fid_var])
                push!(V, -d * pjk)
            end
        end
    end

    Aoff = sparse(I, J, V, n_off, n_off)

    # solve and reconstruct τ
    xoff = idrs(Aoff, b)
    tau_vec = ones(Float64, N^2)     # τ_ii = 1
    tau_vec[off_idx] = xoff
    τ_matrix = reshape(tau_vec, N, N)

    # ---- tau_out = (1/N) * sum_{i,j} p_ij * tau_ij  (sparse-safe) ----
    s = 0.0
    @inbounds for j in 1:N
        for ptr in P.colptr[j]:(P.colptr[j+1]-1)
            i = P.rowval[ptr]
            s += P.nzval[ptr] * τ_matrix[i, j]
        end
    end
    τ_out = s / N



    return τ_matrix, τ_out
end



"""
λ_{k+1}/λ_k = 1
""" 
function tau1(adjancey_matrix::AbstractMatrix, PTE::Float64, mu::Int; row_stochastic::Bool=true)

    n1, n2 = size(adjancey_matrix)
    @assert n1 == n2 "Adjacency must be square"
    N = n1

    # work in Float64 sparse
    M = sparse(Float64.(adjancey_matrix))  # if A already sparse, this keeps it sparse
    deg = sum(M, dims=1)
    P = M ./ deg
    Pcol = sparse(transpose(P))  # column i of Pcol == row i of P

    # coefficients in Eq.(10)
    PRE = (mu > 0) ? PTE * (1 - 1 / N)^mu : 0.0
    λmu = weights_uniform(mu)
    Smu = weighted_sum_f(PTE, N, mu, λmu)
    PREg = PRE * Smu
    #PREg = PRE * (2^(-1/N))^mu
    denom = 1 + PTE + 2 * PREg

    c = (PTE + PREg) / denom
    d = (1 - PTE) / (2 * denom)

    # index in vec(τ) (column-major)
    idx(i, j) = i + (j - 1) * N

    # diagonal mask
    isdiag = falses(N^2)
    @inbounds for i in 1:N
        isdiag[idx(i, i)] = true
    end

    off_idx = findall(!, isdiag)
    n_off = length(off_idx)

    # map full index -> position in off vector
    map = zeros(Int, N^2)
    @inbounds for (pos, fid) in enumerate(off_idx)
        map[fid] = pos
    end

    # RHS
    b = fill(c, n_off)

    # sparse triplets for A_off
    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, n_off * (2d + 1))
    sizehint!(J, n_off * (2d + 1))
    sizehint!(V, n_off * (2d + 1))

    @inline function full_to_ij(fid::Int, N::Int)
        j = (fid - 1) ÷ N + 1
        i = fid - (j - 1) * N
        return i, j
    end

    @inbounds for r in 1:n_off
        fid = off_idx[r]
        i, j = full_to_ij(fid, N)  # equation for τ_ij (i≠j)

        # + τ_ij
        push!(I, r)
        push!(J, r)
        push!(V, 1.0)

        # term: -d * Σ_k p_{ik} τ_{jk}
        ks, vals = findnz(Pcol[:, i])  # row i of P
        for t in eachindex(ks)
            k = ks[t]
            pik = vals[t]
            fid_var = idx(j, k)        # τ_{jk}
            if isdiag[fid_var]
                # τ_{jj}=1 when k=j
                b[r] += d * pik
            else
                push!(I, r)
                push!(J, map[fid_var])
                push!(V, -d * pik)
            end
        end

        # term: -d * Σ_k p_{jk} τ_{ik}
        ks2, vals2 = findnz(Pcol[:, j]) # row j of P
        for t in eachindex(ks2)
            k = ks2[t]
            pjk = vals2[t]
            fid_var = idx(i, k)         # τ_{ik}
            if isdiag[fid_var]
                # τ_{ii}=1 when k=i
                b[r] += d * pjk
            else
                push!(I, r)
                push!(J, map[fid_var])
                push!(V, -d * pjk)
            end
        end
    end

    Aoff = sparse(I, J, V, n_off, n_off)

    # solve and reconstruct τ
    xoff = idrs(Aoff, b)
    tau_vec = ones(Float64, N^2)     # τ_ii = 1
    tau_vec[off_idx] = xoff
    τ_matrix = reshape(tau_vec, N, N)

    # ---- tau_out = (1/N) * sum_{i,j} p_ij * tau_ij  (sparse-safe) ----
    s = 0.0
    @inbounds for j in 1:N
        for ptr in P.colptr[j]:(P.colptr[j+1]-1)
            i = P.rowval[ptr]
            s += P.nzval[ptr] * τ_matrix[i, j]
        end
    end
    τ_out = s / N



    return τ_matrix, τ_out
end



"""
λ_{k+1}/λ_k = 2
""" 
function tau2(adjancey_matrix::AbstractMatrix, PTE::Float64, mu::Int; row_stochastic::Bool=true)

    n1, n2 = size(adjancey_matrix)
    @assert n1 == n2 "Adjacency must be square"
    N = n1

    # work in Float64 sparse
    M = sparse(Float64.(adjancey_matrix))  # if A already sparse, this keeps it sparse
    deg = sum(M, dims=1)
    P = M ./ deg
    Pcol = sparse(transpose(P))  # column i of Pcol == row i of P

    # coefficients in Eq.(10)
    PRE = (mu > 0) ? PTE * (1 - 1 / N)^mu : 0.0
    λmu = weights_geometric(mu, 2.0)
    Smu = weighted_sum_f(PTE, N, mu, λmu)
    PREg = PRE * Smu
    #PREg = PRE * (2^(-1/N))^mu
    denom = 1 + PTE + 2 * PREg

    c = (PTE + PREg) / denom
    d = (1 - PTE) / (2 * denom)

    # index in vec(τ) (column-major)
    idx(i, j) = i + (j - 1) * N

    # diagonal mask
    isdiag = falses(N^2)
    @inbounds for i in 1:N
        isdiag[idx(i, i)] = true
    end

    off_idx = findall(!, isdiag)
    n_off = length(off_idx)

    # map full index -> position in off vector
    map = zeros(Int, N^2)
    @inbounds for (pos, fid) in enumerate(off_idx)
        map[fid] = pos
    end

    # RHS
    b = fill(c, n_off)

    # sparse triplets for A_off
    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, n_off * (2d + 1))
    sizehint!(J, n_off * (2d + 1))
    sizehint!(V, n_off * (2d + 1))

    @inline function full_to_ij(fid::Int, N::Int)
        j = (fid - 1) ÷ N + 1
        i = fid - (j - 1) * N
        return i, j
    end

    @inbounds for r in 1:n_off
        fid = off_idx[r]
        i, j = full_to_ij(fid, N)  # equation for τ_ij (i≠j)

        # + τ_ij
        push!(I, r)
        push!(J, r)
        push!(V, 1.0)

        # term: -d * Σ_k p_{ik} τ_{jk}
        ks, vals = findnz(Pcol[:, i])  # row i of P
        for t in eachindex(ks)
            k = ks[t]
            pik = vals[t]
            fid_var = idx(j, k)        # τ_{jk}
            if isdiag[fid_var]
                # τ_{jj}=1 when k=j
                b[r] += d * pik
            else
                push!(I, r)
                push!(J, map[fid_var])
                push!(V, -d * pik)
            end
        end

        # term: -d * Σ_k p_{jk} τ_{ik}
        ks2, vals2 = findnz(Pcol[:, j]) # row j of P
        for t in eachindex(ks2)
            k = ks2[t]
            pjk = vals2[t]
            fid_var = idx(i, k)         # τ_{ik}
            if isdiag[fid_var]
                # τ_{ii}=1 when k=i
                b[r] += d * pjk
            else
                push!(I, r)
                push!(J, map[fid_var])
                push!(V, -d * pjk)
            end
        end
    end

    Aoff = sparse(I, J, V, n_off, n_off)

    # solve and reconstruct τ
    xoff = idrs(Aoff, b)
    tau_vec = ones(Float64, N^2)     # τ_ii = 1
    tau_vec[off_idx] = xoff
    τ_matrix = reshape(tau_vec, N, N)

    # ---- tau_out = (1/N) * sum_{i,j} p_ij * tau_ij  (sparse-safe) ----
    s = 0.0
    @inbounds for j in 1:N
        for ptr in P.colptr[j]:(P.colptr[j+1]-1)
            i = P.rowval[ptr]
            s += P.nzval[ptr] * τ_matrix[i, j]
        end
    end
    τ_out = s / N



    return τ_matrix, τ_out
end



"""
λ_{k+1}/λ_k = 0.5
""" 
function tau05(adjancey_matrix::AbstractMatrix, PTE::Float64, mu::Int; row_stochastic::Bool=true)

    n1, n2 = size(adjancey_matrix)
    @assert n1 == n2 "Adjacency must be square"
    N = n1

    # work in Float64 sparse
    M = sparse(Float64.(adjancey_matrix))  # if A already sparse, this keeps it sparse
    deg = sum(M, dims=1)
    P = M ./ deg
    Pcol = sparse(transpose(P))  # column i of Pcol == row i of P

    # coefficients in Eq.(10)
    PRE = (mu > 0) ? PTE * (1 - 1 / N)^mu : 0.0
    λmu = weights_geometric(mu, 0.5)
    Smu = weighted_sum_f(PTE, N, mu, λmu)
    PREg = PRE * Smu
    #PREg = PRE * (2^(-1/N))^mu
    denom = 1 + PTE + 2 * PREg

    c = (PTE + PREg) / denom
    d = (1 - PTE) / (2 * denom)

    # index in vec(τ) (column-major)
    idx(i, j) = i + (j - 1) * N

    # diagonal mask
    isdiag = falses(N^2)
    @inbounds for i in 1:N
        isdiag[idx(i, i)] = true
    end

    off_idx = findall(!, isdiag)
    n_off = length(off_idx)

    # map full index -> position in off vector
    map = zeros(Int, N^2)
    @inbounds for (pos, fid) in enumerate(off_idx)
        map[fid] = pos
    end

    # RHS
    b = fill(c, n_off)

    # sparse triplets for A_off
    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, n_off * (2d + 1))
    sizehint!(J, n_off * (2d + 1))
    sizehint!(V, n_off * (2d + 1))

    @inline function full_to_ij(fid::Int, N::Int)
        j = (fid - 1) ÷ N + 1
        i = fid - (j - 1) * N
        return i, j
    end

    @inbounds for r in 1:n_off
        fid = off_idx[r]
        i, j = full_to_ij(fid, N)  # equation for τ_ij (i≠j)

        # + τ_ij
        push!(I, r)
        push!(J, r)
        push!(V, 1.0)

        # term: -d * Σ_k p_{ik} τ_{jk}
        ks, vals = findnz(Pcol[:, i])  # row i of P
        for t in eachindex(ks)
            k = ks[t]
            pik = vals[t]
            fid_var = idx(j, k)        # τ_{jk}
            if isdiag[fid_var]
                # τ_{jj}=1 when k=j
                b[r] += d * pik
            else
                push!(I, r)
                push!(J, map[fid_var])
                push!(V, -d * pik)
            end
        end

        # term: -d * Σ_k p_{jk} τ_{ik}
        ks2, vals2 = findnz(Pcol[:, j]) # row j of P
        for t in eachindex(ks2)
            k = ks2[t]
            pjk = vals2[t]
            fid_var = idx(i, k)         # τ_{ik}
            if isdiag[fid_var]
                # τ_{ii}=1 when k=i
                b[r] += d * pjk
            else
                push!(I, r)
                push!(J, map[fid_var])
                push!(V, -d * pjk)
            end
        end
    end

    Aoff = sparse(I, J, V, n_off, n_off)

    # solve and reconstruct τ
    xoff = idrs(Aoff, b)
    tau_vec = ones(Float64, N^2)     # τ_ii = 1
    tau_vec[off_idx] = xoff
    τ_matrix = reshape(tau_vec, N, N)

    # ---- tau_out = (1/N) * sum_{i,j} p_ij * tau_ij  (sparse-safe) ----
    s = 0.0
    @inbounds for j in 1:N
        for ptr in P.colptr[j]:(P.colptr[j+1]-1)
            i = P.rowval[ptr]
            s += P.nzval[ptr] * τ_matrix[i, j]
        end
    end
    τ_out = s / N



    return τ_matrix, τ_out
end



N = 500
d = 4
PTE = 0.05
adjancey_matrix = create_regular_graph(N, d)


mu_all = 100:20:600
q_vector = zeros(length(mu_all))
for (i, mu) in enumerate(mu_all)
    τ_matrix, τ_out = tau(adjancey_matrix, PTE, mu)
    q_vector[i] = τ_out
end
println("Stationary average IIS probability (last step): ", q_vector)
