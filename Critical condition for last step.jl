using LinearAlgebra

"""
    compute_PIL_finalstep(k, Pt_C, Pomega_C, DPomega_C, Pt_D, Pomega_D, DPomega_D, mu, b, c)

Compute `P_{C→D}^{IL}` and `P_{D→C}^{IL}` (and their derivatives) for the case `λ = (0,...,0,1)`

# Inputs
- `k`          : number of neighbors
- `Pt_C`       : initial distribution of #cooperative neighbors for a focal cooperator
- `Pomega_C`   : neighbor-count transition matrix for a focal cooperator 
- `DPomega_C`  : derivative of `Pomega_C`
- `Pt_D`       : initial distribution for a focal defector (length k+1)
- `Pomega_D`   : neighbor-count transition matrix for a focal defector
- `DPomega_D`  : derivative of `Pomega_D`
- `mu`         : number of trial-and-error steps
- `b, c`       : benefit and cost parameters

# Outputs
- `P_CD_IL`    : probability that a cooperator becomes a defector after trial-and-error
- `P_DC_IL`    : probability that a defector becomes a cooperator after trial-and-error
- `DP_CD_IL`   : derivative of `P_CD_IL`
- `DP_DC_IL`   : derivative of `P_DC_IL`
"""
function compute_PIL_finalstep(
    k::Int,
    Pt_C::AbstractVector,
    Pomega_C::AbstractMatrix,
    DPomega_C::AbstractMatrix,
    Pt_D::AbstractVector,
    Pomega_D::AbstractMatrix,
    DPomega_D::AbstractMatrix,
    mu::Int,
    b::Real,
    c::Real
)
    n = k + 1
    @assert length(Pt_C) == n && length(Pt_D) == n
    @assert size(Pomega_C) == (n, n) && size(Pomega_D) == (n, n)
    @assert size(DPomega_C) == (n, n) && size(DPomega_D) == (n, n)
    @assert mu ≥ 0

    # 1) Distribution after μ trial-and-error steps
    Pomega_C_mu = Pomega_C^mu
    Pomega_D_mu = Pomega_D^mu

    # 2) Build Heaviside matrices G (size (k+1)×(k+1))
    G_C = zeros(Float64, n, n)
    G_D = zeros(Float64, n, n)

    for i in 1:n
        for j in 1:n
            delta_C = b * ((j - 1) - (i - 1)) + (2 * 1 - 1) * c * k  
            delta_D = b * ((j - 1) - (i - 1)) + (2 * 0 - 1) * c * k  
            G_C[i, j] = (delta_C > 0) ? 1.0 : 0.0
            G_D[i, j] = (delta_D > 0) ? 1.0 : 0.0
        end
    end

    # 3) Derivative of matrix power: D(P^μ) = Σ_{i=0}^{μ-1} P^i * DP * P^{μ-1-i}
    DPomega_C_mu = zeros(Float64, n, n)
    DPomega_D_mu = zeros(Float64, n, n)
    if mu > 0
        for i in 0:(mu - 1)
            DPomega_C_mu .+= (Pomega_C^i) * DPomega_C * (Pomega_C^(mu - 1 - i))
            DPomega_D_mu .+= (Pomega_D^i) * DPomega_D * (Pomega_D^(mu - 1 - i))
        end
    end

    # 4) Compute probabilities
    onesvec = ones(Float64, n)
    vC   = (G_C .* Pomega_C_mu) * onesvec
    vD   = (G_D .* Pomega_D_mu) * onesvec
    vC_d = (G_C .* DPomega_C_mu) * onesvec
    vD_d = (G_D .* DPomega_D_mu) * onesvec

    P_CD_IL  = dot(Pt_C, vC)
    P_DC_IL  = dot(Pt_D, vD)
    DP_CD_IL = dot(Pt_C, vC_d)
    DP_DC_IL = dot(Pt_D, vD_d)

    return P_CD_IL, P_DC_IL, DP_CD_IL, DP_DC_IL
end


"""
    compute_PIL_structured(k, mu, b, c, PIL, x_zero, q1, N)

Compute `P_{C→D}^{IL}` and `P_{D→C}^{IL}` (and derivatives) for `λ=(0,...,0,1)`

# Inputs
- `k`      : number of neighbors
- `mu`     : number of trial-and-error steps
- `b, c`   : benefit and cost parameters
- `PIL`    : individual-learning probability
- `q1`     : q^(1) (used to build pairwise conditional probabilities)
- `N`      : population size

# Outputs
- `P_CD_IL`, `P_DC_IL`, `DP_CD_IL`, `DP_DC_IL`
"""
function compute_PIL_structured(
    k::Int,
    mu::Int,
    b::Real,
    c::Real,
    PIL::Real,
    x_zero::Real,
    q1::Real,
    N::Int
)
    n = k + 1
    @assert mu ≥ 0
    @assert N ≥ 2

    # --- Step 1: initial distribution P(t) ---
    q_DD = q1
    q_CD = 1 - q_DD
    q_CC = q1
    q_DC = 1 - q_CC

    Pt_C = zeros(Float64, n)
    Pt_D = zeros(Float64, n)
    for kc in 0:k
        Pt_C[kc + 1] = binomial(k, kc) * (q_CC^kc) * (q_DC^(k - kc))
        Pt_D[kc + 1] = binomial(k, kc) * (q_CD^kc) * (q_DD^(k - kc))
    end

    # --- Step 2: social-learning matrix entries P_SL ---
    P_CC_SL_C = (k - 1) / k * q_CC
    P_CD_SL_C = (k - 1) / k * q_DC + 1 / k
    P_DC_SL_C = (k - 1) / k * q_CD
    P_DD_SL_C = (k - 1) / k * q_DD + 1 / k
    P_CC_SL_D = (k - 1) / k * q_CC + 1 / k
    P_CD_SL_D = (k - 1) / k * q_DC
    P_DC_SL_D = (k - 1) / k * q_CD + 1 / k
    P_DD_SL_D = (k - 1) / k * q_DD

    # --- Step 3: build transition matrices P_omega ---
    Pomega_C = zeros(Float64, n, n)
    Pomega_D = zeros(Float64, n, n)

    PRE = (PIL / N) * (1 - 1 / N)^mu

    for kc in 0:k
        kd = k - kc

        PCC_SL = 1 - (1 - (1 - PIL) * (kd + 1) * kc / k / N)^mu
        PDD_SL = (1 - (1 - PIL) * (kc + 1) * kd / k / N)^mu

        if kc < k
            Pomega_C[kc + 1, kc + 2] =
                kd * ((1 - PIL) / N * P_DC_SL_C + PIL / N) + (kd * PCC_SL - kc * PDD_SL) * PRE
            Pomega_D[kc + 1, kc + 2] =
                kd * ((1 - PIL) / N * P_DC_SL_D + PIL / N) - (kd * PCC_SL - kc * PDD_SL) * PRE
        end
        if kc > 0
            Pomega_C[kc + 1, kc] =
                kc * ((1 - PIL) / N * P_CD_SL_C + PIL / N) - (kd * PCC_SL - kc * PDD_SL) * PRE
            Pomega_D[kc + 1, kc] =
                kc * ((1 - PIL) / N * P_CD_SL_D + PIL / N) + (kd * PCC_SL - kc * PDD_SL) * PRE
        end

        Pomega_C[kc + 1, kc + 1] =
            kc * (1 - PIL) / N * P_CC_SL_C + kd * (1 - PIL) / N * P_DD_SL_C + 1 - k / N
        Pomega_D[kc + 1, kc + 1] =
            kc * (1 - PIL) / N * P_CC_SL_D + kd * (1 - PIL) / N * P_DD_SL_D + 1 - k / N
    end

    # --- Step 4: derivatives of social-learning entries DP_SL ---
    v1 = b * (k - 1) * q_CC - c * k
    v2 = b * (k - 1) * q_CD
    v3 = b * ((k - 1) * q_CC + 1) - c * k
    v4 = b * ((k - 1) * q_CD + 1)

    # Trial from C to D
    DP_CC_SL_C =  (v3 - v4) / k^2 * ((k - 1) * (k - 2) * q_CC * q_DC + (k - 1) * q_CC)
    DP_CD_SL_C = -(v3 - v4) / k^2 * ((k - 1) * (k - 2) * q_CC * q_DC + (k - 1) * q_CC)
    DP_DC_SL_C =  (v1 - v2) / k^2 * ((k - 1) * (k - 2) * q_CD * q_DD + (k - 1) * q_CD)
    DP_DD_SL_C = -(v1 - v2) / k^2 * ((k - 1) * (k - 2) * q_CD * q_DD + (k - 1) * q_CD)

    # Trial from D to C
    DP_CC_SL_D =  (v3 - v4) / k^2 * ((k - 1) * (k - 2) * q_CC * q_DC + (k - 1) * q_DC)
    DP_CD_SL_D = -(v3 - v4) / k^2 * ((k - 1) * (k - 2) * q_CC * q_DC + (k - 1) * q_DC)
    DP_DC_SL_D =  (v1 - v2) / k^2 * ((k - 1) * (k - 2) * q_CD * q_DD + (k - 1) * q_DD)
    DP_DD_SL_D = -(v1 - v2) / k^2 * ((k - 1) * (k - 2) * q_CD * q_DD + (k - 1) * q_DD)

    # --- Step 5: build derivative transition matrices DP_omega ---
    DPomega_C = zeros(Float64, n, n)
    DPomega_D = zeros(Float64, n, n)

    for kc in 0:k
        kd = k - kc
        if kc < k
            DPomega_C[kc + 1, kc + 2] = kd / N * ((1 - PIL) * DP_DC_SL_C)
            DPomega_D[kc + 1, kc + 2] = kd / N * ((1 - PIL) * DP_DC_SL_D)
        end
        if kc > 0
            DPomega_C[kc + 1, kc] = kc / N * ((1 - PIL) * DP_CD_SL_C)
            DPomega_D[kc + 1, kc] = kc / N * ((1 - PIL) * DP_CD_SL_D)
        end
        DPomega_C[kc + 1, kc + 1] = kc / N * ((1 - PIL) * DP_CC_SL_C) + kd / N * ((1 - PIL) * DP_DD_SL_C)
        DPomega_D[kc + 1, kc + 1] = kc / N * ((1 - PIL) * DP_CC_SL_D) + kd / N * ((1 - PIL) * DP_DD_SL_D)
    end

    # --- Step 6: call the final-step routine ---
    return compute_PIL_finalstep(k, Pt_C, Pomega_C, DPomega_C, Pt_D, Pomega_D, DPomega_D, mu, b, c)
end


# =============================================================================
# Figure 2A
# =============================================================================

w     = 0.01        # selection intensity (unused below, kept for completeness)
N     = 500         # population size
k     = 4           # number of neighbors
c     = 1.0         # cost
omega = 0.01        # selection intensity
x_zero = 0.5

PIL   = 0.05        # individual-learning probability
mu_all = 100:20:600
q1_all = [
    0.630050112, 0.631057271, 0.631949063, 0.632737045, 0.633432151, 0.63404462,
    0.634583948, 0.63505885,  0.635477248, 0.635846322, 0.636172508, 0.636461553,
    0.63671856,  0.636948025, 0.637153901, 0.637339648, 0.637508265, 0.637662367,
    0.637804203, 0.637935694, 0.638058509, 0.638174046, 0.638283498, 0.638387875,
    0.638488019, 0.638584632
]

L = length(mu_all)
Der_pc  = zeros(Float64, L)
Der_IL  = zeros(Float64, L)
Der_DIL = zeros(Float64, L)
Der_SL  = zeros(Float64, L)

linjie = zeros(Float64, L)

for (idx, mu) in enumerate(mu_all)
    PRE  = PIL * (1 - 1 / N)^mu            # reflection probability (as in MATLAB)
    q1 = q1_all[idx]

    b = 3.5                                # initial benefit
    while Der_pc[idx] ≤ 0
        b += 0.1

        P_CD_IL, P_DC_IL, DP_CD_IL, DP_DC_IL = compute_PIL_structured(k, mu, b, c, PIL, x_zero, q1, N)


        Der_IL[idx]  = PRE / omega * (P_DC_IL * (1 - x_zero) - P_CD_IL * x_zero)
        Der_DIL[idx] = PRE * (DP_DC_IL * (1 - x_zero) - DP_CD_IL * x_zero)
        Der_SL[idx]  = (1 - PIL) * (k - 1) / k * x_zero * (2q1 - 2q1^2) * (b * (k - 1) * (2q1 - 1) - c * k)

        Der_pc[idx] = Der_SL[idx] + Der_IL[idx] + Der_DIL[idx]
    end

    linjie[idx] = b
end

println("PIL=0.05, ",linjie)



# =============================================================================
# Figure 2E
# =============================================================================

PIL   = 0.10        # individual-learning probability
mu_all = [400, 500, 600]
q1_all = [0.614321194, 0.615299132, 0.616000753]

L = length(mu_all)
Der_pc  = zeros(Float64, L)
Der_IL  = zeros(Float64, L)
Der_DIL = zeros(Float64, L)
Der_SL  = zeros(Float64, L)

linjie = zeros(Float64, L)

for (idx, mu) in enumerate(mu_all)
    PRE  = PIL * (1 - 1 / N)^mu            # reflection probability (as in MATLAB)
    q1 = q1_all[idx]

    b = 3.5                                # initial benefit
    while Der_pc[idx] ≤ 0
        b += 0.1

        P_CD_IL, P_DC_IL, DP_CD_IL, DP_DC_IL = compute_PIL_structured(k, mu, b, c, PIL, x_zero, q1, N)


        Der_IL[idx]  = PRE / omega * (P_DC_IL * (1 - x_zero) - P_CD_IL * x_zero)
        Der_DIL[idx] = PRE * (DP_DC_IL * (1 - x_zero) - DP_CD_IL * x_zero)
        Der_SL[idx]  = (1 - PIL) * (k - 1) / k * x_zero * (2q1 - 2q1^2) * (b * (k - 1) * (2q1 - 1) - c * k)

        Der_pc[idx] = Der_SL[idx] + Der_IL[idx] + Der_DIL[idx]
    end

    linjie[idx] = b
end

println("PIL=0.10, ",linjie)


PIL   = 0.15        # individual-learning probability
mu_all = [400, 500, 600]
q1_all = [0.597611045, 00.598755507, 0.599565577]

L = length(mu_all)
Der_pc  = zeros(Float64, L)
Der_IL  = zeros(Float64, L)
Der_DIL = zeros(Float64, L)
Der_SL  = zeros(Float64, L)
linjie = zeros(Float64, L)

for (idx, mu) in enumerate(mu_all)
    PRE  = PIL * (1 - 1 / N)^mu            # reflection probability (as in MATLAB)
    q1 = q1_all[idx]

    b = 3.5                                # initial benefit
    while Der_pc[idx] ≤ 0
        b += 0.1

        P_CD_IL, P_DC_IL, DP_CD_IL, DP_DC_IL = compute_PIL_structured(k, mu, b, c, PIL, x_zero, q1, N)


        Der_IL[idx]  = PRE / omega * (P_DC_IL * (1 - x_zero) - P_CD_IL * x_zero)
        Der_DIL[idx] = PRE * (DP_DC_IL * (1 - x_zero) - DP_CD_IL * x_zero)
        Der_SL[idx]  = (1 - PIL) * (k - 1) / k * x_zero * (2q1 - 2q1^2) * (b * (k - 1) * (2q1 - 1) - c * k)

        Der_pc[idx] = Der_SL[idx] + Der_IL[idx] + Der_DIL[idx]
    end

    linjie[idx] = b
end

println("PIL=0.15, ",linjie)
