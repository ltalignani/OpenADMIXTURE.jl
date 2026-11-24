# ADMIXTURE-compatible bootstrap standard errors for Q (ADMIXTURE 1.3 manual,
# sections 2.7 and 2.9).
#
# SNPs are resampled with a moving-block bootstrap that never crosses a
# chromosome boundary (blocks keep linkage disequilibrium). Each replicate is
# re-estimated from the point estimates with the looser "minor" termination
# criterion (ADMIXTURE `-c`, default 3 iterations). The standard error of each
# q_ik is the standard deviation of its bootstrap replicates.
#
# ADMIXTURE does not document its block length, but prints it: 10 Mb converted
# to a number of SNPs with the SNP density of the data (`admixture_block`,
# checked against ADMIXTURE 1.3.0 logs). The command line uses that rule; the
# cube root of the number of SNPs (`default_block`) is the fallback when the
# positions span less than two blocks, and the default of `bootstrap_se`.

"""
    default_block(J)

Bootstrap block length (in SNPs) for `J` SNPs when positions are unknown:
`round(cbrt(J))`, the usual rate for block bootstraps.
"""
default_block(J::Integer) = max(1, round(Int, cbrt(J)))

"""
    admixture_block(chroms, positions; mb=10)

ADMIXTURE 1.3's bootstrap block length in SNPs: the number of SNPs in `mb`
megabases at the data's density, `ceil(J * mb * 10^6 / span)`, where `span`
adds up `max - min + 1` of the base-pair positions on each chromosome. When the
data span less than two blocks (ADMIXTURE then crashes, or every replicate is
the data itself), returns `default_block(J)` instead.
"""
function admixture_block(chroms::AbstractVector, positions::AbstractVector{<:Integer}; mb=10)
    J = length(chroms)
    length(positions) == J || throw(DimensionMismatch("chroms and positions differ in length"))
    lo, hi = Dict{eltype(chroms), Int}(), Dict{eltype(chroms), Int}()
    for (c, p) in zip(chroms, positions)
        lo[c] = min(get(lo, c, p), p)
        hi[c] = max(get(hi, c, p), p)
    end
    span = sum(hi[c] - lo[c] + 1 for c in keys(hi))
    block = ceil(Int, J * (mb * 1e6) / span)
    J >= 2block ? block : default_block(J)
end

"""
    block_bootstrap_indices(chroms, block, rng)

SNP indices of one moving-block bootstrap replicate. Each run of consecutive
SNPs on the same chromosome keeps its size and is refilled with blocks of
`block` consecutive SNPs drawn from that same run. `block = 1` resamples SNPs
independently within chromosomes.
"""
function block_bootstrap_indices(chroms::AbstractVector, block::Integer, rng::AbstractRNG)
    block >= 1 || throw(ArgumentError("block length must be ≥ 1"))
    J = length(chroms)
    idx = Int[]
    sizehint!(idx, J)
    start = 1
    while start <= J
        stop = start
        while stop < J && chroms[stop + 1] == chroms[start]
            stop += 1
        end
        b = min(block, stop - start + 1)
        while length(idx) < stop
            s = rand(rng, start:(stop - b + 1))
            for t in s:(s + b - 1)
                length(idx) < stop || break
                push!(idx, t)
            end
        end
        start = stop + 1
    end
    idx
end

"""
    bootstrap_se(g, K, q0, p0; B=200, chroms, block, rng, minor_iter=3,
                 minor_atol=nothing, quiet=false, fit_kwargs...)

Bootstrap standard errors of Q, as an I x K matrix (the layout of a `.Q` file).
`q0` (K x I) and `p0` (K x J) are the point estimates. A replicate stops after
`minor_iter` iterations, or when the log-likelihood changes by less than
`minor_atol` if `minor_iter` is `nothing`.
"""
function bootstrap_se(g::SnpArray, K::Integer, q0::AbstractMatrix, p0::AbstractMatrix;
                      B::Integer=200, chroms::AbstractVector=fill(1, size(g, 2)),
                      block::Integer=default_block(size(g, 2)),
                      rng::AbstractRNG=StableRNG(43),
                      minor_iter::Union{Nothing, Integer}=3,
                      minor_atol::Union{Nothing, Real}=nothing,
                      quiet::Bool=false, fit_kwargs...)
    B >= 2 || throw(ArgumentError("bootstrap needs at least 2 replicates"))
    I, J = size(g)
    length(chroms) == J || throw(ArgumentError("one chromosome per SNP is required"))
    n_iter = minor_iter === nothing ? 1000 : Int(minor_iter)
    # A fixed iteration count must not be cut short by the default relative tolerance.
    rtol = minor_iter === nothing ? 1e-7 : 0.0
    s1 = zeros(Float64, K, I)
    s2 = zeros(Float64, K, I)
    gb = SnpArray(undef, I, J)
    for _ in 1:B
        idx = block_bootstrap_indices(chroms, block, rng)
        for (t, j) in enumerate(idx)
            copyto!(view(gb.data, :, t), view(g.data, :, j))
        end
        fill!(gb.columncounts, 0)  # zero counts are recomputed on demand
        fill!(gb.rowcounts, 0)
        fitrep() = _fit(gb, K; init_q=copy(q0), init_p=p0[:, idx],
            n_iter=n_iter, rtol=rtol, atol=minor_atol,
            rng=StableRNG(rand(rng, UInt32)), progress_bar=false, fit_kwargs...)
        d = quiet ? (@suppress_out fitrep()) : fitrep()
        s1 .+= d.q
        s2 .+= d.q .^ 2
    end
    m = s1 ./ B
    v = max.((s2 .- B .* m .^ 2) ./ (B - 1), 0.0)
    permutedims(sqrt.(v))
end
