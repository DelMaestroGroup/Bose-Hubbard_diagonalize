"""
Calculate the particle entanglement entropy for a subset A, using the SVD.
"""
function particle_entropy(basis:: AbstractIntbasis, Asize::Int, d::Vector{T}, CycleSize::Array{Int64,1}, InvCycles_Id::Vector{Int64}) where {T<:Number}
    basisA, basisB = particle_entropy_bases(basis, Asize)
    # Matrix to SVD
    Amatrix = zeros(T, length(basisA), length(basisB))
    #for i = 1:length(d)
    #    d[i] /= sqrt(CycleSize[i])
    #end
    for (i, braA) in enumerate(basisA)
        for (j, braB) in enumerate(basisB)
            braAVector= getNvector(braA, basisA.L, basisA.N)
            braBVector= getNvector(braB, basisB.L, basisB.N)
            braVector = braAVector + braBVector
            bra= SzToInt(braVector)
            #bra in basis || continue
            sn= serial_num_fast(basis, bra)
            if sn>0
                norm = 1.0
                for k in 1:basis.L
                    norm *= binomial(braVector[k], braAVector[k])
                end
                Amatrix[i, j] = sqrt(norm) * d[InvCycles_Id[sn]] / sqrt(binomial(basis.N, basisA.N))
            end
        end
    end

    S = svdvals(Amatrix)
    err = abs(sum(S.^2) - 1.0)
    S1=0.0
    for k=1:length(S)
        if abs(S[k])>0
            S1 -= abs(S[k])^2*log(abs(S[k])^2)
        end
    end
    S2=-log(sum(S.^4))
    if err > 1e-12
        @warn("RDM eigenvalue error: $(err)")
    end

    S1,S2
end


"""
Generate the bases used for the particle entanglement entropy.
"""
function particle_entropy_bases(basis::Intbasis, Asize::Int)
    Bsize = basis.N - Asize
    basisA = Intbasis(basis.L, Asize)
    basisB = Intbasis(basis.L, Bsize)

    basisA, basisB
end

function particle_entropy_bases(basis::RestrictedIntbasis, Asize::Int)
    Bsize = basis.N - Asize
    if Asize <= basis.n_max
        basisA = Intbasis(basis.L, Asize)
    else
        basisA = RestrictedIntbasis(basis.L, Asize, basis.n_max)
    end
    if Bsize <= basis.n_max
        basisB = Intbasis(basis.L, Bsize)
    else
        basisB = RestrictedIntbasis(basis.L, Bsize, basis.n_max)
    end

    basisA, basisB
end
