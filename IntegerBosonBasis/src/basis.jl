abstract type AbstractIntbasis end

"""
Basis of Bosonic occupation vectors in stars and bars format.
"""
struct Intbasis <: AbstractIntbasis
    "Number of sites."
    L::Int
    "Number of bosons."
    N::Int

    "Number of basis vectors."
    D::Int

    "Occupation vectors (1 by D)."
    vectors::Vector{Int}
end

"""
	Intbasis(L::Int64, N::Int64)

Create a basis for `L` sites and `N` particles.
"""

function Intbasis(L::Int, N::Int)
    L >= 1 || throw(DomainError(L, "At least 1 site is required."))
    N >= 0 || throw(DomainError(N, "At least 0 particles are required."))

    # Basis size.
    D = num_vectors(N,L)
    vectors = Vector{Int}(undef, D)

	# Create initial values
    v = 2^(L+N) - 2^N
    U = collect(Int, N+1:N+L)
    vectors[1] = v

	# Loop to fill basis vectors
	@inbounds for i in 2:D
	    if checkSite(v, U, 1) == 1
			v = hopLeft(v,U,1)
	    else
	        j = findFirstOccupied(v)
			v = hopLeft(v, U, j)
			if checkSite(v, U, j) == 1
				v = emptyToFirst(v, U, j)
			end
	    end
		vectors[i] = v
	end

    Intbasis(L, N, D, vectors)
end


"""
Basis of occupation vectors with a site occupation restriction.
"""
struct RestrictedIntbasis <: AbstractIntbasis
    "Number of sites."
    L::Int
    "Number of particles."
    N::Int
    "Site capacity."
    n_max::Int

    "Number of basis vectors."
    D::Int

    "Occupation vectors (1 by D)."
    vectors::Vector{Int}
end

"""
    RestrictedIntbasis(L::Int, N::Int, n_max::Int)
Create a basis for `L` sites and `N` particles, with no more than `n_max` particles
per site.
"""
function RestrictedIntbasis(L::Int, N::Int, n_max::Int)
    L >= 1 || throw(DomainError(L, "At least 1 site is required."))
    N >= 0 || throw(DomainError(N, "At least 0 particles are required."))
    N <= L * n_max || throw(DomainError(N, "Particles do not fit on the sites."))
    n_max !=1 || @warn("Inefficient for hard-core bosons. Use (to be added later) instead.")

    # Basis size.
    D = num_vectors(N, L, n_max)
    dNn_max = n_max > 0 ? div(N,n_max) : 1
    rNn_max = n_max > 0 ? N-dNn_max*n_max : 0 
    vectors = Vector{Int}(undef, D)

    # Create initial values
    v=0
    U = zeros(Int,L)
    Num= zeros(Int,L)
    for j in 1:dNn_max
        v+= 2^(n_max*j+j-1)
        U[j]=n_max*j+j
        Num[j]=n_max
    end

    if 1 <= (dNn_max + 1) <= L
        v+=2^(N+dNn_max)
        U[dNn_max+1] = N+dNn_max+1
        Num[dNn_max+1]=rNn_max
    end
    if dNn_max+1<L
        v+= 2^(L+N) - 2^(N+dNn_max+1)
    end
    for j in dNn_max+2:L
        U[j]=N+dNn_max+1+j
    end
    vectors[1] = v

    # Loop to fill basis vectors
    @inbounds for i in 2:D
        if Num[1] > 0
            if Num[1] < n_max
                delta = n_max - Num[1]
                Num[1] = n_max
            else
                delta = 0
            end

            j = findfirst(!iszero, Num .< n_max)

            Num[j] += 1
            Num[j-1] -= 1 + delta
        else
            j = findfirst(!iszero, Num)
            k = j + findfirst(!iszero, @view(Num[(j+1):end]) .< n_max)

            Num[k-j] = Num[j] - 1
            Num[k] += 1
            for l in 1:(k-j-1)
                Num[l] = n_max
            end
            # The indices after the first one differ from those in the paper.
            for l in (k-j+1):(k-1)
                Num[l] = 0
            end
        end
        v=SzToInt(Num)
	vectors[i] = v
    end

    RestrictedIntbasis(L, N, n_max, D, vectors)
end
