"""
    num_vectors(N::Int, L::Int)
Compute the number of vectors in a basis with `N` particles and `L` sites.
"""
num_vectors(N::Int, L::Int) = binomial(N+L-1, L-1)
num_vectors(::Intbasis, N::Int, L::Int) = num_vectors(N, L)


# Global cache of generalized Pascal's triangles for restricted num_vectors.
const triangles = Dict{Int, Matrix{Int}}()

"""
    num_vectors(N::Int, L::Int, n_max::Int)
Compute the number of vectors in a basis with `N` particles and `L` sites, and
a limit of `n_max` particles per site.
"""
function num_vectors(N::Int, L::Int, n_max::Int)
    0 <= N <= n_max * L || return 0
    N == L == 0 && return 1

    # Create a new triangle.
    if !haskey(triangles, n_max)
        triangles[n_max] = zeros(Int, 1, 1)
        triangles[n_max][1, 1] = 1
    end

    # Extend an existing triangle.
    if size(triangles[n_max], 1) < L + 1
        t_old = triangles[n_max]
        L_old = size(t_old, 1) - 1
        t = zeros(Int, L+1, n_max*L+1)
        for k in 0:L_old
            for n in 0:(n_max*k)
                t[k+1, n+1] = t_old[k+1, n+1]
            end
        end
        for k in (L_old+1):L
            for n in 0:(n_max*k)
                for m in 0:min(n_max, n)
                    t[k+1, n+1] += t[k, n+1-m]
                end
            end
        end
        triangles[n_max] = t
    end

    triangles[n_max][L+1, N+1]
end

num_vectors(basis::RestrictedIntbasis, N::Int, L::Int) = num_vectors(N, L, basis.n_max)

function getNvector(v::Int,L::Int64,N::Int64)
	Nvector = zeros(Int,L)
        site=1
        for i=1:L+N-1
           Ni=v&1
           site +=Ni
           Nvector[site] +=~Ni&1
           v >>>= 1
	end
	Nvector
end




















"""
	serial_num(L::Int64, N::Int64, v::Int64)

Compute the serial number of occupation vector `v` in a basis with `L` sites
and `N` particles.
"""

#function serial_num(L::Int, N::Int, v::Int)
#    I = 1
#
#	for mu = 1:L
#		s = 0
#		for nu = (mu+1):L
#			s += countParticles(v, nu)
#		end
#		for i = 0:(countParticles(v, mu) - 1)
#			I += num_vectors(N-s-i, mu-1)
#		end
#	end
#
#	I
#end
#serial_num(basis::Intbasis, v::Int64) = serial_num(basis.L, basis.N, v)




function serial_num_fast(basis:: AbstractIntbasis, v::Int)
        Max_Index=basis.D+1
        Min_Index=1
        is_not_found=true
        Index=1
        while is_not_found
            Index= (Max_Index+ Min_Index)÷ 2
            if basis.vectors[Index]>v
                Min_Index=Index
            elseif basis.vectors[Index]<v
                Max_Index=Index
            else 
                is_not_found=false
		Index=-1
            end
        end
	Index
end







"""
	checkSite(v::Int, U::Array, site::Int64)

Return 1 if a given site is occupied, 0 if it is unoccupied.

Shifts the digit right of the index of the site to the least significant bit,
and uses `& 1` to check the state. If `site = 1` and the site is empty,
`U[site] - 2` equals negative one, right shifting by a negative number. This is
the same as a left shift by one, which introduces a new leading zero. `~v` instructs
the function to interpret zeros as empty sites.
"""
function checkSite(v::Int, U::Array, site::Int)::Int64
    (~v >>> (U[site] - 2)) & 1
end

"""
    hopLeft(v::Int, U::Array, site::Int64)

Hop a particle to the left with respect to the binary
representation, to the right with respect to lattice indices.

A swap of the leftmost particle in a site and the bar adjacent by using
`⊻ 0b11` to swap each bit. In the edge case, the leftmost particle is
"lost" in the leading zeros, so the integer is shifted left once.
"""
function hopLeft(v::Int, U::Array, site::Int)
    # checkSite(v, U, site) == 1 || error("No particle found on site ", site) # ! removed because it ran so often it used a significant amount of memory
    U[site] -= 1
    newV = v ⊻ (3 << (U[site] - 1))
    _zeros = leading_zeros(newV) - leading_zeros(v)
	@inbounds for i = 1:size(U,1)
        U[i] += _zeros
    end
    newV <<= _zeros
    newV
end

"""
    emptyToFirst(v::Int, U::Array, site::Int64)

Move all particles from a given site to the first site.

`v & (1 << U[site-1] - 1)` gives the value of the digits to the right of
the occupied site. `v >>> (U[site] - 1)` gives the digits to the left,
which is then shifted left and `|` is used to combine the two values.
The whole integer is shifted left by the number of particles removed.
"""
function emptyToFirst(v::Int, U::Array, site::Int)
    # checkSite(v, U, site) == 1 || error("No particle found on site ", site) # Got rid of because it was using memory also
    @inbounds for i = 1:site-1
        U[i] += U[site] - U[site-1] - 1
    end
    (v & (1 << U[site-1] - 1) | v >>> (U[site] - 1) << (U[site-1] - U[1] + 1)) << (U[1] - 1)
end

"""
    findFirstOccupied(v::Int)

Return the location of the rightmost occupied site.

`trailing_ones` returns the number of bars preceding the rightmost particle,
and `+ 1` returns the position of the adjacent site to the left.
"""
function findFirstOccupied(v::Int)
	trailing_ones(v) + 1
end

"""
   countParticles(v::Int, site::Int64)

Count the number of particles on a given site.

Shifts v to the right as many times as needed to have the site at the far
right, and counts the trailing zeros.
"""
function countParticles(v::Int, site::Int64)
    for i = 1:(site-1)
		v >>>= trailing_zeros(v) + 1
	end
	trailing_zeros(v)
end


"""
	SzToInt(sz::Array{Int})

Utility to convert Jeszenszki representation of a basis vector to an integer representation.
"""
function SzToInt(sz::Array{Int})
	v = 0
	sum = -1
	@inbounds for i = 1:length(sz)
		sum += sz[i]
		v |= 1 << (sum + i)
	end
	v
end


function printV(v::Int)
	L = count_ones(v)
	U = zeros(Int, L)
	index = 0
	@inbounds for i = 1:L
		index += trailing_zeros(v) + 1
		U[i] = index
		v >>>= (trailing_zeros(v) + 1)
	end

	print("┌─", "──"^(U[1] - 1))
	for i = 2:L
		print("┬─", "──"^(U[i] - U[i-1] - 1))
	end
	print("┐\n│ ", "✪ "^(U[1] - 1), "│ ")
	for i = 2:L
		print("✪ "^(U[i] - U[i-1] - 1), "│ ")
	end
	print("\n└─", "──"^(U[1] - 1))
	for i = 2:L
		print("┴─", "──"^(U[i] - U[i-1] - 1))
	end
	println("┘")
	
	nothing
end
