
"""
 set a bit to 1  
"""

function SetToOne(v:: Int64,bit::Int64)
    return ((1 << (bit-1)) | v)
end
"""
 Reverse k bits from right to lift (101100 --> 001101)
 k=N+L
"""

function ReverseKet(v:: Int64,k:: Int64)
           v=Mybswap(v) >>> (64-k+1)
    return SetToOne(v:: Int64,k::Int64)
end
ReverseKet(v::Int64,basis:: AbstractIntbasis) = ReverseKet(v,basis.L+basis.N)


function Mybswap(i::Int64)
           i=(i >>> 32)|(i << 32)
           i=((i& -281470681808896)>>>16)|((i&281470681808895)<<16)
           i=((i& -71777214294589696)>>>8)|((i&71777214294589695)<<8)
           i=((i& -1085102592571150096)>>>4)|((i&1085102592571150095)<<4)
           i=((i& -3689348814741910324)>>>2)|((i&3689348814741910323)<<2)
           return (((i& -6148914691236517206)>>>1)|((i&6148914691236517205)<<1))
end

function Mybswap(i::UInt64)
           i=(i >>> 32)|(i << 32)
           i=((i&0xffff0000ffff0000)>>>16)|((i&0x0000ffff0000ffff)<<16)
           i=((i&0xff00ff00ff00ff00)>>>8)|((i&0x00ff00ff00ff00ff)<<8)
           i=((i&0xf0f0f0f0f0f0f0f0)>>>4)|((i&0x0f0f0f0f0f0f0f0f)<<4)
           i=((i&0xcccccccccccccccc)>>>2)|((i&0x3333333333333333)<<2)
           return (((i&0xaaaaaaaaaaaaaaaa)>>>1)|((i&0x5555555555555555)<<1))
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
    someFunction!(U,_zeros)
    newV <<= _zeros
    newV
end

"""
    hopRight(v::Int, U::Array, site::Int64)

Hop a particle to the right with respect to the binary
representation, to the left with respect to lattice indices.

A swap of the rightmost particle in a site and the bar adjacent by using
`⊻ 0b11` to swap each bit. The site of the bar to flip is equal to `site - 1`,
but if `site = 1` the site of the bar to flip is the leftmost one. In
this case, the integer is shifted right once to remove the first particle.
"""
function hopRight(v::Int, U::Array, site::Int)
    checkSite(v, U, site) == 1 || error("No particle found on site ", site)
    site = length(U) - (length(U) - site + 1) % length(U)
    U[site] += 1
    newV = v ⊻ (3 << (U[site] - 2))
    # someFunction!(U, newV, leading_zeros(newV) - leading_zeros(v))
    for i = 1:size(U,1)
        @inbounds U[i] += leading_zeros(newV) - leading_zeros(v)
    end
    newV <<= leading_zeros(newV) - leading_zeros(v)
end

"""
    someFunction!(U::Array,newV::Int,zeros::Int)

Utility pulled out of hopRight and hopLeft.
"""
function someFunction!(U::Array{Int,1},_zeros::Int)
    #U .+= leading_zeros(newV) - leading_zeros(v)
    @inbounds for i = 1:size(U,1)
        U[i] += _zeros;
    end
    nothing
end


"""
	destroyParticle(v::Int, U::Array{Int64,N} where N, site::Int64)

Destroy a particle on a given site.

`v & (1 << (U[site] - 1) - 1)` gives the digits of v to the right of the particle
to be destroyed. `v >>> U[site] << (U[site] - 1)` gives the digits to the right,
left shifted to the desired position. `|` combines the two, and U is updated to
decrement all indices to the left of the destroyed particle.
"""
function destroyParticle(v::Int, U::Array{Int}, site::Int)
	#checkSite(v, U, site) == 1 || error("No particle found on site ", site)

	@inbounds for i = site:length(U)
		U[i] -= 1
	end
	(v & (1 << (U[site] - 1) - 1)) | v >>> U[site] << (U[site] - 1)
end

"""
	createParticle(v::Int, U::Array{Int64,N} where N, site::Int64)

Create a particle on a given site.

`v & (1 << (U[site] - 2) - 1)` gives the digits to the right of the particle,
including the particle. `v >>> (U[site] - 2) << (U[site] - 1)` gives the digits
to the left, left shifted one digit past the right half. `|` combines the two,
and U is updated to increment all indices left of the created particle.
"""
function createParticle(v::Int, U::Array{Int}, site::Int)
	@inbounds for i = site:length(U)
		U[i] += 1
	end
	(v & (1 << (U[site] - 2) - 1))| v >>> (U[site] - 2) << (U[site] - 1)
end

"""
Used in place of Create and Destroy to move a particle from site 1 to site 2. Does not update U. It's super long, but just a combination of he operations in createParticle and destroyParticle.
"""
function TunnelParticle(v::Int, U::Array{Int}, site1::Int, site2::Int)
    @inbounds for i = site1:length(U)
        U[i] -= 1
    end
    v = (v & (1 << (U[site1] - 1) - 1)) | v >>> U[site1] << (U[site1] - 1)

	(v & (1 << (U[site2] - 1) - 1))| v >>> (U[site2] - 1) << U[site2]
end

# The function occupationShift might exist at some point, which would
# be a more generalized version of this. Doesn't change U, bc whatever. Not well documented.
"Shift the function to the right once, used in translational_symmetry.jl"
#function shiftVector(v::Int, U::Array)
#    oldU = copy(U)
#    for i = 1:(oldU[1]-1)
#        v = hopLeft(v, U, 1)
#    end
#    for i = 2:length(oldU)
#        for j = 1:(oldU[i]-oldU[i-1]-1)
#            v = hopLeft(v, U, i)
#        end
#    end
#    v
#end
#shiftVector(v::Int)  = shiftVector(v, generateU(v))
#function shiftVector(v::Int, i::Int)
#    U = generateU(v)
#    for i = 1:i
#        v = shiftVector(v, U)
#    end
#    v
#end


function shiftV(v::Int,k::Int64)
    v >>>= trailing_zeros(v) + 1
    return SetToOne(v,k)
end

#function shiftV(v::Int,L::Int64,N::Int64)
#    v=~v
#    one=1
#    for bit=1:N
#       one= v&one
#       v >>>= one       
#    end
#    v=~v
#    v >>>= 1
#    return SetToOne(v,L+N)
#end

#function shiftV(v::Int,M::Int64,k::Int64)
#    v >>>= (k + 1)
#    return SetToOne(v,M)
#end


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
            end
        end
	Index
end

function serial_num_Cycles(Cycles_leaders::Array{Int,1}, NumOfCycles::Int, v::Int)
        Max_Index= NumOfCycles+1
        Min_Index=1
        is_not_found=true
        Index=1
        while is_not_found
            Index= (Max_Index+ Min_Index)÷ 2
            if Cycles_leaders[Index]>v
                Min_Index=Index
            elseif Cycles_leaders[Index]<v
                Max_Index=Index
            else 
                is_not_found=false
            end
        end
	Index
end




function serial_num_Cycles_OTF(Cycles_leaders::Array{Int,1}, NumOfCycles::Int, v::Int)
        Max_Index= NumOfCycles+1
        Min_Index=1
        is_not_found=true
        if  v==1
                is_not_found=false

        end 

        Index=1
        while is_not_found
            Index= (Max_Index+ Min_Index)÷ 2
            if Cycles_leaders[Index]<v
                Min_Index=Index
            elseif Cycles_leaders[Index]>v
                Max_Index=Index
            else 
                is_not_found=false
            end
        end
	Index
end






"""
	reverseBasis(v::Int, U::Array{Int64,N} where N, L::Int64)

Reverse a given basis vector.

newU is calculated by placing the indices at their location relative to the
leftmost bar, rather than the least significant bit. v is then generated by
shifting ones to those new locations.
"""
function reverseBasis(v::Int, U::Array{Int}, L::Int)
	newU = Array{Int}(undef, L)
	newU[L] = U[L]
	for i = 1:(L-1)
		newU[i] = U[L] - U[L - i]
	end

	v = 0
	for i = 1:L
		v |= 1 << (newU[i] - 1)
	end
	U = newU
	v
end

"""
    generateU!(U::Array{Int64,N} where N, v::Int, L::Int64)

Generate U for a given vector representation.

Cumulatively counts the number of digits in the represenation, storing
the locations of the ones by repeatedly counting leading zeros and
shifting the integer past its location. Returns nothing but modifies U
in memory.

# Arguments
- `U::Array{Int}`: an arbitrary array of length L.
- `v::I `: a stars-and-bars basis vector.
- `L::Int`: the number of sites.
"""
function generateU!(U::Array{Int}, v::Int, L::Int)
	index = 0
	@inbounds for i = 1:L
		index += trailing_zeros(v) + 1
		U[i] = index
		v >>>= (trailing_zeros(v) + 1)
	end
	nothing
end

"""
	generateU(v::Int)

Generate U from a given basis vector.

An alternative function for generateU!, this function creates U in memory.
"""
function generateU(v::Int)
	U = zeros(Int, count_ones(v))
	generateU!(U, v, count_ones(v))
	U
end

function generateN(U::Array{Int})
	N = copy(U)
	for i = length(N):-1:2
		N[i] = N[i] - N[i-1] - 1
	end
	N[1] = N[1] - 1
	N
end

#function getNvector(v::Int,L::Int64,N::Int64)
#	Nvector = zeros(Int,L)
#        site=1
#        for i=1:L+N-1
#           Ni=v&1
#           site +=Ni
#           Nvector[site] +=~Ni&1
#           v >>>= 1
#	end
#	Nvector
#end

"""
	countN(v::Int)

Count total number of particles in vector `v`.
"""
countN(v::Int)  = count_zeros(v) - leading_zeros(v)


# Flip k bits (0 <--> 1) (101100 --> 010011)
function FlipKet(v::Int, L:: Int64)
    ~v << (64-L) >>> (64-L)
end

"""
	subregion(v::Int, A::Array)

Get a subregion of basis vector `v` of consecutive sites of indices A.

Similar to countParticles.
"""
function subregion(v::Int, A::Array)
	# Right shift to the first zero in in the first site in A
	for i = 1:A[1]-1
		v >>>= (trailing_zeros(v) + 1)
	end
	x = trailing_zeros(v)
	# left shift to the closing bar of the last site in A
	for i = 1:count_ones(v)-length(A)
		v <<= (leading_zeros(v) + 1)
	end
	v <<= leading_zeros(v)
	v  >>> trailing_zeros(v) << x
end
# may eventually want to do this for nonconsecutive sites




#original from A. V. Ponomarev (2009)
#subroutine buildpascal
function buildpascal(lc::Int,nc::Int)
#c lc=number of sites +1
#c nc=number of atoms +1
    jbc=Float64   
    cnkc=zeros(Int,lc,nc)    
    jmaxnc=Int
#c builds the rotated pascal triangle

    for i = 1:lc
        cnkc[i,1] = 1
    end
    for i = 1:lc
        for j = 2:nc
            cnkc[i,j] = 0
        end 
    end 
    for in1 = 2:lc
        cnkc[in1,2] = sum(cnkc[in1-1,1:2])
        if nc-1>1 
            for in2 = 1:nc
                cnkc[in1,in2] = sum(cnkc[in1-1,1:in2])
            end 
        end
    end
    jmax = cnkc[lc,nc]
    return cnkc,jmax
end

#c ---------------------------------------------
#c Returns the many body state bi at position in
#c ---------------------------------------------
#c original from A. V. Ponomarev (2009)
#subroutine in2b(in,bi)
function in2b(inn::Int64,cnkc::Array{Int64,2},jmax::Int64,lc::Int64,nc::Int64,bi::Vector{Int64})
    bi.=0 
    indi = inn-1
    ind_L = lc-1
    ind_N = nc
    while ind_N != 1
        if indi >= cnkc[ind_L,ind_N] 
            indi -= cnkc[ind_L,ind_N]
            bi[lc-ind_L] += 1
            ind_N -= 1
        else
            ind_L -= 1
        end
    end       
end
#c ---------------------------------------------
#c Returns the many body state bi at position in
#c ---------------------------------------------
#c original from A. V. Ponomarev (2009)
#subroutine b2in(bi,in)
function b2in(bi::Vector{Int64},cnkc::Array{Int64,2},jmax::Int64,lc::Int64,nc::Int64)
    inn=1
    for indi=1:lc-2
        for ind_N=0:bi[indi]
            if (bi[indi]-ind_N>0) 
                suma=0
                for k=1:indi-1
                    suma=suma+bi[k]
                end
                if (lc-indi>0) && (nc-ind_N-suma>0) 
                    is=0
                    inn+=cnkc[lc-indi,nc-ind_N-suma]
                end
            end
        end
    end
    return inn
end

function Get_ket_leader(ket::Int64,L::Int64,M::Int64)
    #ket_leader::Int64
    CycleMembs=zeros(Int64, L*2)
    CycleMembs[1]= ket
    ket_rev=ReverseKet(ket,M) #M=L+N
    CycleMembs[L+1]= ket_rev
    for  k= 2:L
        ket = shiftV(ket,M)
        ket_rev = shiftV(ket_rev,M)
        CycleMembs[k]= ket
        CycleMembs[L+k]= ket_rev
    end
    ket_leader = maximum(CycleMembs)
    return ket_leader
end

function Get_ket_leader_OTF(Nv::Vector{Int64},L::Int64,N::Int64,jmax::Int64,cnkc::Array{Int,2})
    lc=L+1
    nc=N+1
    CycleMembs=zeros(Int64, L*2)
    ket =b2in(Nv,cnkc,jmax,lc,nc)
    CycleMembs[1]= ket
    Nv_rev=copy(Nv)
    Nv_temp=copy(Nv)
    reverse!(Nv_rev)
    ket_rev = b2in(Nv_rev,cnkc,jmax,lc,nc)
    CycleMembs[L+1]= ket_rev
    for  k= 2:L
        circshift!(Nv_temp, Nv,1)
        Nv_temp, Nv = Nv,Nv_temp
        ket = b2in(Nv,cnkc,jmax,lc,nc)
        CycleMembs[k]= ket
        circshift!(Nv_temp, Nv_rev,1)
        Nv_temp, Nv_rev = Nv_rev,Nv_temp
        ket_rev = b2in(Nv_rev,cnkc,jmax,lc,nc)
        CycleMembs[L+k]= ket_rev
    end
    ket_leader = minimum(CycleMembs)
    return ket_leader
end

