"""
Create a list of occupation basis for each translational and reflection symmetry cycle for bosonic 1D chains with PBC.
"""
function Translational_Symmetry_Reflection_Cycles_gs(basis::AbstractIntbasis)

    # Makes a Symbol for each, of type Int64
    for x = [:NumOfCycles,:Num_cycles_max,:MemberID,:i_next]
       @eval $x = Int64
    end

    IdStatus = Bool
    D = basis.D  #length(basis)
    L= basis.L
    N= basis.N
    M=L+N
    # Approximate the maximum number of cycles in the basis
    Num_cycles_max = round(Int64,((D/L/2)+6*(D/L/2)^.5))
    # Make an empty array with a row for each cycle and L (# of sites) rows
    Cycles_leaders = zeros(Int64, Num_cycles_max)
    # Make an empty array with number_of_cycles values
    CycleSize = zeros(Int64, Num_cycles_max)
    InvCycles_Id = zeros(Int64, D)
    # Give every vector in the basis a state of true, it's not in a cycle yet
    Status  = trues(D)
    NumOfCycles = 0
    MemberID=0
    # For each vector in the basis
        for i=1: D     
        if Status[i]
            # starting with bra
            i_next = i
            j=2
            while j >0
                # If the vector is not yet in a cycle
                if Status[i_next]
                    if j==2
                        # Increment the number of cycles
                        NumOfCycles += 1
                        # Make bra the first member of a cycle
                        MemberID=0
                        Cycles_leaders[NumOfCycles]= i_next

                    end
                    # IdStatus indicates if the "current" vector is in a cycle
                    IdStatus=true
                    # i_next is the serial number of the "following" vector,
                    # While the "following" vector is not yet in a cycle
                    while IdStatus
                        # The "following" vector is now in a cycle
                        Status[i_next]=false
                        MemberID+=1
                        # Add the "following" vector to Cycles

                        InvCycles_Id[i_next] = NumOfCycles
  
                        # i_next is now the serial number of the vector rotated once
                        i_next = serial_num_fast(basis, shiftV(basis.vectors[i_next],M))
                        IdStatus = Status[i_next]
                    end
                else
                    j=0
                end
                j-=1
                if j==1
                    i_next=serial_num_fast(basis,ReverseKet(basis.vectors[i],M))
                    if ~Status[i_next]
                        j=0
                    end                
                end
                if j==0
                    CycleSize[NumOfCycles]= MemberID
                end
            end
        end
    end
    Cycles_leaders, CycleSize, NumOfCycles, InvCycles_Id
end


"""
Create the sparse translational (and reflection) symmetry block of the hamiltonian of bosonic 1D chains with PBC.
"""
function sparse_Block_Diagonal_Hamiltonian_q0R1_gs(basis::AbstractIntbasis, Cycles_leaders:: Vector{Int64}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, InvCycles_Id:: Vector{Int64}, t::Float64, U::Float64)

    ### added to get to work like in reg sparse
    mus = zeros(basis.L)

    #Creating the block H_(q=0,R=1) of the hamiltonian.
    rows = Int64[]
    cols = Int64[]
    elements = Float64[]
    factor=0.0
    for CycleId = 1:NumOfCycles
        # Diagonal part
        Usum = 0.0
        musum = 0.0
        bra = basis.vectors[Cycles_leaders[CycleId]]
        U_bra = generateU(bra)
        N = generateN(U_bra)
        for j in 1:basis.L
            Usum += N[j] * (N[j] - 1)
            musum += mus[j] * N[j]
        end
        push!(rows, CycleId)
        push!(cols, CycleId)
        push!(elements, U*Usum/2.0 - musum)

        # Off-diagonal part
        CycleId_CycleSize=CycleSize[CycleId]
        for j = 1:basis.L
             j_next = j % basis.L + 1
             # Tunnel right, tunnel left.
             for (site1, site2) in [(j, j_next), (j_next, j)]
                 if N[site1] > 0
                     ket = bra
                     U_bra = generateU(bra)
                     ket = TunnelParticle(ket, U_bra, site1,site2)
                     kId = serial_num_fast(basis, ket)
                     CycleId1 = InvCycles_Id[kId]
                     CycleId1_CycleSize = CycleSize[CycleId1]
                     k1 = CycleId1
                     if CycleId >= k1
                         factor = - t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)*1.0 * sqrt(N[site1]) * sqrt(N[site2]+1)
                         push!(rows, k1)
        	         push!(cols, CycleId)
                         push!(elements, factor)
                         #push!(rows, CycleId)
                         #push!(cols, k1)
                         #push!(elements, factor)
                     end
                 end
             end
         end
    end
        return sparse(Symmetric(sparse(rows, cols, elements, NumOfCycles, NumOfCycles)))
end



"""
Using RestrictedIntbasis
Create the sparse translational (and reflection) symmetry block of the hamiltonian of bosonic 1D chains with PBC.
"""
function sparse_Block_Diagonal_Hamiltonian_q0R1_gs_Rest(basis::AbstractIntbasis, Cycles_leaders:: Vector{Int64}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, InvCycles_Id:: Vector{Int64}, t::Float64, U::Float64)

    ### added to get to work like in reg sparse
    mus = zeros(basis.L)

    #Creating the block H_(q=0,R=1) of the hamiltonian.
    rows = Int64[]
    cols = Int64[]
    elements = Float64[]
    factor=0.0
    for CycleId = 1:NumOfCycles
        # Diagonal part
        Usum = 0.0
        musum = 0.0
        bra = basis.vectors[Cycles_leaders[CycleId]]
        U_bra = generateU(bra)
        N = generateN(U_bra)
        for j in 1:basis.L
            Usum += N[j] * (N[j] - 1)
            musum += mus[j] * N[j]
        end
        push!(rows, CycleId)
        push!(cols, CycleId)
        push!(elements, U*Usum/2.0 - musum)
        n_max=basis.n_max
        # Off-diagonal part
        CycleId_CycleSize=CycleSize[CycleId]
        for j = 1:basis.L
             j_next = j % basis.L + 1
             # Tunnel right, tunnel left.
             for (site1, site2) in [(j, j_next), (j_next, j)]
                 if (N[site1] > 0) && (N[site2] < n_max)
                     ket = bra
                     U_bra = generateU(bra)
                     ket = TunnelParticle(ket, U_bra, site1,site2)
                     kId = serial_num_fast(basis, ket)
                     CycleId1 = InvCycles_Id[kId]
                     CycleId1_CycleSize = CycleSize[CycleId1]
                     k1 = CycleId1
                     if CycleId >= k1
                         factor = - t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)*1.0 * sqrt(N[site1]) * sqrt(N[site2]+1)
                         push!(rows, k1)
        	         push!(cols, CycleId)
                         push!(elements, factor)
                         #push!(rows, CycleId)
                         #push!(cols, k1)
                         #push!(elements, factor)
                     end
                 end
             end
         end
    end
        return sparse(Symmetric(sparse(rows, cols, elements, NumOfCycles, NumOfCycles)))
end

# ===========================================reduced=========================================

"""
Create a list of occupation basis for each translational and reflection symmetry cycle for bosonic 1D chains with PBC.
"""
function Translational_Symmetry_Reflection_Cycles_gs_reduced(basis::AbstractIntbasis)

    # Makes a Symbol for each, of type Int64
    for x = [:NumOfCycles,:Num_cycles_max,:MemberID,:i_next]
       @eval $x = Int64
    end

    IdStatus = Bool
    D = basis.D  #length(basis)
    L= basis.L
    N= basis.N
    M=L+N
    # Approximate the maximum number of cycles in the basis
    Num_cycles_max = round(Int64,((D/L/2)+6*(D/L/2)^.5))
    # Make an empty array with a row for each cycle and L (# of sites) rows
    Cycles_leaders = zeros(Int64, Num_cycles_max)
    # Make an empty array with number_of_cycles values
    CycleSize = zeros(Int64, Num_cycles_max)
    # Give every vector in the basis a state of true, it's not in a cycle yet
    Status  = trues(D)
    Reflection_Extension  = falses(Num_cycles_max)
    NumOfCycles = 0
    MemberID=0
    # For each vector in the basis
        for i=1: D     
        if Status[i]
            # starting with bra
            i_next = i
            j=2
            while j >0
                # If the vector is not yet in a cycle
                if Status[i_next]
                    if j==2
                        # Increment the number of cycles
                        NumOfCycles += 1
                        # Make bra the first member of a cycle
                        MemberID=0
                        Cycles_leaders[NumOfCycles]= basis.vectors[i_next]
                    end
                    # IdStatus indicates if the "current" vector is in a cycle
                    IdStatus=true
                    # i_next is the serial number of the "following" vector,
                    # While the "following" vector is not yet in a cycle
                    while IdStatus
                        # The "following" vector is now in a cycle
                        Status[i_next]=false
                        MemberID+=1  
                        # i_next is now the serial number of the vector rotated once
                        i_next = serial_num_fast(basis, shiftV(basis.vectors[i_next],M))
                        IdStatus = Status[i_next]
                    end
                else
                    j=0
                end
                j-=1
                if j==1
                    i_next=serial_num_fast(basis,ReverseKet(basis.vectors[i],M))
                    if ~Status[i_next]
                        j=0
                    else
                        Reflection_Extension[NumOfCycles]=true
                    end                
                end
                if j==0
                    CycleSize[NumOfCycles]= MemberID
                end
            end
        end
    end
    Cycles_leaders, CycleSize, NumOfCycles, Reflection_Extension
end



"""
Create the sparse translational (and reflection) symmetry block of the hamiltonian of bosonic 1D chains with PBC.
"""
function sparse_Block_Diagonal_Hamiltonian_q0R1_gs_reduced(Cycles_leaders:: Vector{Int64}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, t::Float64, U::Float64, L::Int64, N::Int64)

    ### added to get to work like in reg sparse
    mus = zeros(L)
    U_bra = zeros(Int64,L)
    U_temp = zeros(Int64,L)

    Nv = zeros(Int64,L)
    M=L+N
    LookupTable,NumOnes =buildLookupTable()
    #Creating the block H_(q=0,R=1) of the hamiltonian.
    rows = Int64[]
    cols = Int64[]
    elements = Float64[]
    factor=0.0
    for CycleId = 1:NumOfCycles
        # Diagonal part
        Usum = 0.0
        musum = 0.0
        bra =Cycles_leaders[CycleId]
        bra_temp=copy(bra)
        getU(bra_temp,L, LookupTable,NumOnes, U_bra)
        #U_bra = generateU(bra)
        #Nv = generateN(U_bra)
        U_to_Nv(Nv,U_bra,L)
        for j in 1:L
            Usum += Nv[j] * (Nv[j] - 1)
            musum += mus[j] * Nv[j]
        end
        push!(rows, CycleId)
        push!(cols, CycleId)
        push!(elements, U*Usum/2.0 - musum)
	#CycleMembs=zeros(Int64, L*2)
        # Off-diagonal part
        CycleId_CycleSize=CycleSize[CycleId]
        for j = 1:L
             j_next = j % L + 1
             # Tunnel right, tunnel left.
             for (site1, site2) in [(j, j_next), (j_next, j)]
                 if Nv[site1] > 0
                     ket = bra
                     #U_bra = generateU(bra)
                     U_temp=copy(U_bra)
                     ket = TunnelParticle(ket, U_temp, site1,site2)
                     ket_leader = Get_ket_leader(ket,L,M)

#                     CycleMembs[1]= ket
#                     ket_rev=ReverseKet(ket,M)
#                     CycleMembs[L+1]= ket_rev
#                     for  k= 2:L
#                         ket = shiftV(ket,M)
#                         ket_rev = shiftV(ket_rev,M)
#                         CycleMembs[k]= ket
#                         CycleMembs[L+k]= ket_rev
#                     end
#                     ket_leaders = maximum(CycleMembs)

                     CycleId1 = serial_num_Cycles(Cycles_leaders, NumOfCycles, ket_leader)
                     CycleId1_CycleSize = CycleSize[CycleId1]
                     k1 = CycleId1
                     if CycleId >= k1
                         factor = - t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)*1.0 * sqrt(Nv[site1]) * sqrt(Nv[site2]+1)
                         push!(rows, k1)
        	         push!(cols, CycleId)
                         push!(elements, factor)
                         #push!(rows, CycleId)
                         #push!(cols, k1)
                         #push!(elements, factor)
                     end
                 end
             end
         end
    end
        return Symmetric(sparse(rows, cols, elements, NumOfCycles, NumOfCycles))
#        return  sparse(rows, cols, elements, NumOfCycles, NumOfCycles)

end


"""
Create a list of occupation basis for each translational and reflection symmetry cycle for bosonic 1D chains with PBC.
"""
function Translational_Symmetry_Reflection_Cycles_gs_OTF(basis::AbstractIntbasis)

    # Makes a Symbol for each, of type Int64
    for x = [:NumOfCycles,:Num_cycles_max,:MemberID,:i_next]
       @eval $x = Int64
    end

    IdStatus = Bool
    D = basis.D  #length(basis)
    L= basis.L
    N= basis.N
    M=L+N
    # Approximate the maximum number of cycles in the basis
    Num_cycles_max = round(Int64,((D/L/2)+6*(D/L/2)^.5))
    # Make an empty array with a row for each cycle and L (# of sites) rows
    Cycles_leaders = zeros(Int64, Num_cycles_max)
    # Make an empty array with number_of_cycles values
    CycleSize = zeros(Int64, Num_cycles_max)
    # Give every vector in the basis a state of true, it's not in a cycle yet
    Status  = trues(D)
    Reflection_Extension  = falses(Num_cycles_max)
    NumOfCycles = 0
    MemberID=0
    # For each vector in the basis
        for i=1: D     
        if Status[i]
            # starting with bra
            i_next = i
            j=2
            while j >0
                # If the vector is not yet in a cycle
                if Status[i_next]
                    if j==2
                        # Increment the number of cycles
                        NumOfCycles += 1
                        # Make bra the first member of a cycle
                        MemberID=0
                        Cycles_leaders[NumOfCycles]= i_next
                    end
                    # IdStatus indicates if the "current" vector is in a cycle
                    IdStatus=true
                    # i_next is the serial number of the "following" vector,
                    # While the "following" vector is not yet in a cycle
                    while IdStatus
                        # The "following" vector is now in a cycle
                        Status[i_next]=false
                        MemberID+=1
                        # Add the "following" vector to Cycles

  
                        # i_next is now the serial number of the vector rotated once
                        i_next = serial_num_fast(basis, shiftV(basis.vectors[i_next],M))
                        IdStatus = Status[i_next]
                    end
                else
                    j=0
                end
                j-=1
                if j==1
                    i_next=serial_num_fast(basis,ReverseKet(basis.vectors[i],M))
                    if ~Status[i_next]
                        j=0
                    else
                        Reflection_Extension[NumOfCycles]=true
                    end                
                end
                if j==0
                    CycleSize[NumOfCycles]= MemberID
                end
            end
        end
    end
println("size(Status) = ",Base.summarysize(Status)/1024^3," gb")########################################
    Cycles_leaders, CycleSize, NumOfCycles, Reflection_Extension
end


"""
Create the sparse translational (and reflection) symmetry block of the hamiltonian of bosonic 1D chains with PBC.
"""
function sparse_Block_Diagonal_Hamiltonian_q0R1_gs_OTF(Cycles_leaders:: Vector{Int64}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, t::Float64, U::Float64, L::Int64, N::Int64)
    ### added to get to work like in reg sparse
    mus = zeros(L)
    Nv = zeros(Int64,L)
    Nv_temp = zeros(Int64,L)
   # Nv_rev = zeros(Int64,L)
  #  Nv_old = zeros(Int64,L)
    M=L+N
    nc=N+1
    lc=L+1
    cnkc,jmax=buildpascal(lc,nc)
    #Creating the block H_(q=0,R=1) of the hamiltonian.
    rows = Int64[]
    cols = Int64[]
    elements = Float64[]
    factor=0.0
    for CycleId = 1:NumOfCycles
        # Diagonal part
        Usum = 0.0
        musum = 0.0
        bra =Cycles_leaders[CycleId]
        in2b(bra,cnkc,jmax,lc,nc,Nv)
        for j in 1:L
            Usum += Nv[j] * (Nv[j] - 1)
            musum += mus[j] * Nv[j]
        end
        push!(rows, CycleId)
        push!(cols, CycleId)
        push!(elements, U*Usum/2.0 - musum)
	#CycleMembs=zeros(Int64, L*2)
        # Off-diagonal part
        CycleId_CycleSize=CycleSize[CycleId]
        for j = 1:L
             j_next = j % L + 1
             # Tunnel right, tunnel left.
             for (site1, site2) in [(j, j_next), (j_next, j)]
                 Nv_temp=copy(Nv)
                 if Nv_temp[site1] > 0
                     Nv_temp[site2]+=1
                     Nv_temp[site1]-=1
#                     ket_Leader=Get_ket_leader_OTF(Nv_temp,L,N,jmax,cnkc)

#                    ket =b2in(Nv_temp,cnkc,jmax,lc,nc)
#                    CycleMembs[1]= ket
#                    Nv_rev=copy(Nv_temp)
#                     reverse!(Nv_rev)
#                     ket = b2in(Nv_rev,cnkc,jmax,lc,nc)
#                     CycleMembs[L+1]= ket
#                     for  k= 2:L
#                         circshift!(Nv_old, Nv_temp,1)
#                         Nv_old, Nv_temp = Nv_temp,Nv_old
#                         ket = b2in(Nv_temp,cnkc,jmax,lc,nc)
#                         CycleMembs[k]= ket
#
#                         circshift!(Nv_old, Nv_rev,1)
#                         Nv_old, Nv_rev = Nv_rev,Nv_old
#                         ket = b2in(Nv_rev,cnkc,jmax,lc,nc)
#                         CycleMembs[L+k]= ket
#
#                     end
#                     ket_leaders = minimum(CycleMembs)

                     CycleId1 = serial_num_Cycles_OTF(Cycles_leaders, NumOfCycles, Get_ket_leader_OTF(Nv_temp,L,N,jmax,cnkc))
                     CycleId1_CycleSize = CycleSize[CycleId1]
                     if CycleId >= CycleId1
                     in2b(bra,cnkc,jmax,lc,nc,Nv)

                         factor = - t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)*1.0 * sqrt(Nv[site1]) * sqrt(Nv[site2]+1)
                         push!(rows, CycleId1)
                         push!(cols, CycleId)
                         push!(elements, factor)
                         #push!(rows, CycleId)
                         #push!(cols, CycleId1)
                         #push!(elements, factor)
                     end
                 end
             end
         end
    end
        return Symmetric(sparse(rows, cols, elements, NumOfCycles, NumOfCycles))
#        return  sparse(rows, cols, elements, NumOfCycles, NumOfCycles)
end

function construct_OBDM(Ψ::Vector{Float64},L::Int64,N::Int64,Cycles_leaders::Vector{Int64},CycleSize::Vector{Int64},NumOfCycles::Int64,Reflection_Extension::BitVector)
    U=zeros(Int64,L)
    U_temp=zeros(Int64,L)

    Nv=zeros(Int64,L)
    OBDM=zeros(L)
    M=N+L
    LookupTable,Ones=buildLookupTable()
    for k=1:NumOfCycles
        Ket_Leader=Cycles_leaders[k]
        Ket_temp=Ket_Leader
        getU(Ket_temp,L,LookupTable,Ones,U)
        U_to_Nv(Nv,U,L)
        if Reflection_Extension[k]
            Translational_CycleSize=CycleSize[k]÷2    
        else
            Translational_CycleSize=CycleSize[k]    
        end
        element=0.0
        for i=1:Translational_CycleSize 
            element+=Nv[i]
        end
        OBDM[1]+=element*Ψ[k]*Ψ[k]/Translational_CycleSize
        
        for r=1:L-1
            Ket_temp=Ket_Leader
            element=0.0
            copy!(U_temp,U)
            for i=1:Translational_CycleSize
                if Nv[i]>0
                    Ket_temp2=Ket_temp
                    i_next=(i-1+r)%L+1 
                    bra=destroyParticle(Ket_temp2,U_temp,1)
                    bra=createParticle(bra,U_temp,1+r)
                    CycleId1 = serial_num_Cycles(Cycles_leaders, NumOfCycles,Get_ket_leader(bra,L,M))
                    element+=sqrt(Nv[i]*(Nv[i_next]+1)/(CycleSize[k]*CycleSize[CycleId1]))*Ψ[k]*Ψ[CycleId1]
                end
                Ket_temp=shiftV(Ket_temp,M)
                getU(Ket_temp,L,LookupTable,Ones,U_temp)
            end
            if Reflection_Extension[k]
                Ket_temp=Ket_Leader
                Ket_temp=ReverseKet(Ket_temp,M) 
                getU(Ket_temp,L,LookupTable,Ones,U_temp)
                for i=1:Translational_CycleSize
                    if Nv[L-i+1]>0
                        Ket_temp2=Ket_temp
                        i_next=(i-1+r)%L+1            
                        bra=destroyParticle(Ket_temp2,U_temp,1)
                        bra=createParticle(bra,U_temp,1+r)
                        CycleId1 = serial_num_Cycles(Cycles_leaders,NumOfCycles, Get_ket_leader(bra,L,M))
                        element+=sqrt(Nv[L-i+1]*(Nv[L-i_next+1]+1)/(CycleSize[k]*CycleSize[CycleId1]))*Ψ[k]*Ψ[CycleId1]
                    end
                    Ket_temp=shiftV(Ket_temp,M)
                    getU(Ket_temp,L,LookupTable,Ones,U_temp)
                end
            end
           OBDM[1+r]=element
        end
    end
    return OBDM
end

function construct_OBDM_OTF(Ψ::Vector{Float64},L::Int64,N::Int64,Cycles_leaders::Vector{Int64},CycleSize::Vector{Int64},NumOfCycles::Int64,Reflection_Extension::BitVector)
    nc=N+1
    lc=L+1
    Nv=zeros(Int64,L)
    Nv_temp=zeros(Int64,L)
    Nv_temp2=zeros(Int64,L)


    OBDM=zeros(L)
    cnkc,jmax=buildpascal(lc,nc)
    for k=1:NumOfCycles
        Ket_Leader=Cycles_leaders[k]
        in2b(Ket_Leader,cnkc,jmax,L+1,N+1,Nv)
        if Reflection_Extension[k]
            Translational_CycleSize=CycleSize[k]÷2    
        else
            Translational_CycleSize=CycleSize[k]    
        end
        element=0.0
        for i=1:Translational_CycleSize 
            element+=Nv[i]
        end
        OBDM[1]+=element*Ψ[k]*Ψ[k]/Translational_CycleSize
        for r=1:L-1
            element=0.0
            copy!(Nv_temp,Nv)
            copy!(Nv_temp2,Nv)
            for i=1:Translational_CycleSize
                if Nv[i]>0
                    i_next=(i-1+r)%L+1 
                    Nv_temp[1]-=1
                    Nv_temp[1+r]+=1

                    CycleId1 = serial_num_Cycles_OTF(Cycles_leaders, NumOfCycles,Get_ket_leader_OTF(Nv_temp,L,N,jmax,cnkc))
                    element+=sqrt(Nv[i]*(Nv[i_next]+1)/(CycleSize[k]*CycleSize[CycleId1]))*Ψ[k]*Ψ[CycleId1]
                end
                circshift!(Nv_temp, Nv_temp2,-1)
                copy!(Nv_temp2,Nv_temp)
            end
            if Reflection_Extension[k]
                copy!(Nv_temp,Nv)
                reverse!(Nv_temp)
                copy!(Nv_temp2,Nv_temp)
                for i=1:Translational_CycleSize
                    if Nv[L-i+1]>0
                        i_next=(i-1+r)%L+1  
                        Nv_temp[1]-=1
                        Nv_temp[1+r]+=1
                        CycleId1 = serial_num_Cycles_OTF(Cycles_leaders, NumOfCycles,Get_ket_leader_OTF(Nv_temp,L,N,jmax,cnkc))
                        element+=sqrt(Nv[L-i+1]*(Nv[L-i_next+1]+1)/(CycleSize[k]*CycleSize[CycleId1]))*Ψ[k]*Ψ[CycleId1]
                    end
                circshift!(Nv_temp, Nv_temp2,-1)
                copy!(Nv_temp2,Nv_temp)
                end
            end
           OBDM[1+r]=element
        end
    end
    return OBDM
end

"""
+++++1+++++Create the sparse translational (and reflection) symmetry block of the hamiltonian of bosonic 1D chains with PBC.
"""
function sparse_Block_Diagonal_Hamiltonian_q0R1_gs_reduced1(Cycles_leaders:: Vector{Int64}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, t::Float64, U::Float64, L::Int64, N::Int64)

    ### added to get to work like in reg sparse
    mus = zeros(L)
    U_bra = zeros(Int64,L)
    U_temp = zeros(Int64,L)

    Nv = zeros(Int64,L)
    M=L+N
    LookupTable,NumOnes =buildLookupTable()
    #Creating the block H_(q=0,R=1) of the hamiltonian.
    rows = Int64[]
    cols = Int64[]
    elements = Float64[]
    factor=0.0
    for CycleId = 1:NumOfCycles÷2
        # Diagonal part
        Usum = 0.0
        musum = 0.0
        bra =Cycles_leaders[CycleId]
        bra_temp=copy(bra)
        getU(bra_temp,L, LookupTable,NumOnes, U_bra)
        #U_bra = generateU(bra)
        #Nv = generateN(U_bra)
        U_to_Nv(Nv,U_bra,L)
        for j in 1:L
            Usum += Nv[j] * (Nv[j] - 1)
            musum += mus[j] * Nv[j]
        end
        push!(rows, CycleId)
        push!(cols, CycleId)
        push!(elements, U*Usum/2.0 - musum)
	#CycleMembs=zeros(Int64, L*2)
        # Off-diagonal part
        CycleId_CycleSize=CycleSize[CycleId]
        for j = 1:L
             j_next = j % L + 1
             # Tunnel right, tunnel left.
             for (site1, site2) in [(j, j_next), (j_next, j)]
                 if Nv[site1] > 0
                     ket = bra
                     #U_bra = generateU(bra)
                     U_temp=copy(U_bra)
                     ket = TunnelParticle(ket, U_temp, site1,site2)
                     ket_leader = Get_ket_leader(ket,L,M)

#                     CycleMembs[1]= ket
#                     ket_rev=ReverseKet(ket,M)
#                     CycleMembs[L+1]= ket_rev
#                     for  k= 2:L
#                         ket = shiftV(ket,M)
#                         ket_rev = shiftV(ket_rev,M)
#                         CycleMembs[k]= ket
#                         CycleMembs[L+k]= ket_rev
#                     end
#                     ket_leaders = maximum(CycleMembs)

                     CycleId1 = serial_num_Cycles(Cycles_leaders, NumOfCycles, ket_leader)
                     CycleId1_CycleSize = CycleSize[CycleId1]
                     k1 = CycleId1
                     if CycleId >= k1
                         factor = - t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)*1.0 * sqrt(Nv[site1]) * sqrt(Nv[site2]+1)
                         push!(rows, k1)
        	         push!(cols, CycleId)
                         push!(elements, factor)
                         #push!(rows, CycleId)
                         #push!(cols, k1)
                         #push!(elements, factor)
                     end
                 end
             end
         end
    end
        return Symmetric(sparse(rows, cols, elements, NumOfCycles, NumOfCycles))
#        return  sparse(rows, cols, elements, NumOfCycles, NumOfCycles)

end

"""
+++++2+++++Create the sparse translational (and reflection) symmetry block of the hamiltonian of bosonic 1D chains with PBC.
"""
function sparse_Block_Diagonal_Hamiltonian_q0R1_gs_reduced2(Cycles_leaders:: Vector{Int64}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, t::Float64, U::Float64, L::Int64, N::Int64)

    ### added to get to work like in reg sparse
    mus = zeros(L)
    U_bra = zeros(Int64,L)
    U_temp = zeros(Int64,L)

    Nv = zeros(Int64,L)
    M=L+N
    LookupTable,NumOnes =buildLookupTable()
    #Creating the block H_(q=0,R=1) of the hamiltonian.
    rows = Int64[]
    cols = Int64[]
    elements = Float64[]
    factor=0.0
    for CycleId = NumOfCycles÷2+1:NumOfCycles
        # Diagonal part
        Usum = 0.0
        musum = 0.0
        bra =Cycles_leaders[CycleId]
        bra_temp=copy(bra)
        getU(bra_temp,L, LookupTable,NumOnes, U_bra)
        #U_bra = generateU(bra)
        #Nv = generateN(U_bra)
        U_to_Nv(Nv,U_bra,L)
        for j in 1:L
            Usum += Nv[j] * (Nv[j] - 1)
            musum += mus[j] * Nv[j]
        end
        push!(rows, CycleId)
        push!(cols, CycleId)
        push!(elements, U*Usum/2.0 - musum)
	#CycleMembs=zeros(Int64, L*2)
        # Off-diagonal part
        CycleId_CycleSize=CycleSize[CycleId]
        for j = 1:L
             j_next = j % L + 1
             # Tunnel right, tunnel left.
             for (site1, site2) in [(j, j_next), (j_next, j)]
                 if Nv[site1] > 0
                     ket = bra
                     #U_bra = generateU(bra)
                     U_temp=copy(U_bra)
                     ket = TunnelParticle(ket, U_temp, site1,site2)
                     ket_leader = Get_ket_leader(ket,L,M)

#                     CycleMembs[1]= ket
#                     ket_rev=ReverseKet(ket,M)
#                     CycleMembs[L+1]= ket_rev
#                     for  k= 2:L
#                         ket = shiftV(ket,M)
#                         ket_rev = shiftV(ket_rev,M)
#                         CycleMembs[k]= ket
#                         CycleMembs[L+k]= ket_rev
#                     end
#                     ket_leaders = maximum(CycleMembs)

                     CycleId1 = serial_num_Cycles(Cycles_leaders, NumOfCycles, ket_leader)
                     CycleId1_CycleSize = CycleSize[CycleId1]
                     k1 = CycleId1
                     if CycleId >= k1
                         factor = - t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)*1.0 * sqrt(Nv[site1]) * sqrt(Nv[site2]+1)
                         push!(rows, k1)
        	         push!(cols, CycleId)
                         push!(elements, factor)
                         #push!(rows, CycleId)
                         #push!(cols, k1)
                         #push!(elements, factor)
                     end
                 end
             end
         end
    end
        return Symmetric(sparse(rows, cols, elements, NumOfCycles, NumOfCycles))
#        return  sparse(rows, cols, elements, NumOfCycles, NumOfCycles)

end







# ===========================================reduced=========================================

"""
Create a list of occupation basis for each translational and reflection symmetry cycle for bosonic 1D chains with PBC.
"""
function Translational_Symmetry_Reflection_Cycles_gs_reducedb(L::Int64, N::Int64)

    # Makes a Symbol for each, of type Int64
    for x = [:NumOfCycles,:Num_cycles_max,:MemberID,:i_next]
       @eval $x = Int64
    end
println("======================(Generating the integer basis)======================")########################################
    @time basis = Intbasis(L,N)
println("==================================(end.)==================================")########################################
    IdStatus = Bool
    D = basis.D  #length(basis)
    #L= basis.L
    #N= basis.N
    M=L+N
    # Approximate the maximum number of cycles in the basis
    Num_cycles_max = round(Int64,((D/L/2)+6*(D/L/2)^.5))
    # Make an empty array with a row for each cycle and L (# of sites) rows
    Cycles_leaders = zeros(Int64, Num_cycles_max)
    # Make an empty array with number_of_cycles values
    CycleSize = zeros(Int64, Num_cycles_max)
    # Give every vector in the basis a state of true, it's not in a cycle yet
    Status  = trues(D)
    Reflection_Extension  = falses(Num_cycles_max)
    NumOfCycles = 0
    MemberID=0
    # For each vector in the basis
        for i=1: D     
        if Status[i]
            # starting with bra
            i_next = i
            j=2
            while j >0
                # If the vector is not yet in a cycle
                if Status[i_next]
                    if j==2
                        # Increment the number of cycles
                        NumOfCycles += 1
                        # Make bra the first member of a cycle
                        MemberID=0
                        Cycles_leaders[NumOfCycles]= basis.vectors[i_next]
                    end
                    # IdStatus indicates if the "current" vector is in a cycle
                    IdStatus=true
                    # i_next is the serial number of the "following" vector,
                    # While the "following" vector is not yet in a cycle
                    while IdStatus
                        # The "following" vector is now in a cycle
                        Status[i_next]=false
                        MemberID+=1  
                        # i_next is now the serial number of the vector rotated once
                        i_next = serial_num_fast(basis, shiftV(basis.vectors[i_next],M))
                        IdStatus = Status[i_next]
                    end
                else
                    j=0
                end
                j-=1
                if j==1
                    i_next=serial_num_fast(basis,ReverseKet(basis.vectors[i],M))
                    if ~Status[i_next]
                        j=0
                    else
                        Reflection_Extension[NumOfCycles]=true
                    end                
                end
                if j==0
                    CycleSize[NumOfCycles]= MemberID
                end
            end
        end
    end
#    println("RNumOfCycles=",NumOfCycles)
#    totc=0
#    for tt=1:NumOfCycles
#        if Reflection_Extension[tt]
#            totc+=1
#        end
#    end
#    println("NumOfCycles=",totc*2+NumOfCycles-totc)
    Cycles_leaders, CycleSize, NumOfCycles, Reflection_Extension
end


"""
Create a list of occupation basis for each translational and reflection symmetry cycle for bosonic 1D chains with PBC.
"""
function Translational_Symmetry_Reflection_Cycles_gs_OTFb(L::Int64, N::Int64)

    # Makes a Symbol for each, of type Int64
    for x = [:NumOfCycles,:Num_cycles_max,:MemberID,:i_next]
       @eval $x = Int64
    end
    basis = Intbasis(L,N)
    IdStatus = Bool
    D = basis.D  #length(basis)
    #L= basis.L
    #N= basis.N
    M=L+N
    # Approximate the maximum number of cycles in the basis
    Num_cycles_max = round(Int64,((D/L/2)+6*(D/L/2)^.5))
    # Make an empty array with a row for each cycle and L (# of sites) rows
    Cycles_leaders = zeros(Int64, Num_cycles_max)
    # Make an empty array with number_of_cycles values
    CycleSize = zeros(Int64, Num_cycles_max)
    # Give every vector in the basis a state of true, it's not in a cycle yet
    Status  = trues(D)
    Reflection_Extension  = falses(Num_cycles_max)
    NumOfCycles = 0
    MemberID=0
    # For each vector in the basis
        for i=1: D     
        if Status[i]
            # starting with bra
            i_next = i
            j=2
            while j >0
                # If the vector is not yet in a cycle
                if Status[i_next]
                    if j==2
                        # Increment the number of cycles
                        NumOfCycles += 1
                        # Make bra the first member of a cycle
                        MemberID=0
                        Cycles_leaders[NumOfCycles]= i_next
                    end
                    # IdStatus indicates if the "current" vector is in a cycle
                    IdStatus=true
                    # i_next is the serial number of the "following" vector,
                    # While the "following" vector is not yet in a cycle
                    while IdStatus
                        # The "following" vector is now in a cycle
                        Status[i_next]=false
                        MemberID+=1
                        # Add the "following" vector to Cycles

  
                        # i_next is now the serial number of the vector rotated once
                        i_next = serial_num_fast(basis, shiftV(basis.vectors[i_next],M))
                        IdStatus = Status[i_next]
                    end
                else
                    j=0
                end
                j-=1
                if j==1
                    i_next=serial_num_fast(basis,ReverseKet(basis.vectors[i],M))
                    if ~Status[i_next]
                        j=0
                    else
                        Reflection_Extension[NumOfCycles]=true
                    end                
                end
                if j==0
                    CycleSize[NumOfCycles]= MemberID
                end
            end
        end
    end
#println("size(Status) = ",Base.summarysize(Status)/1024^3," gb")########################################
    Cycles_leaders, CycleSize, NumOfCycles, Reflection_Extension
end




"""
Create the sparse translational (and reflection) symmetry block of the hamiltonian of bosonic 1D chains with PBC.
"""
function sparse_Block_Diagonal_Hamiltonian_q0R1_gs_reduced_comp(t::Float64, U::Float64, L::Int64, N::Int64)

    ### added to get to work like in reg sparse
    mus = zeros(L)
    U_bra = zeros(Int64,L)
    U_temp = zeros(Int64,L)
println("======================(Generating the Cycles)======================")########################################

    @time Cycles_leaders, CycleSize, NumOfCycles,Reflection_Extension=Translational_Symmetry_Reflection_Cycles_gs_reducedb(L,N)

println("==================================(end.)==================================")########################################

    Reflection_Extension =Nothing #################


    Nv = zeros(Int64,L)
    M=L+N
    LookupTable,NumOnes =buildLookupTable()
    #Creating the block H_(q=0,R=1) of the hamiltonian.
    rows = Int64[]
    cols = Int64[]
    elements = Float64[]
    factor=0.0
    for CycleId = 1:NumOfCycles
        # Diagonal part
        Usum = 0.0
        musum = 0.0
        bra =Cycles_leaders[CycleId]
        bra_temp=copy(bra)
        getU(bra_temp,L, LookupTable,NumOnes, U_bra)
        #U_bra = generateU(bra)
        #Nv = generateN(U_bra)
        U_to_Nv(Nv,U_bra,L)
        for j in 1:L
            Usum += Nv[j] * (Nv[j] - 1)
            musum += mus[j] * Nv[j]
        end
        push!(rows, CycleId)
        push!(cols, CycleId)
        push!(elements, U*Usum/2.0 - musum)
	#CycleMembs=zeros(Int64, L*2)
        # Off-diagonal part
        CycleId_CycleSize=CycleSize[CycleId]
        for j = 1:L
             j_next = j % L + 1
             # Tunnel right, tunnel left.
             for (site1, site2) in [(j, j_next), (j_next, j)]
                 if Nv[site1] > 0
                     ket = bra
                     #U_bra = generateU(bra)
                     U_temp=copy(U_bra)
                     ket = TunnelParticle(ket, U_temp, site1,site2)
                     ket_leader = Get_ket_leader(ket,L,M)

#                     CycleMembs[1]= ket
#                     ket_rev=ReverseKet(ket,M)
#                     CycleMembs[L+1]= ket_rev
#                     for  k= 2:L
#                         ket = shiftV(ket,M)
#                         ket_rev = shiftV(ket_rev,M)
#                         CycleMembs[k]= ket
#                         CycleMembs[L+k]= ket_rev
#                     end
#                     ket_leaders = maximum(CycleMembs)

                     CycleId1 = serial_num_Cycles(Cycles_leaders, NumOfCycles, ket_leader)
                     CycleId1_CycleSize = CycleSize[CycleId1]
                     k1 = CycleId1
                     if CycleId >= k1
                         factor = - t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)*1.0 * sqrt(Nv[site1]) * sqrt(Nv[site2]+1)
                         push!(rows, k1)
        	         push!(cols, CycleId)
                         push!(elements, factor)
                         #push!(rows, CycleId)
                         #push!(cols, k1)
                         #push!(elements, factor)
                     end
                 end
             end
         end
    end

    #Cycles_leaders =Nothing
    #CycleSize =Nothing
    #GC.gc()


        return Symmetric(sparse(rows, cols, elements, NumOfCycles, NumOfCycles)), NumOfCycles
#        return  sparse(rows, cols, elements, NumOfCycles, NumOfCycles)

end
