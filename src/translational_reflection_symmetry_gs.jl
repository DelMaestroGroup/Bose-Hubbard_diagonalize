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
