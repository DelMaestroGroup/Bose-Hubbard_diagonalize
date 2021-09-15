"""
Create a list of occupation basis for each translational and reflection symmetry cycle for bosonic 1D chains with PBC.
"""
function Symmetry_Reflection_Cycles_gs(basis::AbstractIntbasis)

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
    Num_cycles_max = round(Int64,((D/2)+6*(D/2)^.5))
    # Make an empty array with a row for each cycle and L (# of sites) rows
    Cycles_leaders = zeros(Int64, Num_cycles_max)
    # Make an empty array with number_of_cycles values
    CycleSize = zeros(Int64, Num_cycles_max)
    InvCycles_Id = zeros(Int64, D)
    # Give every vector in the basis a state of true, it's not in a cycle yet
    Status  = trues(D)
    NumOfCycles = 0
    NumOfCycles_neg = 0
    MemberID=0
    # For each vector in the basis
    for i=1: D    
        # starting with bra
        if Status[i]
            # If the vector is not yet in a cycle
            NumOfCycles += 1
            # Increment the number of cycles
            MemberID=1
            Cycles_leaders[NumOfCycles]= i
            InvCycles_Id[i] = NumOfCycles
            Status[i]=false
            i_next=serial_num_fast(basis,ReverseKet(basis.vectors[i],M))
            if Status[i_next]
                NumOfCycles_neg +=1
                MemberID +=1
                InvCycles_Id[i_next] = NumOfCycles
                Status[i_next]=false
            end  
            CycleSize[NumOfCycles]= MemberID
        end
    end
    Cycles_leaders, CycleSize, NumOfCycles, NumOfCycles_neg, InvCycles_Id
end
