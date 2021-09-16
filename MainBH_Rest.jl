push!(LOAD_PATH, joinpath(dirname(@__FILE__), "src"))
push!(LOAD_PATH, joinpath(dirname(@__FILE__), "IntegerBosonBasis/src"))
using IntegerBasisDiagonalize
using IntegerBosonBasisG
using Serialization
using ArgParse
using Arpack
using Printf
using LinearAlgebra
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
    "L"
        help = "number of sites"
        arg_type = Int
    "N"
        help = "number of particles"
        arg_type = Int
    "n_max"
        help = "Site capacity"
        arg_type = Int
    "--T"
        metavar = "T"
        arg_type = Float64
        default = 1.0
    "--u_min"
        metavar = "U"
        help = "minimum U"
        arg_type = Float64
        default = 0.0
    "--u_max"
        metavar = "U"
        help = "maximum U"
        arg_type = Float64
        default = 1.0
    "--u_step"
        metavar = "U"
        help = "U step"
        arg_type = Float64
    "--u_num"
        metavar = "Num"
        help = "number of U"
        arg_type = Int
    "--u_log"
        help = "use logarithmic scale for U"
        action = :store_true
    "-A"
        help = "size of region A (default: L/2)"
        arg_type = Int64

    "--output", "-o"
        metavar = "FILE"

    end

    return parse_args(s, as_symbols=true)
end

function main()
    parsed_args = parse_commandline()
    L,N,T = parsed_args[:L], parsed_args[:N], parsed_args[:T]
    n_max = parsed_args[:n_max]
    isnothing(parsed_args[:A]) ? A = L ÷ 2 : A = parsed_args[:A]
    output = parsed_args[:output]

    if isnothing(L) && isnothing(N)
        throw("System size must be given.")
    end


    if parsed_args[:u_log] && isnothing(parsed_args[:u_num]) 
        throw("--u-log must be used with --u-num")
    end

    if isnothing(parsed_args[:u_step]) 
        if isnothing(parsed_args[:u_num])
            U_range = parsed_args[:u_min]:0.5: parsed_args[:u_max]
        else
            if parsed_args[:u_log]
                U_range = (10.0).^LinRange(parsed_args[:u_min], parsed_args[:u_max], parsed_args[:u_num])
            else
                U_range = LinRange(parsed_args[:u_min], parsed_args[:u_max], parsed_args[:u_num])
            end
        end
    else
        if isnothing(parsed_args[:u_num]) 
            U_range = parsed_args[:u_min]: parsed_args[:u_step]: parsed_args[:u_max]
        else
            throw("--u-step and --u-num may not both be supplied") 
        end
    end

    U0=U_range[1]
    UM=U_range[length(U_range)]
    NumU=length(U_range)
    println(" L=$(L), N=$(N), ℓ=$(A), T=$(T), U_min=$(U0), U_min=$(UM), u_num=$(NumU), sym gs, n_max=$(n_max)")
 
    if isnothing(output) 
        output = (@sprintf "spatEE_gs_%02d_%02d_%02d_%02d_%+5.3f_%+5.3f_%003d.dat" n_max L N A U0 UM NumU)
    end

    f = open(output, "w")

    write(f, "# L=$(L), N=$(N), T=$(T), sym gs, n_max=$(n_max)\n")
    write(f,@sprintf "#%11s%11s%24s%24s%24s%24s%24s%24s\n" "U" "Eg" "S₁(ℓ=$(A))" "S₁acc(ℓ=$(A))" "S₂(ℓ=$(A))" "S₂acc(ℓ=$(A))" "EN(ℓ=$(A))" "ENacc(ℓ=$(A))")






    basis = RestrictedIntbasis(L,N,n_max)
    println(" Bosonic integer basis generated")
    Cycles_leaders, CycleSize, NumOfCycles, InvCycles_Id = Translational_Symmetry_Reflection_Cycles_gs(basis)
    println("  Symmetries explored")
    Ψ = ones(Float64, NumOfCycles)
    Ψ.= Ψ./sqrt(dot(Ψ,Ψ))
    Eg=0.0
    for U in U_range
        H = sparse_Block_Diagonal_Hamiltonian_q0R1_gs_Rest(basis, Cycles_leaders, CycleSize, NumOfCycles, InvCycles_Id, T, U)
        println("   Hamiltonian generated for U=$(U)")
        E,d = eigs(H,nev=1,which=:SR,tol=1e-13,v0=Ψ)
        println("    Ground state obtained, for U=$(U)")
        H=Nothing
        GC.gc()
        Eg=Real(E[1])
        d=vec(d)
        for k=1:length(d)
            Ψ[k]=Real(d[k])
        end
        S1_sp,S1_acc,S2_sp,S2_acc,Neg,Neg_acc = spatial_entropy_acc_gs_Rest(basis,A,d, CycleSize,InvCycles_Id)
        println("     Entanglement calculated, for U=$(U)")
        write(f, @sprintf "%12.6f%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E\n" U Eg S1_sp S1_acc S2_sp S2_acc Neg Neg_acc)
        flush(f)

    end

end

main()
