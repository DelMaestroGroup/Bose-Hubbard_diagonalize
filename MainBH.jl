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
    "-n"
        help = "particle partition size (default: N/2)"
        arg_type = Int64

    "--output", "-o"
        metavar = "FILE"

    "--getspectrum"
        help = "report the spectrum of the RDM."
        action = :store_true

    "--part"
        help = "calculate particle entanglement."
        action = :store_true

    end

    return parse_args(s, as_symbols=true)
end

function main()
    parsed_args = parse_commandline()
    L,N,T = parsed_args[:L], parsed_args[:N], parsed_args[:T]
    isnothing(parsed_args[:A]) ? A = L ÷ 2 : A = parsed_args[:A]
    isnothing(parsed_args[:n]) ? n = N ÷ 2 : n = parsed_args[:n]
    output = parsed_args[:output]
    getspectrum= parsed_args[:getspectrum]
    part = parsed_args[:part]

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

    if part
        println(" L=$(L), N=$(N), ℓ=$(A), n=$(n), T=$(T), U_min=$(U0), U_min=$(UM), u_num=$(NumU), sym gs")
    else
        println(" L=$(L), N=$(N), ℓ=$(A), T=$(T), U_min=$(U0), U_min=$(UM), u_num=$(NumU), sym gs")
    end

    if isnothing(output) 
        if part
            output = (@sprintf "spatEE_gs_%02d_%02d_%02d_%02d_%+5.3f_%+5.3f_%003d.dat" L N A n U0 UM NumU)
        else
            output = (@sprintf "spatEE_gs_%02d_%02d_%02d_%+5.3f_%+5.3f_%003d.dat" L N A U0 UM NumU)
        end
    end
    outputPn = (@sprintf "spatEE_gs_Pn_%02d_%02d_%02d_%+5.3f_%+5.3f_%003d.dat" L N A U0 UM NumU)
    outputS1 = (@sprintf "spatEE_gs_S1_%02d_%02d_%02d_%+5.3f_%+5.3f_%003d.dat" L N A U0 UM NumU)
    outputS2 = (@sprintf "spatEE_gs_S2_%02d_%02d_%02d_%+5.3f_%+5.3f_%003d.dat" L N A U0 UM NumU)
    outputSNeg = (@sprintf "spatEE_gs_SNeg_%02d_%02d_%02d_%+5.3f_%+5.3f_%003d.dat" L N A U0 UM NumU)

    f = open(output, "w")
    fPn = open(outputPn, "w")
    fS1 = open(outputS1, "w")
    fS2 = open(outputS2, "w")
    fSNeg = open(outputSNeg, "w")

    write(f, "# L=$(L), N=$(N), T=$(T), sym gs\n")
    if part
        write(f,@sprintf "#%11s%11s%24s%24s%24s%24s%24s%24s%24s%24s\n" "U" "Eg" "S₁(ℓ=$(A))" "S₁acc(ℓ=$(A))" "S₂(ℓ=$(A))" "S₂acc(ℓ=$(A))" "EN(ℓ=$(A))" "ENacc(ℓ=$(A))" "S₁(n=$(n))" "S₂(n=$(n))")
    else
        write(f,@sprintf "#%11s%11s%24s%24s%24s%24s%24s%24s\n" "U" "Eg" "S₁(ℓ=$(A))" "S₁acc(ℓ=$(A))" "S₂(ℓ=$(A))" "S₂acc(ℓ=$(A))" "EN(ℓ=$(A))" "ENacc(ℓ=$(A))")
    end


    write(fPn, "# L=$(L), N=$(N), T=$(T), ℓ=$(A), sym gs\n")
    write(fPn, "# Pn\n")
    write(fPn,@sprintf "#%11s" "U")
    for i=0:N
        write(fPn,@sprintf "%24s" "n=$(i)")
    end    
    write(fPn, "\n")

    write(fS1, "# L=$(L), N=$(N), T=$(T), ℓ=$(A), sym gs\n")
    write(fS1, "# S1\n")
    write(fS1,@sprintf "#%11s" "U")
    for i=0:N
        write(fS1,@sprintf "%24s" "n=$(i)")
    end    
    write(fS1, "\n")

    write(fS2, "# L=$(L), N=$(N), T=$(T), ℓ=$(A), sym gs\n")
    write(fS2, "# S2\n")
    write(fS2,@sprintf "#%11s" "U")
    for i=0:N
        write(fS2,@sprintf "%24s" "n=$(i)")
    end    
    write(fS2, "\n")

    write(fSNeg, "# L=$(L), N=$(N), T=$(T), ℓ=$(A), sym gs\n")
    write(fSNeg, "# SNeg\n")
    write(fSNeg,@sprintf "#%11s" "U")
    for i=0:N
        write(fSNeg,@sprintf "%24s" "n=$(i)")
    end    
    write(fSNeg, "\n")






    basis = Intbasis(L,N)
    println(" Bosonic integer basis generated")
    Cycles_leaders, CycleSize, NumOfCycles, InvCycles_Id = Translational_Symmetry_Reflection_Cycles_gs(basis)
    println("  Symmetries explored")
    Ψ = ones(Float64, NumOfCycles)
    Ψ.= Ψ./sqrt(dot(Ψ,Ψ))
    Eg=0.0
    for U in U_range
        H = sparse_Block_Diagonal_Hamiltonian_q0R1_gs(basis, Cycles_leaders, CycleSize, NumOfCycles, InvCycles_Id, T, U)
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
        S1_sp,S1_acc,S2_sp,S2_acc,Neg,Neg_acc,Pn,S1,S2,SNeg = spatial_entropy_op_Ts_gs(basis, A, d, CycleSize,InvCycles_Id, getspectrum, U)


        if part
println("======================(Calculating The entanglement)======================")########################################
            @time S1_part,S2_part=particle_entropy(basis, n, d, CycleSize, InvCycles_Id)
println("==================================(end.)==================================")########################################
        end

        println("     Entanglement calculated, for U=$(U)")
        if part
            write(f, @sprintf "%12.6f%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E\n" U Eg S1_sp S1_acc S2_sp S2_acc Neg Neg_acc S1_part S2_part)
            flush(f)
        else
            write(f, @sprintf "%12.6f%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E\n" U Eg S1_sp S1_acc S2_sp S2_acc Neg Neg_acc)
            flush(f)
        end

        write(fPn,@sprintf "%12.6f" U)
        for i=1:N+1
            write(fPn,@sprintf "%24.12E" Pn[i])
        end    
        write(fPn, "\n")
        flush(fPn)

        write(fS1,@sprintf "%12.6f" U)
        for i=1:N+1
            write(fS1,@sprintf "%24.12E" S1[i])
        end    
        write(fS1, "\n")
        flush(fS1)

        write(fS2,@sprintf "%12.6f" U)
        for i=1:N+1
            write(fS2,@sprintf "%24.12E" S2[i])
        end    
        write(fS2, "\n")
        flush(fS2)

        write(fSNeg,@sprintf "%12.6f" U)
        for i=1:N+1
            write(fSNeg,@sprintf "%24.12E" SNeg[i])
        end    
        write(fSNeg, "\n")
        flush(fSNeg)


    end
    close(f)
    close(fPn)
    close(fS1)
    close(fS2)
    close(fSNeg)

end

main()
