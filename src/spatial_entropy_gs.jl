# From https://github.com/MelkoCollective/BH_diagonalize/tree/master/src

"""
	spatial_entropy_op_Ts_gs(basis::AbstractIntbasis, Asize::Int, d::Vector{T}, CycleSize::Array{Int64,1}, InvCycles_Id::Vector{Int64}, getspectrum::Bool, U::Float64) where {T<:Number}

Calculate both the spatial and the operational entanglement entropies of a
region A, using the SVD. The latter is the "entanglement of particles"
introduced by Wiseman and Vaccaro in 2003.

Asize is an integer representing the number of sites in the first region, so `Asize < basis.L`, converted into an array of integers from 1 to Asize. d is a vector of the wavefunction coefficients, i.e. for sparse hamiltonian Ψ, d is vec(eigs(Ψ, nev=1)[2]), the coefficients for the wavefunctions.
"""


using Printf
"""
Spatial entropy using symmetries with S2_sp and Neg.
"""
function spatial_entropy_op_Ts_gs(basis::AbstractIntbasis, Asize::Int, d::Vector{T}, CycleSize::Array{Int64,1}, InvCycles_Id::Vector{Int64}, getspectrum::Bool, U::Float64) where {T<:Number}

	A,B = 1:Asize, Asize+1:basis.L

	# Normalize accounting for symmetries
	for i = 1:length(d)
	    d[i] /= sqrt(CycleSize[i])
	end
        L=basis.L
        N=basis.N
        M=L+N
	# Matrices to SVD
        norms = zeros(Float64, N+1)
        S1 = zeros(Float64, N+1)
        S2 = zeros(Float64, N+1)
        SNeg= zeros(Float64, N+1)
	if getspectrum
		output = (@sprintf "spatEE_gs_spectrum_%02d_%02d_%02d_%+5.3f.dat" L N Asize U)
		fs = open(output, "w")
		write(fs, "# L=$(L), N=$(N), U=$(U),ℓ=$(Asize), sym gs\n")
                write(fs,@sprintf "#%4s%24s\n" "n" "eigenvalues")
	end
	for i=0:basis.N
		#Amatrices = []
		basisA = Intbasis(Asize,i)
		basisB = Intbasis(L-Asize,N-i) 
		Cycles_leadersA, CycleSizeA, NumOfCyclesA, NumOfCycles_negA, InvCycles_IdA =Symmetry_Reflection_Cycles_gs(basisA) 
		Cycles_leadersB, CycleSizeB, NumOfCyclesB, NumOfCycles_negB, InvCycles_IdB =Symmetry_Reflection_Cycles_gs(basisB)           
		DimA = basisA.D
		DimB = basisB.D
		Amatrix1 =zeros(Float64, NumOfCyclesA, NumOfCyclesB)
		Amatrix2 =zeros(Float64, NumOfCycles_negA, NumOfCycles_negB)
		MA=basisA.L+basisA.N
                iA2=0
		for iA1 =1: NumOfCyclesA
                	braA= basisA.vectors[Cycles_leadersA[iA1]]
			if CycleSizeA[iA1]>1
				iA2+=1
			end
			iB2=0
			for iB1 =1: NumOfCyclesB
               			braB= basisB.vectors[Cycles_leadersB[iB1]] 
				if CycleSizeB[iB1]>1
					iB2+=1
				end
				bra1=(braB<<(basisA.L+ basisA.N))|braA
				bra2=(braB<<(basisA.L+ basisA.N))| ReverseKet(braA,MA)
				Amatrix1[iA1,iB1]=(d[InvCycles_Id[serial_num_fast(basis, bra1)]]+d[InvCycles_Id[serial_num_fast(basis, bra2)]])/2.0*(CycleSizeA[iA1]* CycleSizeB[iB1])^0.5
				if (CycleSizeA[iA1]+CycleSizeB[iB1])>3
					Amatrix2[iA2,iB2]=(d[InvCycles_Id[serial_num_fast(basis, bra1)]]-d[InvCycles_Id[serial_num_fast(basis, bra2)]])
				end
			end
		end
#println("size(Amatrix1(n=$(i))) = ",Base.summarysize(Amatrix1)/1024^3," gb")########################################
#println("size(Amatrix2(n=$(i))) = ",Base.summarysize(Amatrix2)/1024^3," gb")########################################
		Sn1=svdvals!(Amatrix1)
		Sn2=svdvals!(Amatrix2)
                Sn12=[s for s in [Sn1,Sn2]]
                Sn=vcat(Sn12...)
                Amatrix1 =Nothing
		Amatrix2 =Nothing
                GC.gc()

		norms[i+1]=sum(Sn.^2)
                S2[i+1]=sum(Sn.^4)
                SNeg[i+1]=sum(Sn) 
		for k=1:length(Sn)
			if abs(Sn[k])>0
				S1[i+1] -= abs(Sn[k])^2*log(abs(Sn[k])^2)
			end
		end
		Sn=Sn.^2
		if getspectrum

			sort!(Sn,rev=true)
			write(fs,@sprintf "%24s" "$(i)")
			for lambda in Sn
				write(fs, @sprintf "%24.12E" lambda)
			end
			write(fs, "\n")
			flush(fs)
		end
	end

        if getspectrum
		close(fs)
	end

	norm_err = abs(sum(norms) - 1.0)
	if norm_err > 1e-12
		@warn("norm error: $(norm_err)")
	end

	# Spatial.
        S1_sp=sum(S1)
	S2_sp = -log(sum(S2))
        Neg=2.0*log(sum(SNeg))
	# Operational.
	S1_op = S1_sp
	for n in norms
		if n>0
			S1_op += log(n) * n
		end
	end
        S2_op=-2.0*log(sum(S2.^(0.5)))
        Neg_op=log(sum(SNeg.^(2.0)))


	S1_sp,S1_op,S2_sp,S2_op,Neg,Neg_op,norms,S1,S2,SNeg
end


#######################################################################
#######################################################################
#######################################################################
#######################################################################

"""
Spatial entropy using symmetries with S2_sp and Neg(scaling).
"""
function spatial_entropy_op_Ts_gs_scaling(basis::AbstractIntbasis, Asize::Int, d::Vector{T}, CycleSize::Array{Int64,1}, InvCycles_Id::Vector{Int64}) where {T<:Number}
fileout="out.dat"
	A,B = 1:Asize, Asize+1:basis.L

	# Normalize accounting for symmetries
	for i = 1:length(d)
	    d[i] /= sqrt(CycleSize[i])
	end
        L=basis.L
        N=basis.N
        M=L+N
	# Matrices to SVD
        norms = zeros(Float64, N+1)
        S1 = zeros(Float64, N+1)
        S2 = zeros(Float64, N+1)
        SNeg= zeros(Float64, N+1)
	for i=0:basis.N
		#Amatrices = []
		basisA = Intbasis(Asize,i)
		basisB = Intbasis(L-Asize,N-i) 
		Cycles_leadersA, CycleSizeA, NumOfCyclesA, NumOfCycles_negA, InvCycles_IdA =Symmetry_Reflection_Cycles_gs(basisA) 
		Cycles_leadersB, CycleSizeB, NumOfCyclesB, NumOfCycles_negB, InvCycles_IdB =Symmetry_Reflection_Cycles_gs(basisB)   
		DimA = basisA.D
		DimB = basisB.D
		Amatrix1 =zeros(Float64, NumOfCyclesA, NumOfCyclesB)
		Amatrix2 =zeros(Float64, NumOfCycles_negA, NumOfCycles_negB)
		MA=basisA.L+basisA.N
                iA2=0
		for iA1 =1: NumOfCyclesA
                	braA= basisA.vectors[Cycles_leadersA[iA1]]
			if CycleSizeA[iA1]>1
				iA2+=1
			end
			iB2=0
			for iB1 =1: NumOfCyclesB
               			braB= basisB.vectors[Cycles_leadersB[iB1]]
				if CycleSizeB[iB1]>1
					iB2+=1
				end
				bra1=((braB<<(basisA.L+ basisA.N))|braA)
				bra2=((braB<<(basisA.L+ basisA.N))|(ReverseKet(braA,MA)))
				Amatrix1[iA1,iB1]=(d[InvCycles_Id[serial_num_fast(basis, bra1)]]+d[InvCycles_Id[serial_num_fast(basis, bra2)]])/2.0*(CycleSizeA[iA1]* CycleSizeB[iB1])^0.5
				if (CycleSizeA[iA1]+CycleSizeB[iB1])>3
					Amatrix2[iA2,iB2]=(d[InvCycles_Id[serial_num_fast(basis, bra1)]]-d[InvCycles_Id[serial_num_fast(basis, bra2)]])
				end
			end
		end
println("size(Amatrix1(n=$(i))) = ",Base.summarysize(Amatrix1)/1024^3," gb")########################################
println("size(Amatrix2(n=$(i))) = ",Base.summarysize(Amatrix2)/1024^3," gb")########################################
		@time Sn1=svdvals!(Amatrix1)
		@time Sn2=svdvals!(Amatrix2)
                Sn12=[s for s in [Sn1,Sn2]]
                Sn=vcat(Sn12...)
                Amatrix1 =Nothing
		Amatrix2 =Nothing
                #GC.gc()
		norms[i+1]=sum(Sn.^2)
                S2[i+1]=sum(Sn.^4)
                SNeg[i+1]=sum(Sn) 
		for k=1:length(Sn)
			if abs(Sn[k])>0
				S1[i+1] -= abs(Sn[k])^2*log(abs(Sn[k])^2)
			end
		end
	end
	norm_err = abs(sum(norms) - 1.0)
	if norm_err > 1e-12
		@warn("norm error: $(norm_err)")
	end

	# Spatial.
        S1_sp=sum(S1)
	S2_sp = -log(sum(S2))
        Neg=2.0*log(sum(SNeg))
	# Operational.
	S1_op = S1_sp
	for n in norms
		if n>0
			S1_op += log(n) * n
		end
	end
        S2_op=-2.0*log(sum(S2.^(0.5)))
        Neg_op=log(sum(SNeg.^(2.0)))


	S1_sp,S1_op,S2_sp,S2_op,Neg,Neg_op,norms,S1,S2,SNeg
end

############################
############################
############################
"""
Using RestrictedIntbasis
Spatial entropy using symmetries with S2_sp and Neg.
"""
function spatial_entropy_op_Ts_gs_Rest(basis::AbstractIntbasis, Asize::Int, d::Vector{T}, CycleSize::Array{Int64,1}, InvCycles_Id::Vector{Int64}) where {T<:Number}

	A,B = 1:Asize, Asize+1:basis.L

	# Normalize accounting for symmetries
	for i = 1:length(d)
	    d[i] /= sqrt(CycleSize[i])
	end
        n_max=basis.n_max
        L=basis.L
        N=basis.N
        M=L+N
	# Matrices to SVD
        norms = zeros(Float64, N+1)
        S1 = zeros(Float64, N+1)
        S2 = zeros(Float64, N+1)
        SNeg= zeros(Float64, N+1)
	for i=0:N
		if  (i<=n_max* Asize)&&((N-i)<=n_max*(L-Asize))
			#Amatrices = []
			basisA = RestrictedIntbasis(Asize,i, n_max)
			basisB = RestrictedIntbasis(L-Asize,N-i, n_max) 
			Cycles_leadersA, CycleSizeA, NumOfCyclesA, NumOfCycles_negA, InvCycles_IdA =Symmetry_Reflection_Cycles_gs(basisA) 
			Cycles_leadersB, CycleSizeB, NumOfCyclesB, NumOfCycles_negB, InvCycles_IdB =Symmetry_Reflection_Cycles_gs(basisB)           
			DimA = basisA.D
			DimB = basisB.D
			Amatrix1 =zeros(Float64, NumOfCyclesA, NumOfCyclesB)
			Amatrix2 =zeros(Float64, NumOfCycles_negA, NumOfCycles_negB)
			MA=basisA.L+basisA.N
			iA2=0
			for iA1 =1: NumOfCyclesA
        	        	braA= basisA.vectors[Cycles_leadersA[iA1]]
				if CycleSizeA[iA1]>1
					iA2+=1
				end
				iB2=0
				for iB1 =1: NumOfCyclesB
               				braB= basisB.vectors[Cycles_leadersB[iB1]] 
					if CycleSizeB[iB1]>1
						iB2+=1
					end
					bra1=(braB<<(basisA.L+ basisA.N))|braA
					bra2=(braB<<(basisA.L+ basisA.N))| ReverseKet(braA,MA)
					Amatrix1[iA1,iB1]=(d[InvCycles_Id[serial_num_fast(basis, bra1)]]+d[InvCycles_Id[serial_num_fast(basis, bra2)]])/2.0*(CycleSizeA[iA1]* CycleSizeB[iB1])^0.5
					if (CycleSizeA[iA1]+CycleSizeB[iB1])>3
						Amatrix2[iA2,iB2]=(d[InvCycles_Id[serial_num_fast(basis, bra1)]]-d[InvCycles_Id[serial_num_fast(basis, bra2)]])
					end
				end
			end
#println("size(Amatrix1(n=$(i))) = ",Base.summarysize(Amatrix1)/1024^3," gb")########################################
#println("size(Amatrix2(n=$(i))) = ",Base.summarysize(Amatrix2)/1024^3," gb")########################################
			Sn1=svdvals!(Amatrix1)
			Sn2=svdvals!(Amatrix2)
                	Sn12=[s for s in [Sn1,Sn2]]
 	        	Sn=vcat(Sn12...)
        	        Amatrix1 =Nothing
			Amatrix2 =Nothing
                	GC.gc()

			norms[i+1]=sum(Sn.^2)
                	S2[i+1]=sum(Sn.^4)
                	SNeg[i+1]=sum(Sn) 
			for k=1:length(Sn)
				if abs(Sn[k])>0
					S1[i+1] -= abs(Sn[k])^2*log(abs(Sn[k])^2)
				end
			end
		end
	end
	norm_err = abs(sum(norms) - 1.0)
	if norm_err > 1e-12
		@warn("norm error: $(norm_err)")
	end

	# Spatial.
        S1_sp=sum(S1)
	S2_sp = -log(sum(S2))
        Neg=2.0*log(sum(SNeg))
	# Operational.
	S1_op = S1_sp
	for n in norms
		if n>0
			S1_op += log(n) * n
		end
	end
        S2_op=-2.0*log(sum(S2.^(0.5)))
        Neg_op=log(sum(SNeg.^(2.0)))


	S1_sp,S1_op,S2_sp,S2_op,Neg,Neg_op,norms,S1,S2,SNeg
end


