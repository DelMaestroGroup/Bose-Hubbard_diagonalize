module IntegerBasisDiagonalize

using SparseArrays
using LinearAlgebra: svdvals!,Symmetric, svdvals
using Arpack: eigs
using IntegerBosonBasisG

export
    generateU!,
    generateU,
    generateN,
    countN,
    FlipKet,
    subregion,
    hopLeft,
    hopRight,
    someFunction!,
    destroyParticle,
    createParticle,
    TunnelParticle,
    emptyToFirst,
    shiftVector,
    reverseBasis,
    ReverseKet,
    shiftV,
    Get_ket_leader,
    Get_ket_leader_OTF,
    construct_OBDM,
    construct_OBDM_OTF,
buildpascal,
in2b,
b2in,

    Translational_Symmetry_Reflection_Cycles_gs,
    Translational_Symmetry_Reflection_Cycles_gs_reduced,
    Translational_Symmetry_Reflection_Cycles_gs_reducedb,   
    Translational_Symmetry_Reflection_Cycles_gs_OTF,
    Translational_Symmetry_Reflection_Cycles_gs_OTFb,
    sparse_Block_Diagonal_Hamiltonian_q0R1_gs,
    sparse_Block_Diagonal_Hamiltonian_q0R1_gs_reduced,
    sparse_Block_Diagonal_Hamiltonian_q0R1_gs_OTF,
    Symmetry_Reflection_Cycles_gs,
    serial_num_Cycles,
    serial_num_Cycles_OTF

include("utilities.jl")
include("translational_reflection_symmetry_gs.jl")

end
