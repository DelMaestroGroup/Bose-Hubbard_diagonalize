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

    Translational_Symmetry_Reflection_Cycles_gs,
    sparse_Block_Diagonal_Hamiltonian_q0R1_gs,
    sparse_Block_Diagonal_Hamiltonian_q0R1_gs_Rest,
    Symmetry_Reflection_Cycles_gs,
    spatial_entropy_acc_gs,
    spatial_entropy_acc_gs_Rest,
    particle_entropy

include("utilities.jl")
include("translational_reflection_symmetry_gs.jl")
include("reflection_symmetry_gs.jl")
include("spatial_entropy_gs.jl")
include("particle_entropy.jl")

end
