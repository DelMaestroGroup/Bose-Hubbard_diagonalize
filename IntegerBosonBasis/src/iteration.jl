## Occupation vector iteration and array properties.

function Base.iterate(basis::AbstractIntbasis, state=1)
    state > basis.D && return nothing
    basis.vectors[state], state+1
end

Base.eltype(::Type{AbstractIntbasis}) = Int64
Base.length(basis::AbstractIntbasis) = basis.D

function Base.in(v::Int64, basis::Intbasis)
    count_ones(v) == basis.L && count_zeros(v) - leading_zeros(v) == basis.N
end

Base.getindex(basis::AbstractIntbasis, i::Int) = basis.vectors[i]
