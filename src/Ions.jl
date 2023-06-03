abstract type AbstractIon end

Base.isless(a::AbstractIon, b::AbstractIon) = a.z < b.z || a.z == b.z && a.mz < b.mz

struct Ion <: AbstractIon
    mz::Float64
    z::Int
    m::Float64
    Ion(mz, z) = new(mz, z, mz * z)
end
