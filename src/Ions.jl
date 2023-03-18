abstract type AbstractIon end

Base.isless(a::AbstractIon, b::AbstractIon) = a.z < b.z || a.z == b.z && a.mz < b.mz

struct Ion <: AbstractIon
    m::Float64
    z::Int
    mz::Float64
    Ion(mz, z) = new(mz * z, z, mz)
    Ion(m, z, mz) = new(m, z, mz)
end
