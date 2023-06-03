abstract type AbstractIon end

Base.isless(a::AbstractIon, b::AbstractIon) = a.z < b.z || a.z == b.z && a.mz < b.mz

struct Ion <: AbstractIon
    mz::Float64
    z::Int
    m::Float64
    Ion(mz, z) = new(mz, z, mz * z)
end

Base.show(io::IO, i::Ion) = begin
    if i.z > 0 write(io, "(m/z=$(i.mz))$(i.z == 1 ? "" : supscript(string.(i.z)))⁺")
    elseif i.z < 0 write(io, "(m/z=$(i.mz))$(i.z == -1 ? "" : supscript(string.(-i.z)))⁻")
    else write(io, "(m/z=$(i.mz))⁰")
    end
    return nothing
end
