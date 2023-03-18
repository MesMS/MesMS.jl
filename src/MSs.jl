abstract type AbstractMS end
abstract type AbstractTandemMS <: AbstractMS end

Base.isless(a::AbstractMS, b::AbstractMS) = a.id < b.id

Base.@kwdef struct MS1 <: AbstractMS
    id::Int = 0
    retention_time::Float64 = 0.0
    injection_time::Float64 = 0.0
    peaks::Vector{<:AbstractPeak} = Peak[]
end

Base.@kwdef struct MS2 <: AbstractTandemMS
    id::Int = 0
    pre::Int = 0
    retention_time::Float64 = 0.0
    injection_time::Float64 = 0.0
    activation_center::Float64 = 0.0
    isolation_width::Float64 = 0.0
    ions::Vector{<:AbstractIon} = Ion[]
    peaks::Vector{<:AbstractPeak} = Peak[]
end

dict_by_id(M::AbstractArray{<:AbstractMS}) = Dict([m.id => m for m in M])
