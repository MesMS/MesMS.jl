import CSV

"`input` accept the same types as `CSV.File` including file path and `IO`"
read_precursor(input) = begin
    I = Dict{Int, Vector{Ion}}()
    for r in CSV.File(input)
        push!(get!(I, r.scan, Ion[]), Ion(r.mz, r.z))
    end
    return I
end

"`input` accept the same types as `CSV.File` including file path and `IO`"
read_cross_link_pair(input) = begin
    T_v = typeof((idx=0, mz=0.0, z=0, pairs=[(0.0, 0.0)]))
    P = Dict{Int, Vector{T_v}}()
    for r in CSV.File(input)
        v = get!(P, r.scan, T_v[])
        i = findfirst(x -> x.idx == r.prec, v)
        if isnothing(i)
            push!(v, (; idx=r.prec, r.mz, r.z, pairs=typeof((0.0, 0.0))[]))
            i = length(v)
        end
        push!(v[i].pairs, (r.a, r.b))
    end
    return P
end
