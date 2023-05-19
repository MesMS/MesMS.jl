import CSV

"`input` accept the same types as `CSV.File` including file path and `IO`"
read_precursor(input; verbose=true) = begin
    verbose && isa(input, AbstractString) && @info "precursor ion reading from " * input
    I = Dict{Int, Vector{Ion}}()
    for r in CSV.File(input)
        push!(get!(I, r.scan, Ion[]), Ion(r.mz, r.z))
    end
    return I
end

write_precursor(path, ions) = CSV.write(path, map(i -> (; i.scan, i.mz, i.z), ions))

"`input` accept the same types as `CSV.File` including file path and `IO`"
read_cross_link_pair(input; verbose=true) = begin
    verbose && isa(input, AbstractString) && @info "cross link pair reading from " * input
    T_v = typeof((scan=0, xl=0, mz=0.0, z=0, pairs=[(0.0, 0.0)]))
    P = Dict{Int, Vector{T_v}}()
    for r in CSV.File(input)
        ions = get!(P, r.scan, T_v[])
        i = findfirst(x -> x.xl == r.xl, ions)
        if isnothing(i)
            push!(ions, (; r.scan, r.xl, r.mz, r.z, pairs=typeof((0.0, 0.0))[]))
            i = length(ions)
        end
        push!(ions[i].pairs, (r.a, r.b))
    end
    return P
end

write_cross_link_pair(path, ions) = begin
    CSV.write(path, [(; i.scan, i.xl, i.mz, i.z, a=p[1], b=p[2]) for i in ions for p in i.pairs])
end
