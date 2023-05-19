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

write_precursor(path, ions) = begin
    ions = map(i -> (; i.scan, i.mz, i.z), ions)
    @info "precursor ion saving to " * path
    CSV.write(path * "~", ions)
    mv(path * "~", path; force=true)
end

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
    rows = []
    for ion in ions, pair in ion.pairs
        push!(rows, (; ion.scan, ion.xl, ion.mz, ion.z, a=pair[1], b=pair[2]))
    end
    @info "cross-link pair saving to " * path
    CSV.write(path * "~", rows)
    mv(path * "~", path; force=true)
end
