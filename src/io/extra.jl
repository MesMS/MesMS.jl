import CSV

"`input` accept the same types as `CSV.File` including file path and `IO`"
read_precursor(input) = begin
    I = Dict{Int, Vector{Ion}}()
    for r in CSV.File(input)
        push!(get!(I, r.scan, Ion[]), Ion(r.mz, r.z))
    end
    return I
end
