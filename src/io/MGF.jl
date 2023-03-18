using Base: Filesystem

read_mgf(io::IO) = begin
    M = MS2[]
    mz = 0.0
    z = 0
    peaks = Peak[]
    while !eof(io)
        line = readline(io)
        if length(line) == 0
            continue
        elseif line == "BEGIN IONS"
            mz = 0.0
            z = 0
            peaks = Peak[]
        elseif line == "END IONS"
            push!(M, MS2(; ions=[Ion(mz, z)], peaks))
        elseif startswith(line, "PEPMASS=")
            mz = parse(Float64, line[9:end])
        elseif startswith(line, "CHARGE=")
            z = line[end] == '+' ? parse(Int, line[8:end-1]) : -parse(Int, line[8:end-1])
        elseif !occursin('=', line)
            m, i = split(line)
            push!(peaks, Peak(parse(Float64, m), parse(Float64, i)))
        end
    end
    return M
end

read_mgf(fname::AbstractString) = open(read_mgf, fname)

write_mgf(io::IO, m::AbstractTandemMS, title="$(m.id)") = begin
    for ion in m.ions
        write(io, "BEGIN IONS\n")
        write(io, "TITLE=$(title)\n")
        write(io, "PEPMASS=$(ion.mz)\n")
        write(io, "CHARGE=$(abs(ion.z))$(ion.z >= 0 ? '+' : '-')\n")
        foreach(p -> write(io, "$(p.mz) $(p.inten)\n"), m.peaks)
        write(io, "END IONS\n")
    end
end

write_mgf(io::IO, M::AbstractArray{<:AbstractTandemMS}, title=m -> "$(m.id)") = foreach(m -> write_mgf(io, m, title(m)), M)
write_mgf(fname::AbstractString, args...) = open(io -> write_mgf(io, args...), fname; write=true)
