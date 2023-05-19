using Base: Filesystem

read_ms1(io::IO) = begin
    M = MS1[]
    id = 0
    total_ion_current = 0.0
    base_peak_intensity = 0.0
    base_peak_mass = 0.0
    retention_time = 0.0
    injection_time = 0.0
    peaks = Peak[]
    while !eof(io)
        line = readline(io)
        if length(line) == 0
            continue
        elseif line[1] == 'H'
            continue
        elseif line[1] == 'S'
            !isempty(peaks) && push!(M, MS1(;
                id, total_ion_current, base_peak_intensity, base_peak_mass, retention_time, injection_time, peaks,
            ))
            id = parse(Int, split(line[3:end])[1])
            total_ion_current = 0.0
            base_peak_intensity = 0.0
            base_peak_mass = 0.0
            retention_time = 0.0
            injection_time = 0.0
            peaks = Peak[]
        elseif line[1] == 'I'
            k, v = split(line[3:end])
            if k == "TotalIonCurrent"
                total_ion_current = parse(Float64, v)
            elseif k == "BasePeakIntensity"
                base_peak_intensity = parse(Float64, v)
            elseif k == "BasePeakMass"
                base_peak_mass = parse(Float64, v)
            elseif k == "RetentionTime" || k == "RetTime"  || k == "RTime"
                retention_time = parse(Float64, v)
            elseif k == "IonInjectionTime" || k == "InjectionTime"
                injection_time = parse(Float64, v)
            end
        elseif !isspace(line[2])
            m, i = split(line)
            push!(peaks, Peak(parse(Float64, m), parse(Float64, i)))
        end
    end
    !isempty(peaks) && push!(M, MS1(;
        id, total_ion_current, base_peak_intensity, base_peak_mass, retention_time, injection_time, peaks,
    ))
    return M
end

read_ms2(io::IO) = begin
    M = MS2[]
    id = 0
    pre = 0
    total_ion_current = 0.0
    base_peak_intensity = 0.0
    base_peak_mass = 0.0
    retention_time = 0.0
    injection_time = 0.0
    activation_center = 0.0
    isolation_width = 0.0
    ions = Ion[]
    peaks = Peak[]
    while !eof(io)
        line = readline(io)
        if length(line) == 0
            continue
        elseif line[1] == 'H'
            continue
        elseif line[1] == 'S'
            !isempty(peaks) && push!(M, MS2(; id, pre,
                total_ion_current, base_peak_intensity, base_peak_mass,
                retention_time, injection_time, activation_center, isolation_width, ions, peaks,
            ))
            items = split(line[3:end])
            id = parse(Int, items[1])
            pre = 0
            total_ion_current = 0.0
            base_peak_intensity = 0.0
            base_peak_mass = 0.0
            retention_time = 0.0
            injection_time = 0.0
            activation_center = length(items) >= 3 ? parse(Float64, items[3]) : 0.0
            isolation_width = 0.0
            ions = Ion[]
            peaks = Peak[]
        elseif line[1] == 'I'
            k, v = split(line[3:end])
            if k == "PrecursorScan"
                pre = parse(Int, v)
            elseif k == "TotalIonCurrent"
                total_ion_current = parse(Float64, v)
            elseif k == "BasePeakIntensity"
                base_peak_intensity = parse(Float64, v)
            elseif k == "BasePeakMass"
                base_peak_mass = parse(Float64, v)
            elseif k == "RetentionTime" || k == "RetTime"  || k == "RTime"
                retention_time = parse(Float64, v)
            elseif k == "IonInjectionTime" || k == "InjectionTime"
                injection_time = parse(Float64, v)
            elseif k == "ActivationCenter"
                activation_center = parse(Float64, v)
            elseif k == "IsolationWidth"
                isolation_width = parse(Float64, v)
            end
        elseif line[1] == 'Z'
            z_, mh_ = split(line[3:end])
            z = parse(Int, z_)
            mz = mh_to_mz(parse(Float64, mh_), z)
            push!(ions, Ion(mz, z))
        elseif !isspace(line[2])
            m, i = split(line)
            push!(peaks, Peak(parse(Float64, m), parse(Float64, i)))
        end
    end
    !isempty(peaks) && push!(M, MS2(; id, pre,
        total_ion_current, base_peak_intensity, base_peak_mass,
        retention_time, injection_time, activation_center, isolation_width, ions, peaks,
    ))
    return M
end

read_ms1(fname::AbstractString; verbose=true) = begin
    verbose && @info "MS1 reading from " * fname
    return open(read_ms1, fname)
end

read_ms2(fname::AbstractString; verbose=true) = begin
    verbose && @info "MS2 reading from " * fname
    return open(read_ms2, fname)
end

write_ms1(io::IO, m::AbstractMS) = begin
    write(io, "S\t$(m.id)\t$(m.id)\n")
    write(io, "I\tTotalIonCurrent\t$(m.total_ion_current)\n")
    write(io, "I\tBasePeakIntensity\t$(m.base_peak_intensity)\n")
    write(io, "I\tBasePeakMass\t$(m.base_peak_mass)\n")
    write(io, "I\tRetentionTime\t$(m.retention_time)\n")
    write(io, "I\tIonInjectionTime\t$(m.injection_time)\n")
    foreach(p -> write(io, "$(p.mz) $(p.inten)\n"), m.peaks)
end

write_ms2(io::IO, m::AbstractTandemMS) = begin
    write(io, "S\t$(m.id)\t$(m.id)\t$(m.activation_center)\n")
    write(io, "I\tPrecursorScan\t$(m.pre)\n")
    write(io, "I\tTotalIonCurrent\t$(m.total_ion_current)\n")
    write(io, "I\tBasePeakIntensity\t$(m.base_peak_intensity)\n")
    write(io, "I\tBasePeakMass\t$(m.base_peak_mass)\n")
    write(io, "I\tRetentionTime\t$(m.retention_time)\n")
    write(io, "I\tIonInjectionTime\t$(m.injection_time)\n")
    write(io, "I\tActivationCenter\t$(m.activation_center)\n")
    foreach(i -> write(io, "Z\t$(i.z)\t$(mz_to_mh(i.mz, i.z))\n"), m.ions)
    foreach(p -> write(io, "$(p.mz) $(p.inten)\n"), m.peaks)
end

write_ms1(io::IO, M::AbstractArray{<:AbstractMS}) = foreach(m -> write_ms1(io, m), M)
write_ms2(io::IO, M::AbstractArray{<:AbstractTandemMS}) = foreach(m -> write_ms2(io, m), M)
write_ms1(fname::AbstractString, args...) = open(io -> write_ms1(io, args...), fname; write=true)
write_ms2(fname::AbstractString, args...) = open(io -> write_ms2(io, args...), fname; write=true)

count_msx(path, ext; verbose=true) = begin
    verbose && @info "counting $(path)"
    s = sum(match_path(path, ext)) do f
        n = count(l -> startswith(l, "S\t"), readlines(f))
        verbose && @info "$(basename(f)):\t$(n)"
        return n
    end
    verbose && @info "total:\t$(s)"
    return s
end
