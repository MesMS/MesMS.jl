import CSV

read_peak(io::IO) = begin
    M = map(_ -> read!(io, Array{Float64}(undef, read(io, Int64))), 1:read(io, Int64))
    I = map(_ -> read!(io, Array{Float64}(undef, read(io, Int64))), 1:read(io, Int64))
    return map(((m, i),) -> Peak.(m, i), zip(M, I))
end

write_peak(io::IO, P::AbstractArray{<:AbstractArray{<:AbstractPeak}}) = begin
    write(io, Int64(length(P)))
    foreach(ps -> write(io, Int64(length(ps)), getfield.(ps, :mz)), P)
    write(io, Int64(length(P)))
    foreach(ps -> write(io, Int64(length(ps)), getfield.(ps, :inten)), P)
end

write_peak(io::IO, M::AbstractArray{<:AbstractMS}) = begin
    write_peak(io::IO, getfield.(M, :peaks))
end

read_peak(path::AbstractString; verbose=true) = begin
    verbose && @info "peak list reading from " * path
    return open(io -> read_peak(io), path)
end

write_peak(path::AbstractString, X) = open(io -> write_peak(io, X), path; write=true)

read_ms(path; split=true, verbose=true) = begin
    P = read_peak(path; verbose)
    verbose && @info "scan meta reading from " * splitext(path)[1] * ".scan.csv"
    meta = CSV.File(splitext(path)[1] * ".scan.csv")
    n = length(P)
    st = meta.ScanType .|> Symbol
    id = meta.ScanID
    peak = hasproperty(meta, :PeakIndex) ? meta.PeakIndex : Vector(1:n)
    pre = hasproperty(meta, :PrecursorScan) ? meta.PrecursorScan : zeros(Int, n)
    tic = hasproperty(meta, :TotalIonCurrent) ? meta.TotalIonCurrent : zeros(Float64, n)
    bpi = hasproperty(meta, :BasePeakIntensity) ? meta.BasePeakIntensity : zeros(Float64, n)
    bpm = hasproperty(meta, :BasePeakMass) ? meta.BasePeakMass : zeros(Float64, n)
    rt = hasproperty(meta, :RetentionTime) ? meta.RetentionTime : zeros(Float64, n)
    it = hasproperty(meta, :IonInjectionTime) ? meta.IonInjectionTime : zeros(Float64, n)
    ac = hasproperty(meta, :ActivationCenter) ? meta.ActivationCenter : zeros(Float64, n)
    iw = hasproperty(meta, :IsolationWidth) ? meta.IsolationWidth : zeros(Float64, n)
    mz = hasproperty(meta, :PrecursorMZ) ? meta.PrecursorMZ : zeros(Float64, n)
    z = hasproperty(meta, :PrecursorCharge) ? meta.PrecursorCharge : zeros(Int, n)
    ions = [z[i] != 0 && mz[i] != 0 ? [Ion(mz[i], z[i])] : Ion[] for i in 1:n]
    M = map(1:n) do i
        if st[i] == :MS1
            return MS1(; id=id[i],
                total_ion_current=tic[i], base_peak_intensity=bpi[i], base_peak_mass=bpm[i],
                retention_time=rt[i], injection_time=it[i], peaks=P[peak[i]],
            )
        elseif st[i] == :MS2
            return MS2(; id=id[i], pre=pre[i],
                total_ion_current=tic[i], base_peak_intensity=bpi[i], base_peak_mass=bpm[i],
                retention_time=rt[i], injection_time=it[i], peaks=P[peak[i]],
                activation_center=ac[i], isolation_width=iw[i],
                ions=ions[i],
            )
        else
            @warn "skipping unknown scan type: " * string(st[i])
        end
    end
    if split
        return (;
            MS1=filter(m -> typeof(m) == MS1, M) |> Vector{MS1},
            MS2=filter(m -> typeof(m) == MS2, M) |> Vector{MS2},
        )
    else
        return M
    end
end

write_scanmeta(io::IO, M) = begin
    write(io,
        "ScanType,ScanID,PeakIndex,TotalIonCurrent,BasePeakIntensity,BasePeakMass,\
        RetentionTime,IonInjectionTime,\
        PrecursorScan,ActivationCenter,IsolationWidth,PrecursorMZ,PrecursorCharge\n"
    )
    for (i, m) in enumerate(M)
        if typeof(m) == MS1
            write(io,
                "MS1,$(m.id),$(i),$(m.total_ion_current),$(m.base_peak_intensity),$(m.base_peak_mass),\
                $(m.retention_time),$(m.injection_time),\
                0,0.0,0.0,0.0,0\n"
            )
        elseif typeof(m) == MS2
            write(io,
                "MS2,$(m.id),$(i),$(m.total_ion_current),$(m.base_peak_intensity),$(m.base_peak_mass),\
                $(m.retention_time),$(m.injection_time),\
                $(m.pre),$(m.activation_center),$(m.isolation_width),\
                $(join([i.mz for i in m.ions], '|')),$(join([i.z for i in m.ions], '|'))\n"
            )
        end
    end
end

write_scanmeta(path::AbstractString, M) = open(io -> write_scanmeta(io, M), path; write=true)
