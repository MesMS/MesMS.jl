import CSV

read_mes(io::IO; split=true) = begin
    format = read!(io, Array{UInt8}(undef, 4))
    if format != b"MES\n" error("invalid .mes data") end
    version = read(io, UInt32)
    if version != 0 error("unsupported .mes data version") end
    n = read(io, UInt64)
    M = map(_ -> read!(io, Array{Float64}(undef, read(io, UInt64))), 1:n)
    I = map(_ -> read!(io, Array{Float64}(undef, read(io, UInt64))), 1:n)
    P = map(((m, i),) -> Peak.(m, i), zip(M, I))
    s = read!(io, Array{UInt8}(undef, 2^13))
    if any(s .!= UInt8('\n')) error("invalid .mes data") end
    meta = CSV.File(io)
    st = meta.ScanType .|> Symbol
    id = meta.ScanID
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
                retention_time=rt[i], injection_time=it[i], peaks=P[i],
            )
        elseif st[i] == :MS2
            return MS2(; id=id[i], pre=pre[i],
                total_ion_current=tic[i], base_peak_intensity=bpi[i], base_peak_mass=bpm[i],
                retention_time=rt[i], injection_time=it[i], peaks=P[i],
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

write_mes(io::IO, M) = begin
    write(io, b"MES\n", UInt32(0), UInt64(length(M)))
    foreach(m -> write(io, UInt64(length(m.peaks)), getfield.(m.peaks, :mz)), M)
    foreach(m -> write(io, UInt64(length(m.peaks)), getfield.(m.peaks, :inten)), M)
    write(io, fill(UInt8('\n'), 2^13))
    CSV.write(io, [(;
            ScanType=string(m.type),
            ScanID=m.id,
            PrecursorScan= typeof(m) == MS2 ? m.pre : 0,
            TotalIonCurrent=m.total_ion_current,
            BasePeakIntensity=m.base_peak_intensity,
            BasePeakMass=m.base_peak_mass,
            RetentionTime=m.retention_time,
            IonInjectionTime=m.injection_time,
            ActivationCenter=typeof(m) == MS2 ? m.activation_center : 0.0,
            IsolationWidth=typeof(m) == MS2 ? m.isolation_width : 0.0,
            PrecursorMZ=typeof(m) == MS2 ? m.ions[1].mz : 0.0,
            PrecursorCharge=typeof(m) == MS2 ? m.ions[1].z : 0,
        ) for m in M]; append=true, writeheader=true,
    )
end

read_mes(path::AbstractString; split=true, verbose=true) = begin
    verbose && @info "MES reading from " * path
    return open(io -> read_mes(io; split), path)
end

write_mes(path::AbstractString, M) = open(io -> write_mes(io, M), path; write=true)
