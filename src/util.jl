using Dates

"convert a char to its subscript form"
subscript(c::Char) = begin
    if '0' ≤ c ≤ '9' return Char(0x2080 + c - '0')
    elseif c == '+' return '₊'
    elseif c == '-' return '₋'
    elseif c == '=' return '₌'
    elseif c == '(' return '₍'
    elseif c == ')' return '₎'
    else return c
    end
end
subscript(s::AbstractString) = join(subscript.(collect(s)))

"convert a char to its supscript form"
supscript(c::Char) = begin
    if c == '0' return '⁰'
    elseif c == '1' return '¹'
    elseif c == '2' return '²'
    elseif c == '3' return '³'
    elseif c == '4' return '⁴'
    elseif c == '5' return '⁵'
    elseif c == '6' return '⁶'
    elseif c == '7' return '⁷'
    elseif c == '8' return '⁸'
    elseif c == '9' return '⁹'
    elseif c == '+' return '⁺'
    elseif c == '-' return '⁻'
    elseif c == '=' return '⁼'
    elseif c == '(' return '⁽'
    elseif c == ')' return '⁾'
    else return c
    end
end
supscript(s::AbstractString) = join(supscript.(collect(s)))

split(str::AbstractString, args...; kwargs...) = Base.split(str, args...; kwargs...)
split(::Missing, args...; kwargs...) = SubString{AbstractString}[]

match_path(path, ext=""::Union{AbstractString, AbstractArray{AbstractString}};
    join=true, self=true,
) = begin
    if isa(ext, AbstractString) ext = [ext] end
    files = filter(readdir(dirname(path))) do f
        (self && (f == basename(path))) ||
        (startswith(f, basename(path)) && (isempty(ext) || any(e -> endswith(f, e), ext)))
    end
    return join ? joinpath.(dirname(path), files) : files
end

read_all(reader, path, ext=""; verbose=true, single=false, key=first∘splitext∘basename) = begin
    paths = match_path(path, ext)
    if single && length(paths) > 1
        @warn "multiple files found: path=$(path), ext=$(ext)"
    end
    return map(paths) do p
        verbose && @info "reading from " * p
        key(p) => reader(p)
    end |> Dict
end

parse_range(::Type{T}, s::AbstractString) where T = begin
    vs = map(x -> parse(T,  x),  split(s, ':'))
    if length(vs) <= 2 return range(vs[begin], vs[end])
    elseif length(vs) == 3 return range(vs[begin], vs[end]; step=vs[begin+1])
    else return nothing
    end
end

date_mark(date=today()) = replace(string(date), "-"=>"")

open_url(url; verbose=true) = begin
    verbose && @info "opening " * url
    if Sys.isapple()
        run(`open $(url)`)
    elseif Sys.iswindows()
        run(`cmd /k start $(url) "&&" exit`)
    else
        run(`open $(url)`)
    end
end
