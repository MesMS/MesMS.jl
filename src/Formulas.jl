abstract type AbstractFormula end

struct Formula <: AbstractFormula
    data::NamedTuple
end

Formula(; kwargs...) = Formula(NamedTuple(kwargs))

Base.propertynames(f::AbstractFormula, private::Bool=false) = propertynames(getfield(f, :data), private)
Base.getproperty(f::AbstractFormula, name::Symbol) = getproperty(getfield(f, :data), name)
Base.getindex(f::AbstractFormula, name::Symbol) = getproperty(f, name)

Base.parse(T, ::Type{Formula}, s::AbstractString) = begin
    es = map(e -> split(e, '('), split(s, ')'; keepempty=false))
    es = map(((e, n),) -> (e=Symbol(e), n=parse(T, n)), es)
    d = Dict{Symbol, T}()
    foreach(x -> d[x.e] = get(d, x.e, 0) + x.n, es)
    return Formula(; d...)
end

Base.parse(::Type{Formula}, s::AbstractString) = parse(Int, Formula, s)

Base.show(io::IO, f::MesMS.Formula) = begin
    write(io, join(map(p -> "$(p)$(subscript.(string(getproperty(f, p))))", propertynames(f))))
    return nothing
end

mass(f::Formula, tab::NamedTuple) = sum(f[n] * tab[n] for n in propertynames(f); init=0.0)
