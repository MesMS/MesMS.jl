unify_aa_seq(seq::AbstractString; itol=true) = strip(itol ? replace(seq, 'I' => 'L') : seq)
unify_aa_seq(seq::Missing; itol=true) = ""

read_fasta(io::IO; itol=true) = begin
    seqs = []
    for line in readlines(io)
        if startswith(line, '>') push!(seqs, (; name=strip(line[2:end]), lines=[]))
        else push!(seqs[end].lines, unify_aa_seq(line; itol))
        end
    end
    return Dict{String, String}(map(s -> (s.name, join(s.lines)), seqs))
end

read_fasta(fname::AbstractString; verbose=true, itol=true) = begin
    verbose && @info "FASTA reading from " * fname
    return open(io -> read_fasta(io; itol), fname)
end

write_fasta(io::IO, S::Dict; n=60, itol=true) = begin
    for (k, v) in S
        v = unify_aa_seq(v; itol)
        write(io, ">$k\n")
        for i in 1:n:length(v)
            write(io, v[i:min(i + n - 1, end)])
            write(io, "\n")
        end
    end
    return io
end

write_fasta(fname::AbstractString, S::Dict; n=60, itol=true) = open(io -> write_fasta(io, S; n, itol), fname; write=true)
