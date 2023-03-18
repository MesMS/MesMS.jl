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

read_fasta(fname::AbstractString; itol=true) = open(io -> read_fasta(io; itol), fname)
