module BISC195Bioinformatics1

export normalizeDNA,
        composition,
        gc_content,
        complement,
        reverse_complement,
        parse_fasta

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("Invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end


# Your code here.
# Don't forget to export your functions!

"""
    composition(sequence)

Counts the number of each type of base (including 'N')
in a DNA sequence and returns a Dict

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> composition("ACCGGGTTTTN")
    Dict{Char,Int64} with 5 entries:
        'A' => 1
        'G' => 3
        'T' => 4
        'N' => 1
        'C' => 2
        
    julia> composition("AAX")
    ERROR: Invalid base, X
"""
function composition(sequence)
    sequence = normalizeDNA(sequence)
    comp = Dict('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0, 'N' => 0)
    for base in sequence
        if base == 'A'
            comp['A'] = comp['A'] + 1
        elseif base == 'C'
            comp['C'] = comp['C'] + 1
        elseif base == 'G'
            comp['G'] = comp['G'] + 1
        elseif base == 'T'
            comp['T'] = comp['T'] + 1
        elseif base == 'N'
            comp['N'] = comp['N'] + 1
        end
    end
    return comp
end

"""
    gc_content(sequence)

Calculates the GC ratio of a DNA sequence.
The GC ratio is the total number of G and C bases divided by the total length of the sequence.

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> gc_content("AATG")
    0.25

    julia> gc_content("cccggg") * 100
    100.0

    julia> gc_content("ATta")
    0.0

    julia> gc_content("ATty")
    Error: Invalid base Y encountered

    julia> gc_content("ATCGN")
    0.4
"""
function gc_content(sequence)
    basecomp = composition(sequence)
    return (basecomp['G'] + basecomp['C']) / length(sequence)
end

"""
    complement(base::Char)

Get the DNA complement of the provided base:

    A <-> T
    G <-> C

Accepts uppercase or lowercase `String` or `Char`,
but always returns an uppercase `Char`.
If a valid base is not provided, the function throws an error.

Examples
≡≡≡≡≡≡≡≡≡≡

    julia> complement('C')
    'G': ASCII/Unicode U+0047 (category Lu: Letter, uppercase)
"""
function complement(base::Char)
    complements = Dict("A" => 'T',
                       "T" => 'A',
                       "G" => 'C',
                       "C" => 'G',
                       "N" => 'N')
    base = uppercase(string(base))
    !(base in keys(complements)) && error("Invalid base $base")
    return complements[base]
end

"""
    complement(sequence::AbstractString)

Takes a DNA sequence and returns the complement
of that sequence.

Takes lowercase or uppercase sequences,
but always returns uppercase.

Examples
≡≡≡≡≡≡≡≡≡≡

    julia> complement("ATTN")
    "TAAN"

    julia> complement("ATTAGC")
    "TAATCG"
"""
function complement(sequence::AbstractString)
    seq = normalizeDNA(sequence)
    comp = map(complement, seq)
    return comp
end

"""
    reverse_complement(sequence)

Takes a DNA sequence and returns the reverse complement
of that sequence.

Takes lowercase or uppercase sequences,
but always returns uppercase.

Examples
≡≡≡≡≡≡≡≡≡≡

    julia> reverse_complement("AAATTT")
    "AAATTT"

    julia> reverse_complement("GCAT")
    "ATGC"

    julia> rc = reverse_complement("TTGGG");

    julia> println(rc)
    CCCAA

    julia> reverse_complement("ATTAGC")
    "GCTAAT"

    julia> reverse_complement("ATN")
    "NAT"
"""
function reverse_complement(sequence)
    seq = normalizeDNA(sequence)
    comp = map(complement, seq)
    reversecomp = reverse(comp)
    return reversecomp
end

"""
    function parse_fasta(path)

Reads a fasta-formated file and returns 2 vectors,
one containing headers,
the other containing the entire sequence as a `String`.

Examples
≡≡≡≡≡≡≡≡≡

    julia> ex1 = parse_fasta("data/ex1.fasta");

    julia> ex1[1]
    2-element Array{String,1}:
    "ex1.1 | easy"
    "ex1.2 | multiline"

    julia> ex1[2]
    2-element Array{String,1}:
    "AATTATAGC"
    "CGCCCCCCAGTCGGATT"

    julia> ex2 = parse_fasta("data/ex2.fasta");
    ERROR: invalid base H
"""
function parse_fasta(path)
    headers = String[]
    sequences = String[]
    sequence = []
    for line in eachline(path)
        if occursin('>', line)
            push!(headers, chop(line, head = 1, tail = 0))
            push!(sequences, join(sequence))
            sequence = []
        else
            push!(sequence, line)
        end
    end
    push!(sequences, join(sequence))
    popfirst!(sequences)
    normalizeDNA.(sequences)
    return (headers, sequences)
end

end # module BISC195Bioinformatics1