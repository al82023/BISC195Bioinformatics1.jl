module Assignment07

export normalizeDNA,
        composition,
        gc_content

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
        # add 1 to each base as it occurs
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
For more info about GC content, see here:
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

end # module Assignment07