# title          :lac_lib_short.py
# description    :Generates 6400 unique inducible lac promoters with varied operator sites
# author         :timcyu
# date           :2/23/18
# =============================================================================================
def reverse_complement(seq):
    """
    Return the reverse complement of a nucleotide string
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 't': 'a', 'a': 't', 'c': 'g', 'g': 'c'}
    rc = ''.join([complement[nt] for nt in seq[::-1]])
    return rc

def best_A_content(oligo):
    '''
    Choose the strand with the lowest A content because A's are harder to
        synthesize.
    '''
    rc_oligo = reverse_complement(oligo)

    oligo_As = sum( [1 for nt in oligo if nt == 'A'] )
    rc_As = sum( [1 for nt in rc_oligo if nt == 'A'] )

    if oligo_As < rc_As:
            final_oligo = oligo
    else:
            final_oligo = rc_oligo

    return final_oligo

if __name__ == '__main__':
    # Upstream and downstream background sequences
    five_prime_bg = "TCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCG"
    three_prime_bg = "TCACACAGGAAACAGCTATGA"

    # Core promoter sequences (contains 4 different -10/-35 element pairs)
    core_lacUV5 = "AATGTAAGTTAGCTCATTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGG"
    core_1 = "AATGTAAGTTAGCTCATTCATTAGGCACCCCAGGCTTTTACATTTATGCTTCCGGCTCGTATAAAGTGTGG"
    core_2 = "AATGTAAGTTAGCTCATTCATTAGGCACCCCAGGCTTTGCAATTTATGCTTCCGGCTCGAATAATGTGTGG"
    core_3 = "AATGTAAGTTAGCTCATTCATTAGGCACCCCAGGCTTTGACATTTATGCTTCCGGCTCGGATAATGTGTGG"

    # Lac operator site sequences
    lacO_1 = "AATTGTGAGCGGATAACAATT"
    lacO_2 = "AAATTGTAGCGAGTAACAACC"
    lacO_3 = "GGCAGTGAGCGCAACGCAATT"
    lacO_sym = "AAATTGTGAGCGCTCACAATT"

    # Store the three regions in lists for convenience
    bg = [five_prime_bg, three_prime_bg]
    core = [core_lacUV5, core_1, core_2, core_3]
    operator = [lacO_1, lacO_2, lacO_3, lacO_sym]

    # Store their names in lists for naming
    core_name = ["lacUV5", "core1", "core2", "core3"]
    operator_name = ["O1", "O2", "O3", "Osym"]

    # Generates a list of all possible 5' upstream elements with 1 operator (16 total) and 2 operators (64 total).
    # These operators shift loc base pairs upstream relative to the 5' end of the core promoter sequence.
    loc = [7, 9, 18, 22]
    up1 = []
    name_up1 = []
    up2 = []
    name_up2 = []
    temp = ""

    for i in loc:
        for j, k in zip(operator, operator_name):
            temp = five_prime_bg[::-1][0:i] + j[::-1] + five_prime_bg[::-1][i:43]
            up1.append(temp[::-1])
            name_up1.append("lacI_" + k + "_offset" + str(i) + "_")

    for i in loc:
        for j, m in zip(operator, operator_name):
            for k, n in zip(operator, operator_name):
                temp = five_prime_bg[::-1][0:i] + j[::-1] + k[::-1] + five_prime_bg[::-1][i + 21:43]
                up2.append(temp[::-1])
                name_up2.append("lacI_" + n + "_" + m + "_offset" + str(i) + "_")

    # These lists contain all 80 upstream elements and their names
    totUpstream = up1 + up2
    totUpstream_names = name_up1 + name_up2

    # Generates all possible 3' downstream elements with 1 operator (4 total)
    # Also generates all scenarios with 2 operators (16 total)
    down1 = []
    name_down1 = []
    down2 = []
    name_down2 = []

    for i, j in zip(operator, operator_name):
        down1.append(i + three_prime_bg)
        name_down1.append("_" + j)

    for i, m in zip(operator, operator_name):
        for j, n in zip(operator, operator_name):
            down2.append(i + j)
            name_down2.append("_" + m + "_" + n)

    # These lists contain all 20 downstream elements and their names
    totDownstream = down1 + down2
    totDownstream_names = name_down1 + name_down2

    # Generates temporary list containing combinations of cores with downstream elements
    temp = []
    temp_names = []
    for i, m in zip(core, core_name):
        for j, n in zip(totDownstream, totDownstream_names):
            temp.append(i + j)
            temp_names.append(m + n)

    # Generates full list containing library of all possible upstream, core, and downstream combinations
    library = []
    library_names = []

    # Forward and Reverse primers for cloning purposes. Restriction sites included
    fwd_primer = "ACCTGTAATTCCAAGCGTCTCGAG"
    rev_primer = "GCTAGCGGTGTTTAGTTAGCATCC"

    count = 0
    for i, m in zip(totUpstream, totUpstream_names):
        for j, n in zip(temp, temp_names):
            library.append(best_A_content(fwd_primer + i + j + rev_primer))
            library_names.append(m + n)
            count = count + 1

    print ("This library contains " + str(count) + " unique sequences.")

    for i in library:
        print len(i)

    # Creates a fasta format file containing all sequences and names
    file = open("lac_lib_short.txt", "w")
    for i in range(len(library)):
        file.write(">" + library_names[i] + "\n" + library[i] + "\n")
    file.close()
