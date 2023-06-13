import pandas as pd

# define code blocks
# for substitution miscall filter
def miscall():                                                                  # if substitution identified in junction flanks, append junction analysis info into miscall lists
    jxn = read[lpos:rpos + 10 + MH_len]
    miscall_jxns.append(jxn)
    sequences_miscall.append(read)
    left10_miscall.append(left_xmers[ldel])
    right10_miscall.append(right_xmers[rdel])
    left_deletion_miscall.append(ldel)
    right_deletion_miscall.append(rdel)
    MH_size_miscall.append(MH_len)
    MH_ID_miscall.append(MH_nt)
    total_deletion_miscall.append(ldel + rdel + MH_len)
    insert_size_miscall.append(insert_len)
    insert_ID_miscall.append(insert_nt)

def no_miscall():                                                               # if no substitution identified in junction flanks, append junction analysis info into data output lists
    sequences.append(read)
    junctions.append(jxn)
    left_deletion.append(ldel)
    left10.append(left_xmers[ldel])
    right_deletion.append(rdel)
    right10.append(right_xmers[rdel])
    MH_size.append(MH_len)
    MH_ID.append(MH_nt)
    insert_size.append(insert_len)
    insert_ID.append(insert_nt)
    total_deletion.append(tdel)

# define locus - there are different approaches for defining left and right flanks below, remove # signs to activate the code block you wish to use
# R26-MHD (R26A) and R26-TINS (R26E) are predefined below because we use them frequently
# manual input can be formatted with a 5-NGG-3 PAM, a 5-CCN-3 PAM, or with manual input of the left and right flanks with respect to the breaksite

# R26-MHD (R26A); note orientation of R26-MHD, sequence reads in csv input files are the strand containing the 5-NGG-3 PAM
# In this orientation, 5-NGG-3 PAM is downstream (3') relative to the gRNA target sequence
#locus = 'AGTTCTCTGCTGCCTCCTGGCTTCTGAGGACCGCCCTGGGCCTGGGAGAATCCCTTCCCCCTCTTCCCTCGTGATCTGCAACTCCAGTCTTTCTAGAAGATGGGCGGGAGTCTTCTGGGCAGGCTTAAAGGCTAACCTGGTGTGTGGGCGTTGTCCTGCAGGGGAATTGAACAGGTGTAAAATTGGAGGGACAAGACTTCCCACAGA'
#gRNA = 'ACTCCAGTCTTTCTAGAAGA'
#x = locus.find (sgRNA)                            # identify 5' location of gRNA annealing site, cut site is 17 bp from 5' location (5' location + 17)
#left_direct = locus[0:x + 17]                      # define left flank as the sequence upstream (left) of the cut site
#right_direct = locus[x + 17:len(locus)]            # define right flank as the sequence downstream (right) of the cut site

#R26-TINS (R26E); note orientation of R26-TINS, sequence reads in csv input files are the strand containing the 5-CCN-3 PAM
#locus = 'TGGCTTATCCAACCCCTAGACAGAGCATTGGCATTTTCCCTTTCCTGATCTTAGAAGTCTGATGACTCATGAAACCAGACAGATTAGTTACATACACCACAAATCGAGGCTGTAGCTGGGGCCTCAACACTGCAGTTCTTTTATAACTCCTTAGTACACTTTTTGTTGATCCTTTGCCTTGATCCTTAATTTTCAGTGTCTATCACCTCTCCCGTCAGGTGGTGTTCCACA'
#sgRNA = 'CCACAAATCGAGGCTGTAGCTGG'
#y = locus.find (sgRNA)
#left_direct = locus[0:y + 6]
#right_direct = locus[y + 6:len(locus)]

#R26D; note orientation of R26-TINS, sequence reads in csv input files are the strand containing the 5-CCN-3 PAM
#locus = 'TATTGCTTATCTTTAGTTCCCAACACTTGGGAGGCAGAGGCCAGCCAGGGCTATGTGACAAAAACCTTGTCTAGAGGAGAAACTTCATAGCTTATTTCCTATTCACGTAACCAGGTTAGCAAAATTTACCAGCCAGAGATGAAGCTAACAGTGTCCACTATATTTGTAGTGTTTTAAGTCAATTTTTTAAATATACTTAATAGAATTAAAGCTATGGTGAACCAAGTACAA'
#sgRNA = 'CCTATTCACGTAACCAGGTTAGC'
#y = locus.find (sgRNA)
#left_direct = locus[0:y + 6]
#right_direct = locus[y + 6:len(locus)]

# LBR (MEF204); note orientation of R26-MHD, sequence reads in csv input files are the strand containing the 5-NGG-3 PAM
# In this orientation, 5-NGG-3 PAM is downstream (3') relative to the gRNA target sequence
#locus = 'ATGGGTTTGGGCTTGGAGTAGTGGGGAGGAGCCCACAAGTTTTCGGGATGAAGTCAGGCCTGGCTGAGTCACCCGAGGGGCACTTAGGGAGCAGGGGGGCAGTGAACACCTCTGCATGAGCAGGGGCATAAAAACGGAAGGTGACAACATCCAGCTTTGCCTATTTCAACATTTAGCTCAGAGCCTCCAAGTACAAAGAAAGAGGAAGGAAATGTACCCTCCTTCTCTCTTCTCTTG'
#sgRNA = 'GAACACCTCTGCATGAGCAG'
#x = locus.find (sgRNA)                            # identify 5' location of gRNA annealing site, cut site is 17 bp from 5' location (5' location + 17)
#left_direct = locus[0:x + 17]                      # define left flank as the sequence upstream (left) of the cut site
#right_direct = locus[x + 17:len(locus)]            # define right flank as the sequence downstream (right) of the cut site

# HPRT exon7 gRNA from Tjisterman, et al (Cell Reports), sequence reads in csv input files are the strand containing the 5-NGG-3 PAM
# In this orientation, 5-NGG-3 PAM is downstream (3') relative to the gRNA target sequence
#locus = 'AATTTATTACACTCCAAATTCATAATCAAGTTCTAATTTTTTTAGAATTCTAATTTTTTATCTAAAAAAATTCAACCTTTTTGATTCCACAGGGACAAACACCGTTGTTTCCACGTATCTTCGGCCATGAAGCTGGAGGGTAATAGAAACACTAATCTTCTTTGCTTCGTTTTGGATATTTTTAAGGTTTTAGAGATTCAAGGTCGTTTTTTTTGTTGTTGTGTAGGATTGTTGAGAGTGTTGGAGAA'
#sgRNA = 'ATCTTCGGCCATGAAGCTGG'
#x = locus.find (sgRNA)                            # identify 5' location of gRNA annealing site, cut site is 17 bp from 5' location (5' location + 17)
#left_direct = locus[0:x + 17]                      # define left flank as the sequence upstream (left) of the cut site
#right_direct = locus[x + 17:len(locus)]            # define right flank as the sequence downstream (right) of the cut site

# Manual input of locus and gRNA; note the orientation of your locus of interest, use this block if the sequence reads in the input file contain 5'-NGG-3' PAM
locus = input ('Input locus sequence containing 5-NGG-3 PAM:')
sgRNA = input ('Input gRNA sequence containing 5-NGG-3 PAM:')
x = locus.find (sgRNA)
left_direct = locus[0:x + 17]
right_direct = locus[x + 17:len(locus)]

# Manual input of locus and gRNA; note the orientation of your locus of interest, use this block if the sequence reads in the input file contain 5'-CCN-3' PAM
#locus = input ('Input locus sequence containing 5-CCN-3 PAM:')
#sgRNA = input ('Input gRNA sequence containing 5-CCN-3 PAM:')
#y = locus.find (sgRNA)
#left_direct = locus[0:y + 6]
#right_direct = locus[y + 6:len(locus)]

# Manual input of left and right flank germline sequences
#left_direct = input ('Input left flank sequence 5-3prime (breaksite at 3prime end):')
#right_direct = input ('Input right flank sequence 5-3 prime (breaksite at 5prime end):')
#locus = left_direct + right_direct

# defining list of left and right 1 nt walking 10mers for junction comparison
left_xmers = {}                                                                 # open empty dictionary for list of left_xmers
for i in range(0, len(left_direct)-10):                                          # iterate loop for possible walking 1 nt deletions, stop loop when you reach last possible 10mer
    left_xmers.update ({i : left_direct[len(left_direct)-10-i:len(left_direct)-i]})    # append 10mers of left flank corresponding to iterative deletions

right_xmers = {}                                                                # repeat same process as above for the right flank and generate list of right flank walking  10mers
for i in range(0, len(right_direct)-10):
    (right_xmers.update ({i : right_direct[i:i+10]}))

germline_jxn = left_direct[-10:] + right_direct[:10]                              # define germline junction of 10 nt left and 10 nt right, for use in frequency calculations later

# to analyze multiple files in a single experiment (single locus), input .csv file name with list of .csv read files to analyze
files0 = input ('Input CSV with list of all read files to analyze:')            # provide name of .csv file with list of files to analyze
files1 = pd.read_csv(files0, names = ['files'])
files2 = files1['files'].tolist()                                               # converts provided list of .csv names in iterable list of files

library = []
initial_reads = []
jxns_w_ambig = []
failed_matches = []
jxns_w_subs = []
filtered_jxns = []

# start of junction analysis module
for file in files2:                                                             # iterate junction analysis loop for each file in .csv list
    print('start', file)
    df = pd.read_csv(file + '.csv', names = ["ID", "blank", "reads"])               # interpret .csv file with raw reads in column 3, this is where our .csv output lists reads, generated in CLC Genomics Workbench
    raw_reads = df['reads'].tolist()                                                # place raw reads from .csv file into list of reads

    # open lists of all outputs for csv
    # for junction analysis (reads with accurate left/right 10mer matches and no base substitutions will be placed here)
    sequences = []
    left10 = []
    right10 = []
    junctions = []
    left_deletion = []
    right_deletion = []
    MH_size = []
    MH_ID = []
    total_deletion = []
    insert_size = []
    insert_ID = []
    ambiguity = []
    failed = 0                                                                      # open empty count of reads which fail to map either a left or right flank 10mer match

    # for ambiguity filter (reads from above lists will pass through ambiguity filter and be placed here if no ambiguities occur in the defined junction)
    final_jxns = []
    sequences_2 = []
    left10_2 = []
    right10_2 = []
    junctions_2 = []
    left_deletion_2 = []
    right_deletion_2 = []
    MH_size_2 = []
    MH_ID_2 = []
    total_deletion_2 = []
    insert_size_2 = []
    insert_ID_2 = []

    # for substitution miscall filter (reads from miscall/substitution filter that display evidence of a base substitution in the junction flanks (10 nt to either side) will be placed here)
    substitutions = 0                                                               # open empty count of junctions with evident base substitution
    miscall_jxns = []
    sequences_miscall = []
    left10_miscall = []
    right10_miscall = []
    junctions_miscall = []
    left_deletion_miscall = []
    right_deletion_miscall = []
    MH_size_miscall = []
    MH_ID_miscall = []
    total_deletion_miscall = []
    insert_size_miscall = []
    insert_ID_miscall = []

    # define mutational signatures
    for read in raw_reads:                                                          # iterate over junction analysis loop for each read in the raw reads list
        # Deletion identifier
        # define left and right deletion sizes by aligning pre-defined xmers
        ldel = 0                                                                        # define left and right deletion initially as 0 for each read, then walk 1 nt at a time until a match occurs int eh read in question
        rdel = 0
        MH_len = 0                                                                      # define microhomology length as 0 initially, will be redefined if MH is found in the junction
        while read.find(left_xmers[ldel]) == -1:                                        # starting at a deletion size of 0, if the left 10mer match fails, iterate left deletion size in 1 nt increments until match is found
            if ldel == len(left_xmers) - 1:                                                 # breaks loop if you run out of 10mers
                lpos = -1                                                                   # define left match position as -1, indicating a failed match
                failed += 1                                                                 # count if left match fails (if you run out of 10mers to test), count reflects reads that are removed from final analysis due to failed flank match
                break
            else:                                                                       # if the current attempted match fails
                ldel += 1                                                                   # increase left deletion size by one and check next left 10mer
        else:                                                                           # if the current left 10mer match DOES NOT fail (i.e. if it is found somewhere in the read)
            lpos = int(read.find(left_xmers[ldel]))                                         # define left position as the match location and left deletion size as the current iteration count (ldel), terminate left search loop
        while read.find(right_xmers[rdel]) == -1:                                       # same loop as above, now identifying the location of right 10mer match and right deletion size
            if rdel == len(right_xmers) - 1:
                rpos = -1
                failed += 1
                break
            else:
                rdel += 1

        else:
            rpos = int(read.find(right_xmers[rdel]))

        jxn = read[lpos:rpos+10]                                                        # define junction as the sequence between the left match position and the right match position + 10 (accounting for the 10mer match length)
                                                                                        # junction will be 20 nt if simple deletion, 10 from each flank), >20 if an insertion is present, and <20 if microhomology is present due to the inherent MH-dependent overlap in left/right junctions

        # Microhomology identifier
        # jxn defined as 20mer (left10 + right10); less than 20 means left/right xmer overlap, i.e. a microhomology
        if lpos != -1 and rpos != -1:                                                   # only consider reads where both left and right sides of junction found a match, ignore if either flank match failed (-1)
            if (len(jxn)) < 20:                                                             # MH present if jxn length is <20
                MH_len = 20 - len(jxn)                                                          # expecting 20 nt jxn, MH length corresponds to the difference from 20
                MH_nt = right_xmers[rdel][0:MH_len]                                             # define MH sequence as the first x nts of the right flank, where x is the length of MH
                MH_extend = read[rpos+10:rpos+10+MH_len]                                        # extend the length of MH-containing jxns so that all junctions are uncer the same influence of amplification-induced base substitution (see substitution filter below)
                locus_position = locus.find(right_xmers[rdel]) + 10                             # define location of MH extension to compare each MH read to the germline (expected) extension sequence

                # when MH forces junction extension, filter out reads with base substitution in right flank extension
                if MH_extend == locus[locus_position:locus_position+MH_len]:                    # if MH extension sequence matches corresponding germline sequence, read is valid
                    jxn = read[lpos:lpos + 20]                                                      # redefine MH-containing jxn to include extended nts, 20 nt in length
                    check_MH_extend = 1                                                             # jxn passes MH substitution filter (defined as ==1)
                else:                                                                               # if MH extension DOES NOT match corresponding germline sequence, there is a substition adjacent to the right 10mer match
                    check_MH_extend = 0                                                         # jxn fails to pass MH substitution filter, has a substitution (defined as ==0)
                    miscall()                                                                       # append jxn anaylsis info to miscall lists, ignore for final data output file
                    substitutions += 1                                                              # count substituted read, count reflects reads that are removed from final analysis due to flanking base substitution

            else:                                                                           # if jxns are >=20, MH is not present and read will undergo substitution filter in next code block
                check_MH_extend = 1                                                             # jxn passes MH substitution filter (defined as ==1)
                MH_len = 0                                                                      # if jxn length >=20, no MH present
                MH_nt = ''

            tdel = ldel + rdel + MH_len                                                     # define total deletion size as the sum of left deletion, right deletion, and MH-mediated deletion

            # Insertion identifier
            # jxn defined as 20mer (left10 + right10); greater than 20 is insertion
            if (len(jxn)) > 20:                                                             # if junction >20, insertion present between left and right 10mer matches
                insert_len = len(jxn) - 20                                                      # expecting 20 nt jxn, insertion length corresponds to the difference from 20
                insert_nt = jxn[10:10 + insert_len]                                             # define inserted nts as the nucleotides between left and right 10mer matches
            else:                                                                           # if junction <= 20, no insertion present
                insert_len = 0                                                                  # no insertion, thus insertion length = 0
                insert_nt = ''                                                                  # define insertion identity as blank

            # filter out reads with substitution (base miscall) in 10 flanking nts to left and right
            # substitution in left or right flank will push the left or right 10mer match further into the flank until it is past the substitution
            # this issue will incorrectly identify an insertion between the left and right flanks with length corresponding to the location of the substitution
            if len(insert_nt) < 3:                                                          # limit substitution searching to within 3 nt of the flank end, if pseudo-insertion is <3, we cannot properly determine if a substitution is confounding interpretation of the read
                check_left = 1                                                                  # jxn passes left and right substitution filters (defined as ==1)
                check_right = 1

            if len(insert_nt) >= 3 and ldel >= 3:                                           # if insertion is >= 3 and left deletion is >= 3, pseudo insertion is possible
                left_kmer = left_xmers[ldel-3]                                                  # walk the left 10mer 3 positionscloser to the junction to include the potential substitution base and two nts upstream (5') of this position, will be used for comparison to insert
                l1 = left_kmer[-2:]                                                             # isolate two nts upstream of the potential substitution

                if insert_nt[1:3] != l1:                                                        # if 2nd and 3rd nts of insert DO NOT match the 2 nts isolated from the adjusted left 10mer (skipping the potentially substituted base), indel call is accurate and no substitution is present
                    check_left = 1                                                                  # jxn passes left substitution filter (defined as ==1)
                else:                                                                           # if 2nd and 3rd nts of the insert match the 2 nts isolated from the adjusted left 10mer (skipping the potentially substituted base), substitution is a valid explanation for generation of a pseudo-insertion
                    miscall()                                                                       # append jxn anaylsis info to miscall lists, ignore for final data output file
                    check_left = 0                                                                  # jxn fails to pass left substitution filter, has a substitution (defined as ==0)
                    substitutions += 1                                                              # count substituted read, count reflects reads that are removed from final analysis due to flanking base substitution
                    continue
            else:                                                                           # if insert >= 3, but left deletion < 3, substitution is not possible
                check_left = 1                                                                  # jxn passes left substitution filter (defined as ==1)

            if len(insert_nt) >= 3 and rdel >= 3:                                           # same analysis as above, adjusted for the right flank
                right_kmer = right_xmers[rdel-3]
                r1 = right_kmer[:2]

                if insert_nt[-3:-1] != r1:
                    check_right = 1
                else:
                    miscall()
                    check_right = 0
                    substitutions += 1
                    continue
            else:
                check_right = 1

            if check_left == 1 and check_right == 1 and check_MH_extend == 1:
                no_miscall()                                                                    # if jxn passes MH, left, and right substitution filters, append anaylsis info to data output file to be included in final analysis

    # filter out reads with ambiguities in defined junction (e.g. N, R, Y, K, S, W, etc)
    ambiguity_jxns = 0                                                                  # open ambiguity jxn count at 0
    for i in range(0,len(junctions)-1):                                                 # iterate over the list of jxns
        s1 = junctions[i]                                                                   # identify jxn analyzed for this loop
        ambig = s1.replace("A","").replace("C","").replace("G","").replace("T","")          # delete all ACGT bases, leave ambiguity calls behind

        if len(ambig) == 0:                                                                 # if no bases remain after removal of ACGT, write unambiguous sequence and corresponding analysis info to new lists for .csv output
            final_jxns.append(s1)
            sequences_2.append(sequences[i])
            left10_2.append(left10[i])
            right10_2.append(right10[i])
            left_deletion_2.append(left_deletion[i])
            right_deletion_2.append(right_deletion[i])
            MH_size_2.append(MH_size[i])
            MH_ID_2.append(MH_ID[i])
            total_deletion_2.append(total_deletion[i])
            insert_size_2.append(insert_size[i])
            insert_ID_2.append(insert_ID[i])

        else:                                                                           # if any base calls reamin after removal of ACGT, they are ambiguities
            miscall()
            ambiguity_jxns += 1                                                             # count reflects reads that are removed from final analysisdue to one or more ambiguous base calls
            continue

    # start of TINS (templated insertions) module
    TINS = []                                                                       # open empty lists for accumulation of TINS data
    TINS_type = []
    TINS_initiation = []

    left_rev = left_direct[::-1]                                                    # define left reverse flank (reading sequence from right to left)
    left_revcom = left_rev.replace("A", "x").replace("C", "y").replace("G", "C").replace("T", "A").replace("x", "T").replace("y", "G")              # define left reverse complement, replace A to T, C to G, G to C, and T to A in reverse sequence from last step
    right_rev = right_direct[::-1]                                                  # same as above, adjusted for right flank
    right_revcom = right_rev.replace("A", "x").replace("C", "y").replace("G", "C").replace("T", "A").replace("x", "T").replace("y", "G")

    x = -1                                                                          # establish x = -1 to initiate loop at list item 0

    for insert in insert_ID_2:                                                          # iterate through list of inserts
        x += 1                                                                          # initiate first iteration of loop at index 0 in list of inserts
        right_jmer = right10_2[x]                                                   # establish right flank as the list item corresponding to the current insert in question
        left_jmer = left10_2[x]                                                     # establish left flank as the list item corresponding to the current insert in question

        if len(insert) < 3:                                                             # only consider insert >= 3 for potential TINS, any insert < 3 is skipped over
            TINS.append(' ')
            TINS_type.append(' ')
            TINS_initiation.append(' ')

        elif len(insert) >= 5:                                                          # if insert length is >= 5, scan for template
            s1 = insert[:5]                                                                 # define s1 as first 5 nt (5' end of insert) for comparison to flanks
            s2 = insert[-5:]                                                                # define s2 as last 5 nt (3' end of insert)
            s3 = s2[::-1]                                                                   # define reverse of s2 (reading from right to left)
            s4 = s3.replace("A", "x").replace("C", "y").replace("G", "C").replace("T", "A").replace("x", "T").replace("y", "G")                     # define complement of s3 (s4 is reverse complement of s1)

            # all TINS are limited to initiation within 30 nts of the junction
            if 0 <= right_direct.find(s1) <= 30:                                            # top strand insert, direct TINS (compare s1 against right direct)
                TINS.append(s1)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(right_direct.find(s1))
                TINS_type.append('DR')
            elif 0 <= left_revcom.find(s1) <= 30:                                           # top strand insert, inverted TINS (compare s1 against left revcom)
                TINS.append(s1)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(left_revcom.find(s1))
                if left_deletion_2[x] >= TINS_initiation[x]:
                    TINS_type.append('SS')                                                       # define TINS class and append to list for later data output, applicable for all TINS loops
                else:
                    TINS_type.append('IR')
            elif 0 <= left_revcom.find(s4) <= 30:                                           # bottom strand insert, direct TINS (compare s4 against left revcom)
                TINS.append(s4)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(left_revcom.find(s4))
                TINS_type.append('DR')
            elif 0 <= right_direct.find(s4) <= 30:                                          # bottom strand insert, inverted TINS (compare s4 against right direct)
                TINS.append(s4)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(right_direct.find(s4))
                if right_deletion_2[x] >= TINS_initiation[x]:
                    TINS_type.append('SS')                                                       # define TINS class and append to list for later data output, applicable for all TINS loops
                else:
                    TINS_type.append('IR')
            else:                                                                           # if none of the above match, insert is non-templated
                TINS_type.append('not templated')
                TINS.append(' ')
                TINS_initiation.append(' ')

        elif len(insert) == 4:                                                          # if insert length is 4, infer 1 nt to define 5mer for TINS search
            s2l = left_jmer[-1] + insert                                                   # infer 1 nt in left flank (5' of insertion) for search using "bottom strand" of insert
            s3l = s2l[::-1]                                                                 # define reverse complement of s2l
            s4l = s3l.replace("A", "x").replace("C", "y").replace("G", "C").replace("T", "A").replace("x", "T").replace("y", "G")

            s2r = insert + right_jmer[0]                                                   # infer 1 nt in right flank (3' of insertion) for search using "top strand" of insert

            if 0 <= right_direct.find(s2r) <= 30:                                            # top strand insert, direct TINS (compare s1 against right direct)
                TINS.append(s2r)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(right_direct.find(s2r))
                TINS_type.append('DR')
            elif 0 <= left_revcom.find(s2r) <= 30:                                           # top strand insert, inverted TINS (compare s1 against left revcom)
                TINS.append(s2r)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(left_revcom.find(s2r))
                if left_deletion_2[x] >= TINS_initiation[x]:
                    TINS_type.append('SS')                                                       # define TINS class and append to list for later data output, applicable for all TINS loops
                else:
                    TINS_type.append('IR')
            elif 0 <= left_revcom.find(s4l) <= 30:                                           # bottom strand insert, direct TINS (compare s4 against left revcom)
                TINS.append(s4l)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(left_revcom.find(s4l))
                TINS_type.append('DR')
            elif 0 <= right_direct.find(s4l) <= 30:                                          # bottom strand insert, inverted TINS (compare s4 against right direct)
                TINS.append(s4l)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(right_direct.find(s4l))
                if right_deletion_2[x] >= TINS_initiation[x]:
                    TINS_type.append('SS')                                                       # define TINS class and append to list for later data output, applicable for all TINS loops
                else:
                    TINS_type.append('IR')
            else:                                                                           # if none of the above match, insert is non-templated
                TINS_type.append('not templated')
                TINS.append(' ')
                TINS_initiation.append(' ')

        elif len(insert) == 3:                                                          # if insert length is 3, infer 2 nts to define 5mer for TINS search
            s2l = left_jmer[-2:] + insert                                                  # infer 2 nts in left flank (5' of insertion) for search using "bottom strand" of insert
            s3l = s2l[::-1]                                                                 # define reverse complement of s2l
            s4l = s3l.replace("A", "x").replace("C", "y").replace("G", "C").replace("T", "A").replace("x", "T").replace("y", "G")

            s2r = insert + right_jmer[:2]                                                  # infer 2 nts in right flank (3' of insertion) for search using "top strand" of insert

            if 0 <= right_direct.find(s2r) <= 30:                                            # top strand insert, direct TINS (compare s1 against right direct)
                TINS.append(s2r)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(right_direct.find(s2r))
                TINS_type.append('DR')
            elif 0 <= left_revcom.find(s2r) <= 30:                                           # top strand insert, inverted TINS (compare s1 against left revcom)
                TINS.append(s2r)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(left_revcom.find(s2r))
                if left_deletion_2[x] >= TINS_initiation[x]:
                    TINS_type.append('SS')                                                       # define TINS class and append to list for later data output, applicable for all TINS loops
                else:
                    TINS_type.append('IR')
            elif 0 <= left_revcom.find(s4l) <= 30:                                           # bottom strand insert, direct TINS (compare s4 against left revcom)
                TINS.append(s4l)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(left_revcom.find(s4l))
                TINS_type.append('DR')
            elif 0 <= right_direct.find(s4l) <= 30:                                          # bottom strand insert, inverted TINS (compare s4 against right direct)
                TINS.append(s4l)                                                                 # append TINS match as templated nts, applicable for all TINS loops
                TINS_initiation.append(right_direct.find(s4l))
                if right_deletion_2[x] >= TINS_initiation[x]:
                    TINS_type.append('SS')                                                       # define TINS class and append to list for later data output, applicable for all TINS loops
                else:
                    TINS_type.append('IR')
            else:                                                                           # if none of the above match, insert is non-templated
                TINS_type.append('not templated')
                TINS.append(' ')
                TINS_initiation.append(' ')

        else:
            x += 1                                                                          # increase index on insert list by 1
            continue

    data = {'reads': sequences_2,                                                   # compile lists of jxn analysis data from jxns passing all filters, organize into data frame
            'left flank': left10_2, 'right flank': right10_2, 'reconstructed junction': final_jxns,
            'left deletion': left_deletion_2, 'right deletion': right_deletion_2,
            'MH size': MH_size_2, 'MH ID': MH_ID_2, 'total deletion': total_deletion_2,
            'insert size': insert_size_2, 'insert ID': insert_ID_2,
            'TINS sequence': TINS, 'TINS initiation': TINS_initiation,'TINS class': TINS_type}

    miscall_data = {'reads': sequences_miscall,                                     # compile lists of jxn analysis data from jxns passing all filters, organize into data frame
                    'left flank': left10_miscall, 'right flank': right10_miscall, 'reconstructed junction': miscall_jxns,
                    'left deletion': left_deletion_miscall, 'right deletion': right_deletion_miscall,
                    'MH size': MH_size_miscall, 'MH ID': MH_ID_miscall, 'total deletion': total_deletion_miscall,
                    'insert size': insert_size_miscall, 'insert ID': insert_ID_miscall}


    analysis_df = pd.DataFrame(data)
    analysis_df.to_csv(file + '_jxn_analysis.csv')                                  # output data frame of analysis for jxns passing filters into .csv

#    miscall_analysis_df = pd.DataFrame(miscall_data)                                # output data frame of analysis for jxns failing filters into .csv
#    miscall_analysis_df.to_csv(file + '_miscall_results.csv')

    print('reads:', len(raw_reads))                                                 # print initial number of raw reads
    print('jxns with ambiguity:', ambiguity_jxns)                                   # print number of jxns filtered out due to ambiguity
    print('failed matches:', failed)                                                # print number of jxns filtered out due to failed flank match
    print('jxns with substitution:', substitutions)                                 # print number of jxns filtered out due to substitution (MH extension, left flank, or right flank)
    print('final junctions:', len(final_jxns))                                      # print final number of jxns to be used in analysis, passing all above filters

    library.append(file)
    initial_reads.append(len(raw_reads))
    jxns_w_ambig.append(ambiguity_jxns)
    failed_matches.append(failed)
    jxns_w_subs.append(substitutions)
    filtered_jxns.append(len(final_jxns))

    print("fin jxn analysis -", file, '\n')                                       # identifies completion of the TINS analysis code block for each file

# start of frequency module
frames = []                                                                     # open empty list for compilation of jxn analysis from all analyzed files in list from "files2"

# establish list of unique jxns and corresponding jxn analysis data
for file in files2:                                                             # iterate through all files in files2 list
    file_data = pd.read_csv(file + '_jxn_analysis.csv', usecols = ['reads', 'left flank', 'right flank',               # read jxn analysis .csv file from each file in list, extracting columns listed
                                                                    'reconstructed junction', 'left deletion', 'right deletion',
                                                                    'MH size', 'MH ID', 'total deletion', 'insert size', 'insert ID',
                                                                    'TINS sequence', 'TINS initiation', 'TINS class'],
                                                                    keep_default_na = False)
    frames.append(file_data)                                                        # append all data extracted from jxn analysis files into single data frame

compiled_df = pd.concat(frames)                                                 # concatenate data frames from all jxn analysis files
print('compiled', len(compiled_df))                                             # print number of compiled reads from all libraries

compiled_df.drop_duplicates(subset ="reconstructed junction", keep = 'first', inplace = True)                           # removes duplicates from unique_df
unique_df = compiled_df.reset_index(drop = True)                                # reset index after removing duplicates
print('unique', len(unique_df))                                                 # print number of unique reads from all libraries (duplicates removed)

unique_jxns = unique_df['reconstructed junction'].tolist()                      # establish list of unique reconstructed jxns for frequency counting

for file in files2:                                                             # iterate through all files in files2 list
    data_df = pd.read_csv(file + '_jxn_analysis.csv')                              # read data from file in question
    jxns = data_df['reconstructed junction'].tolist()                               # extract list of reconstructed jxns from file in question

    jxn_count = []                                                                  # open empty lists for frequency data accumulation
    raw_jxn_freq = []
    repair_jxn_freq = []

    for i in range(0,len(unique_jxns)):                                             # iterate over each jxn from compiled unique list
        if germline_jxn.find(unique_jxns[i]) != -1:                                     # if jxn matches germline definition, omit from repair adjusted frequency
            count = jxns.count(unique_jxns[i])                                              # counting unique jxn matches in full jxn list
            frequency = count / len(jxns)                                                   # calculate the raw frequency of each unique jxn
            jxn_count.append(count)                                                         # append jxn count to list for output
            raw_jxn_freq.append(frequency)                                                  # append raw frequency to list for output
            repair_jxn_freq.append(' ')                                                     # append blank cell for repair adjusted frequency (because read is defined as germline
        else:                                                                           # for reads the DO NOT match germline jxn definition
            count = jxns.count(unique_jxns[i])                                              # counting unique jxn matches in full jxn list
            frequency = count / len(jxns)                                                   # calculate frequency of each unique jxn
            repair_frequency = frequency / (1-jxns.count(germline_jxn) / len(jxns))         # adjust jxn frequency for repair frequency (denominator adjusted to omit germline frequency)
            jxn_count.append(count)                                                         # append jxn count to list for output
            raw_jxn_freq.append(frequency)                                                  # append raw frequency to list for output
            repair_jxn_freq.append(repair_frequency)                                        # append repair adjusted frequency to list for output

    col1 = file + 'counts'                                                          # define new column as jxn count for each library
    col2 = file + 'raw freq'                                                        # define new column as raw frequency for each library
    col3 = file + 'repair freq'                                                     # define new column as repair adjusted frequency for each library
    compiled_df[col1] = jxn_count                                                   # append jxn count list to output data frame
    compiled_df[col2] = raw_jxn_freq                                                # append raw frequency list to output data frame
    compiled_df[col3] = repair_jxn_freq                                             # append repair adjusted frequency to output data frame

    print('fin frequency calc -', file)                                             # identifies completion of jxn frequency code block for each file

summary_data = {"library": library, "initial reads":initial_reads, "jxns with ambiguity":jxns_w_ambig,
                "failed matches":failed_matches, "jxns with substitution":jxns_w_subs, "filtered jxns":filtered_jxns}
summary_df = pd.DataFrame(summary_data)

summary_df.to_csv(str(files0[:-4]) + '_summary.csv')
compiled_df.to_csv(str(files0[:-4]) + '_full_compiled_freq.csv')                # generates final compiled output .csv containing jxn analysis, TINS classifications, and jxn frequencies for every library
print('fin fin')                                                                # identifies completion of full code process