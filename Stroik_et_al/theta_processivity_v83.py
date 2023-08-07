import pandas as pd

# WT_83 = sequencing initiation-->"GATCCTTTGCTGATTTACTACTTTGATACTTTATATGTTTGCGAGTTTCTATGTTTAGTATGCAG"<--MH/junction (synthesis proceeds right to left)

wt_xmers = ['GATCCTTTGCTGA','GCTGATTTACTAC','ACTACTTTGATAC','GATACTTTATATG','ATATGTTTGCGAG',            # establish unmutated sequences
         'GCGAGTTTCTATG','CTATGTTTAGTAT']
del1_xmers = ['GATCCTTGCTGA','GCTGATTACTAC','ACTACTTGATAC','GATACTTATATG','ATATGTTGCGAG',               # establish sequences with single nt deletion in TTT repeat
              'GCGAGTTCTATG','CTATGTTAGTAT']
ins1_xmers = ['GATCCTTTTGCTGA','GCTGATTTTACTAC','ACTACTTTTGATAC','GATACTTTTATATG','ATATGTTTTGCGAG',     # establish sequences with single nt insertion in TTT repeat
              'GCGAGTTTTCTATG','CTATGTTTTAGTAT']

TTT_location = [1,2,3,4,5,6,7]                                                                          # define TTT positions relative to sequencing reads, junctions is nearest to position 7
data_compiled = {'TTT_position': TTT_location, 'TTT_barcode': wt_xmers}                                 # establish dictionary for compiled data output at all TTT positions

# input CSV name with list of CSV read files to analyze
files0 = input ('Input CSV with list of all read files to analyze:')                                    # input adn read .csv with all file names to run and compile a full experiment at once
files1 = pd.read_csv(files0, names=['files'])
files2 = files1['files'].tolist()

for file in files2:
    df = pd.read_csv(file + '.csv', names=["ID", "blank", "reads"])
    raw_reads = df['reads'].tolist()                                                                    # pull seqeunce reads from each .csv, our reads are palced in column 3 based on our output workflow

    filter_reads_1 = []                                                                                 # establish empty lists for filtering reads
    filter_reads_2 = []

    for read in raw_reads:
        rev = read[::-1]                                                                                # reverse complement to read strand with TTT repeats
        revcom = rev.replace("A", "x").replace("C", "y").replace("G", "C").replace("T", "A").replace("x", "T").replace("y", "G")

        if revcom.find("AGTATGCAG") != -1:                                                              # filter for contaminating reads, all control and cellular synthesis samples of this substrate will contain this barcode
            filter_reads_1.append(revcom)

    if str(file).find("83") != -1:                                                                      # for files named as cellular synthesis ('83-1.csv' for example)
        for read in filter_reads_1:
            if read.find("AGTATGCAGATCT") != -1:                                                        # specific extended barcode for cellular synthesis samples, filtering out contaminating control reads
                filter_reads_2.append(read)
    if str(file).find("cntrl") != -1:                                                                   # for files named as cellular synthesis ('cntrl-1.csv' for example)
        for read in filter_reads_1:
            if read.find("AGTATGCAGTTTT") != -1:                                                        # specific extended barcode for control synthesis samples, filtering out contaminating cellular reads
                filter_reads_2.append(read)

    lst2 = [[] for x in range(0,7)]
    reads_2 = []

    for k in range (0,7):                                                                               # for each TTT position, starting at sequencing initiation, ending at repair junction
        wt = 0                                                                                          # establish 0 value for wt, 1 bp deletion, 1 bp insertion; will add to this as matches are identified
        del1 = 0
        ins1 = 0
        no_match = 0                                                                                    # estab;iosh 0 value for reads which don't match queries (wt_xmers, del1_xmers, ins1_xmers)

        for read in filter_reads_2:                                                                     # for each read in filtered list
            if read.find(wt_xmers[k]) != -1:                                                            # check for matches to wt xmers (containing TTT at each position)
                lst2[k].append("wt")                                                                    # add identity of read at each position to list
                wt += 1                                                                                 # if match found, add 1 to count of wt matches at that position

            if read.find(del1_xmers[k]) != -1:                                                          # check for matches to 1 bp deletion xmers (containing TT at each position)
                lst2[k].append("del1")                                                                  # add identity of read at each position to list
                del1 += 1                                                                               # if match found, add 1 to count of deletion matches at that position

            if read.find(ins1_xmers[k]) != -1:                                                          # check for matches to 1 bp insert xmers (containing TTTT at each position)
                lst2[k].append("ins1")                                                                  # add identity of read at each position to list
                ins1 += 1                                                                               # if match found, add 1 to count of insertion matches at that position


            kmers = [wt_xmers[k], del1_xmers[k], ins1_xmers[k]]                                         # compile identity of each read
            if any(x in read for x in kmers):                                                           # scan compiled read identities; if any are classified, continue
                continue
            else:                                                                                       # if identity is unspecified at all 7 positions...
                lst2[k].append("no match")                                                              # define read as 'no match'
                no_match += 1                                                                           # add 1 to count of no match reads

        indel = del1 + ins1                                                                             # total mutation occurrences, addition of 1 bp insertions and deletions at each position

        print("TTT location:", k)                                                                       # for checking progress of code, printing outputs at each TTT position for each library
        print("no match:", no_match)
        print("no mutation:", wt)
        print("del1:", del1)
        print("ins1:", ins1)
        print("indel:", indel)
        print("")

# TTT_counting
    wt_counts = []                                                                                      # open empty list for compiling wt vs indel counts and frequencies
    del1_counts = []
    ins1_counts = []
    no_match_counts = []
    indel_counts = []
    indel_frequency = []

    for k in range (0,7):                                                                               # for each TTT position
        wt_count = 0                                                                                    # establish 0 value for counts of wt, insertion, deletion, and failed matches
        del1_count = 0
        ins1_count = 0
        no_match_count = 0

        for unit in lst2[k]:                                                                            # scan list of mutation status at each position in each read, add to corresponding value
            if (unit == "wt"):
                wt_count += 1
            if (unit == "del1"):
                del1_count += 1
            if (unit == "ins1"):
                ins1_count += 1

            if (unit == "no match"):
                no_match_count += 1

        indel_count = del1_count + ins1_count                                                           # add insertion and deletion counts together, indel frequency is of interest for pol theta mutation signatures

        wt_counts.append(wt_count)                                                                      # append counts of all identities at each positions for each library
        del1_counts.append(del1_count)
        ins1_counts.append(ins1_count)
        indel_counts.append(indel_count)
        no_match_counts.append(no_match_count)

        if indel_count + wt_count != 0:                                                                 # if indel count is non-zero
            indel_freq = indel_count / (indel_count + wt_count)                                         # calculate indel frequency, using total matched reads as denominator
            indel_frequency.append(indel_freq)                                                          # append indel frequency to list at each position, organized for output
        else:                                                                                           # if indel count is zero, it will not already be established as such
            indel_freq = 0                                                                              # force establish of indel count = 0 in this case
            indel_frequency.append(indel_freq)                                                          # append 0 value to indel count in this case
        k += 1                                                                                          # iterate to next TTT position

    data2 = {'TTT_position': TTT_location, 'TTT_barcode': wt_xmers,                                     # establish disctionary for DataFrame output to .csv
             'wt': wt_counts, 'del1': del1_counts, 'ins1':ins1_counts, 'no match': no_match_counts,     # including all identity counts and total indel frequency for each individual sample input
             'indels': indel_counts, 'indel frequency': indel_frequency}

    analysis_df = pd.DataFrame(data2)                                                                   # convert dict to DataFrame
    analysis_df_transposed = analysis_df.T                                                              # transpose DataFrame so TTT positions read left to right, instead of top to bottom
    analysis_df_transposed.to_csv(file + '_indel_counting.csv')                                         # export transposed DataFrame to .csv for each input file with matching name

    data_compiled[str(file)] = indel_frequency                                                          # establish indel frequency at each TTT position for unique file names, add as a new column to 'data_compiled' dictionary established line 13

    print("fin", file)                                                                                  # signal completion of each file and transition to the next file as listed in files2 (see lines 16-18)

compiled_df = pd.DataFrame(data_compiled)                                                               # convert 'data_compiled' dictionary to DataFrame
compiled_df.to_csv('compiled_indel_frequency.csv')                                                      # output compiled DataFrame to .csv

print("fin fin")                                                                                        # signal completion of all loops for all files without error