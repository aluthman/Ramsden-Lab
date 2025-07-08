#sbatch --nodes=1 --ntasks=1 --cpus-per-task=9 --time=0-12:00:00 --mem=80g --wrap="python3 /nas/longleaf/home/aluthman/bubble_multi_v5.py"
#sacct -o JobID,JobName,state,CPUtime,maxRSS

import pandas as pd
import time
import multiprocessing as mp

def process_read(partition, local_lst):
    local_tags = []
    failed = 0
    # Define GAG substrate flanks
#    left_flank = 'AACTTCGTGAGGACACACTGTAAGCGATGCTCTCA?????TTCTGTTAGACTTGTGGTGGATGACAGAAGCAGGGTAGCCAGTCTGAGAGCTGAGGTC'
#    right_flank = 'TCCACCTTAGCTGTATAGTCACCCTGCAGATCTTCACTCTCACACCCATCAATAGCACGATTCACTCTGTTCCATGCAGCGGACTGTCCAAGTTTCAGATTGCGTA??????????CTCGACAACGTCACCGAAG'

    # Define GCGT substrate flanks
    left_flank = 'AACTTCGTGAGGACACACTGTAAGCGATGCTCTCA?????TTCTGTTAGACTTGTGGTGGATGACAGAAGCAGGGTAGCCAGTCTGAGAGCTGAGGCAT'
    right_flank = 'TAGCTTAGCTGTATAGTCACCCTGCAGATCTTCACTCTCACACCCATCAATAGCACGATTCACTCTGTTCCATGCAGCGGACTGTCCAAGTTTCAGATTGCGTA??????????CTCGACAACGTCACCGAAG'

    locus = left_flank + right_flank
    germline_jxn = left_flank[-10:] + right_flank[:10]

    # defining list of left and right 1 nt walking 10mers for junction comparison
    left_xmers = {}                                                                     # open empty dictionary for list of left_xmers
    for i in range(0, 20):                                                              # iterate loop for possible walking 1 nt deletions, stop loop when you reach last possible 10mer
        left_xmers.update ({i : left_flank[len(left_flank)-10-i:len(left_flank)-i]})    # append 10mers of left flank corresponding to iterative deletions

    right_xmers = {}                                                                    # repeat same process as above for the right flank and generate list of right flank walking  10mers
    for i in range(0, 20):
        right_xmers.update ({i : right_flank[i:i+10]})

    top_IDs = ['CCGAA','ACGAA','GCGAA','TCGAA','CAGAA','CGGAA','CTGAA','CCAAA','CCCAA','CCTAA','CCGCA','CCGGA','CCGTA','CCGAC','CCGAG','CCGAT']
    bottom_IDs = ['GGCTT','AGCTT','CGCTT','TGCTT','GACTT','GCCTT','GTCTT','GGATT','GGGTT','GGTTT','GGCAT','GGCCT','GGCGT','GGCTA','GGCTC','GGCTG']

    for read in local_lst:                                                          # iterate over junction analysis loop for each read in the raw reads list
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

        jxn = read[lpos:rpos+10]            # define junction as the sequence between the left match position and the right match position + 10 (accounting for the 10mer match length)
        # filter out reads with ambiguities in defined junction (e.g. N, R, Y, K, S, W, etc)
        ambiguity_jxns = 0                                                                  # open ambiguity jxn count at 0
        jxn_ambig = jxn.replace("A","").replace("C","").replace("G","").replace("T","")          # delete all ACGT bases, leave ambiguity calls behind
        if len(jxn_ambig) > 0:                                                                 # if no bases remain after removal of ACGT, write unambiguous sequence and corresponding analysis info to new lists for .csv output                                                                      # if any base calls reamin after removal of ACGT, they are ambiguities
            ambiguity_jxns += 1                                                             # count reflects reads that are removed from final analysisdue to one or more ambiguous base calls

        strand_ID = read[read.find('GATGCTCTCA')+10:read.find('GATGCTCTCA')+15]
        if strand_ID in top_IDs:
            strand_class = 'top'
            check_strand = 1
        elif strand_ID in bottom_IDs:
            strand_class = 'bottom'
            check_strand = 1
        else:
            strand_class = 'N/A'
            check_strand = 0

        randomer = read[read.find('GCGTA')+5:read.find('CTCGACAA')]
        rand_ambig = randomer.replace("A","").replace("C","").replace("G","").replace("T","")
        if len(randomer) == 10 and len(rand_ambig) == 0:
            tagged_jxn = ' -' + strand_class + '--' + jxn + '---' + randomer
            check_UMI = 10
        elif len(randomer) != 10 or len(rand_ambig) != 0:
            check_UMI = 0

        if check_strand == 1 and len(jxn_ambig) == 0 and check_UMI == 10:
            local_tags.append(tagged_jxn)
    return local_tags

if __name__ == '__main__':
    full_st = time.time()
    top_strands = []
    bottom_strands = []
    strand_ID = []
    barcodes = []
    frames = []
    consensus_seqs = []

    fGT = ['GT2-1b','GT2-2b','GT2-3b','GT30-1b','GT30-2b','GT30-3b']
    fTC = ['TC2-1b','TC2-2b','TC2-3b','TC30-1b','TC30-2b','TC30-3b']
    GAG2 = ['GAG2-1b','GAG2-2b','GAG2-3b']
    GAG30 = ['GAG30-1b','GAG30-2b','GAG30-3b']
    
    library = []
    initial_reads = []
    jxns_w_ambig = []
    failed_matches = []
    jxns_w_subs = []
    filtered_jxns = []
    compiled_jxns = []

    # start of junction analysis module
    for file in fTC:                       # swap between fGT and fTC to analyze different data sets
        print('start', file)
        t0 = time.time()
        df = pd.read_csv(file + '.csv', names = ["ID", "blank", "reads"])               # interpret .csv file with raw reads in column 3, this is where our .csv output lists reads, generated in CLC Genomics Workbench
        raw_reads = df['reads'].tolist()                                                # place raw reads from .csv file into list of reads

        # open lists of all outputs for csv
        # for junction analysis (reads with accurate left/right 10mer matches and no base substitutions will be placed here)
        UMIs = []
        UMI_reps = []
        UMI_IDs = []
        tagged_junctions = []
        ambiguity = []
        failed = 0                                                                      # open empty count of reads which fail to map either a left or right flank 10mer match
        final_jxns = []

        # for substitution miscall filter (reads from miscall/substitution filter that display evidence of a base substitution in the junction flanks (10 nt to either side) will be placed here)
        substitutions = 0                                                               # open empty count of junctions with evident base substitution
        umis = []
        tagged_reads = []
        key_lengths = []
        key_lengths2 = []
        t1 = 0
        t2 = 0
        t3 = 0

        umi_lists = {'AA':[],'AC':[],'AG':[],'AT':[],
                     'CA':[],'CC':[],'CG':[],'CT':[],
                     'GA':[],'GC':[],'GG':[],'GT':[],
                     'TA':[],'TC':[],'TG':[],'TT':[]
                  }

        for read in raw_reads:
            if len(read) >150:
                randomer = read[read.find('GCGTA')+5:read.find('CTCGACAA')]
                rand_ambig = randomer.replace("A","").replace("C","").replace("G","").replace("T","")
                if len(randomer) == 10 and len(rand_ambig) == 0:
                    key = str(randomer[:2])
                    umi_lists[key].append(read)
                    umis.append(randomer)
        library_length = len(umis)
        sorted_unique_UMIs_2 = sorted(list(set(umis)))
        print(file, 'UMIs:', len(umis))
        print(file, 'unique UMIs:', len(sorted_unique_UMIs_2))

        t1 = int(time.time() - t0)
        print('t1: ', t1, 'seconds')

        with mp.Pool(processes=mp.cpu_count()) as pool:
            args = [(partition, umi_lists[partition]) for partition in umi_lists]
            results = pool.starmap(process_read, args)

        for local_tags in results:
            tagged_junctions.extend(local_tags)

        t2 = int(time.time() - t0)
        print('t2: ', t2, 'seconds')

        d = {x:[] for x in umis}
        
        for query in tagged_junctions:
            if query.find('----') == -1:
                UMI_tag = query[query.find('---') + 3:]
                d[str(UMI_tag)].append(query)
        
        for keylen in d:
            read_count = len(d[keylen])
            UMI_reps.append(read_count)
        UMI_reps.sort(reverse=True)
        threshold = 0.02
        thresh_os = library_length*threshold/5000
        print(file, "thresh_os:",threshold*100,'% =', thresh_os)
        
        for i in d.copy():
            if len(d[i]) < thresh_os:
                d.pop(i)

        u = {x:[] for x in list(d.keys())}
        c = {x:[] for x in list(d.keys())}

        for key3 in d:
            tags = sorted(set(d[key3]))
            u[key3] = tags

        for key4 in d:
            key_lengths2.append(len(u[key4]))
        cap2 = max(key_lengths2)
        
        for key4 in u:
            while len(u[key4]) < cap2:
                u[key4].append('-')

        for key4 in d:
            tag_lst = sorted(set(d[key4]),reverse = True)
            temp = []
            i = 0
            con_seq = ''
            c[key4].append('x')
            c[key4].append('y')
            for tag in tag_lst:
                cnt = d[key4].count(tag)
                temp.append(cnt)
            max_strand = max(temp)

            for tag in tag_lst:
                cnt = d[key4].count(tag)
                if cnt >= 0.2*len(d[key4]):
                    c[key4].append(tag)
                    c[key4].append(str(cnt))
                    new_fill = tag[:tag.find('---')]
                    con_seq = con_seq + new_fill + '---'
            if len(con_seq) > 2:
                c[key4][0] = con_seq
                c[key4][1] = len(con_seq)
                consensus_seqs.append(con_seq)
            else:
                c.pop(key4)
                
        for key4 in c:
            while len(c[key4]) < 12:
                c[key4].append('')

        paired_df = pd.DataFrame(c)
        transposed_paired_df = paired_df.T
        transposed_paired_df.to_csv(file + '_2umi_20strand_uc.csv')
        t3 = int(time.time() - t0)
        print('t3: ', t3, 'seconds')

        fin_file = int(time.time() - t0)
        print('Time to run program: ', file, fin_file, 'seconds')

    unique_consensus_seqs = list(set(consensus_seqs))
    compiled = {'consensus_seqs':unique_consensus_seqs}

    for file in fTC:                        # swap between fGT and fTC to analyze different data sets
        consensus_counts = []
        consensus_freq = []
        file_data = pd.read_csv(file + '_2umi_20strand_uc.csv', usecols = ['0'],
                                keep_default_na = False)
        pulled_data = file_data['0'].tolist()
#        pulled_data = transposed_paired_df[0].tolist()
        
        for consensus in unique_consensus_seqs:
            count = pulled_data.count(consensus)
            consensus_counts.append(count)
            freq = count/len(pulled_data)
            consensus_freq.append(freq)

        col1 = file + '_counts'
        col2 = file + '_freq'
        compiled[col1] = consensus_counts
        compiled[col2] = consensus_freq

    all_files_df = pd.DataFrame(compiled)
    all_files_df.to_csv('fTC_2umi_20strand_freq.csv')                # change output file name according to file set analyzed (fGT or fTC)

    fin = int(time.time() - full_st)
    print('Time to run program: ', fin, 'seconds')
    print('fin fin')                                                                # identifies completion of full code process