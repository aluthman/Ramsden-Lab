import pandas as pd
def revcom(x):          #reverse complement of input sequence
    x_rev = x[::-1]     #reverse sequence
    x_revcom = x_rev.replace("A", "x").replace("C", "y").replace("G", "C").replace("T", "A").replace("x", "T").replace("y", "G")        #complement sequence
    return x_revcom     #feed revcom sequence back to code body

def MH_search (l,r,flank,ins):
    for ldel in range (0, 16):      #defines junction matrix for all ldels and rdels within 15 nt of breaksite
        for rdel in range(0, 16):
            i = 2                                       #start at MH size 8 (10-i)
            #code block for TINS in the left flank
            if inv_MH_len != -1 and flank == 'left':
                left_ins = l + ins                      #append insert sequence to left flank for resolving MH search
                right = 'NNNNNNNNNN' + r                #append 10 Ns upstream of right flank for resolving MH search

                ld = left_ins[len(left_ins) - 10 - ldel:len(left_ins) - ldel]       #define left search term as last 10 nts of left sequence
                rd = right[rdel:rdel + 10]                                          #define right search term as first 10 nts of right sequence

                while ld[i:] != rd[i:] and i <=8:                                   #for MH sizes 2-8 nts (i=8 iterating down to 2), if MH of ld doesn't match MH of rd
                    i += 1                                                          #decrease MH size by one
                if ld[i:] == rd[i:] and i <= 8:                                     #if MH of ld matches MH of rd, MH found of size 10 - i
                    left_match = left_flank[0:len(left_flank) - left_ss_deletion]   #remove deleted nts from left flank sequence
                    right_match = right_flank[0 + rdel:len(right_flank)]            #remove deleted nts from right flank sequence
                    MH_len = 10 - i                                                 #define MH size
                    MH_nt = rd[i:]                                                  #define MH sequence based on ld/rd match and MH size
                    insert_len = ins.find(MH_nt)                                    #locate resolving MH within the insert
                    insert = ins[:insert_len]                                       #redefine insert without resolving MH, as MH will not appear in final product
                    junction = left_match + insert + right_match                    #define junction with left and right deletions, and redefined insert

                    #append all calculated values to corresponding lists for output
                    if insert_len != -1:
                        MH_size.append(MH_len)
                        MH_ID.append(MH_nt)
                        rdel2 = rdel - MH_len                                           #redefine right deletion to exclude the size of the MH
                        left_flank_list.append(left_match)
                        right_flank_list.append(right_match)
                        rdel_size_MH.append(rdel2)
                        junction_list.append(junction)

                        seed_seq.append(left_seed_MH_revcom)
                        seed_location.append(left_seed_pos)
                        hairpins.append(left_hairpin)
                        seed_size.append(inv_MH_len)
                        template_flank.append("left")
                        TINS_seq.append(insert)
                        TINS_len.append(insert_len)

                        #separate list appending blocks for inv TINS and ss TINS
                        #changes to deletion size in the TINS flank and insert identity
                        #for inv TINS, TINS flank deletion corresponds to location of priming MH
                        #for ss TINS, TINS flank deletion corresponds to location of priming MH match deeper in flank (MH location+hairpin+MH size)
                        if ins_ID == "inv":
                            TINS_class.append("inv")
                            ldel_size.append(left_seed_pos)
                            tdel = rdel2 + left_seed_pos
                            tdel_size.append(tdel)
                        if ins_ID == "ss":
                            TINS_class.append("ss")
                            ldel_size.append(left_ss_deletion)
                            tdel = rdel2 + left_ss_deletion
                            tdel_size.append(tdel)
                        break

            #code block for TINS in the right flank
            #same code block as left flank TINS, with modification outlined...
            #left and right flanks are reverse complemented, so appending to ldel and rdel lists needs to be swapped...
            #and inserts need to be reverse complemented to read the sequenced strand
            if inv_MH_len != -1 and flank == 'right':
                left_ins = l + ins
                right = 'NNNNNNNNNN' + r

                ld = left_ins[len(left_ins) - 10 - ldel:len(left_ins) - ldel]
                rd = right[rdel:rdel + 10]

                while ld[i:] != rd[i:] and i <=8:
                    i += 1
                if ld[i:] == rd[i:] and i <= 8:
                    left_match = left_flank[0:len(left_flank) - left_ss_deletion]
                    right_match = right_flank[0 + rdel:len(right_flank)]
                    MH_len = 10 - i
                    MH_nt = rd[i:]
                    insert_len = ins.find(MH_nt)
                    insert = ins[:insert_len]
                    insert_revcom = revcom(insert)                      #reverse complement of insert to call from the sequenced strand
                    junction = left_match + insert_revcom + right_match

                    if insert_len != -1:
                        MH_size.append(MH_len)
                        MH_nt_revcom = revcom(MH_nt)
                        MH_ID.append(MH_nt_revcom)
                        rdel2 = rdel - MH_len
                        left_flank_list.append(left_match)
                        right_flank_list.append(right_match)
                        ldel_size.append(rdel2)                             #flip left and right deletions for output to account for reverse complement of input flanks
                        junction_list.append(junction)

                        seed_seq.append(left_seed_MH_revcom)
                        seed_location.append(left_seed_pos)
                        hairpins.append(left_hairpin)
                        seed_size.append(inv_MH_len)
                        template_flank.append("right")
                        TINS_seq.append(insert_revcom)
                        TINS_len.append(insert_len)

                        if ins_ID == "inv":
                            TINS_class.append("inv")
                            rdel_size_MH.append(left_seed_pos)              #flip left and right deletions for output to account for reverse complement of input flanks
                            tdel = rdel2 + left_seed_pos
                            tdel_size.append(tdel)
                        if ins_ID == "ss":
                            TINS_class.append("ss")
                            rdel_size_MH.append(left_ss_deletion)           #flip left and right deletions for output to account for reverse complement of input flanks
                            tdel = rdel2 + left_ss_deletion
                            tdel_size.append(tdel)
                        break

            #code block for MH.deletion products, where inv_MH_len was previously established as null (-1)
            #same code block as left flank TINS...
            #but there is no insertion, so MH search occurs in native left and right flanks
            if inv_MH_len == -1:
                left = l
                right = 'NNNNNNNNNN' + r

                ld = left[len(left) - 10 - ldel:len(left) - ldel]
                rd = right[rdel:rdel + 10]

                while ld[i:] != rd[i:] and i <=8:
                    i += 1
                if ld[i:] == rd[i:] and i <= 8:
                    tdel = ldel + rdel
                    left_match = left_flank[0:len(left_flank) - ldel]
                    right_match = right_flank[0 + rdel:len(right_flank)]
                    junction = left_match + right_match
                    junction_list.append(junction)
                    MH_len = 10 - i
                    MH_nt = rd[i:]

                    MH_size.append(MH_len)
                    MH_ID.append(MH_nt)
                    rdel2 = rdel - MH_len
                    ldel_size.append(ldel)
                    tdel_size.append(tdel)
                    left_flank_list.append(left_match)
                    right_flank_list.append(right_match)
                    rdel_size_MH.append(rdel2)
                    seed_seq.append("MH.del")
                    seed_location.append("")
                    hairpins.append("")
                    seed_size.append("")
                    template_flank.append("")
                    TINS_seq.append("")
                    TINS_len.append(0)
                    TINS_class.append("MH.del")

locus_name = input ("Enter target locus name:")             #request locus input, if locus name matches any of our common loci below, it will use that sequence

#R26-TINS (R26E); note orientation of R26-TINS, sequence reads in csv input files are the strand containing the 5-CCN-3 PAM
if locus_name == 'R26E':
    locus = 'TGGCTTATCCAACCCCTAGACAGAGCATTGGCATTTTCCCTTTCCTGATCTTAGAAGTCTGATGACTCATGAAACCAGACAGATTAGTTACATACACCACAAATCGAGGCTGTAGCTGGGGCCTCAACACTGCAGTTCTTTTATAACTCCTTAGTACACTTTTTGTTGATCCTTTGCCTTGATCCTTAATTTTCAGTGTCTATCACCTCTCCCGTCAGGTGGTGTTCCACA'
    sgRNA = 'CCACAAATCGAGGCTGTAGCTGG'
    y = locus.find (sgRNA)
    left_flank = locus[0:y + 6]
    right_flank = locus[y + 6:len(locus)]

#R26-TINS (R26E_v2); note orientation of R26-TINS, sequence reads in csv input files are the strand containing the 5-CCN-3 PAM
if locus_name == 'R26E_v2':
    locus = 'TGGCTTATCCAACCCCTAGACAGAGCATTGGCATTTTCCCTTTCCTGATCTTAGAAGTCTGATGACTCATGAAACCAGACAGATTAGTTACATACACCACAAATCGAGGCTGTAGCTGGGGCCTCAACACTGCAGTTCTTTTATAACTCCTTAGTACACTTTTTGTTGATCCTTTGCCTTGATCCTTAATTTTCAGTGTCTATCACCTCTCCCGTCAGGTGGTGTTCCACA'
    sgRNA = 'CCTCAACACTGCAGTTCTTTTAT'
    y = locus.find (sgRNA)
    left_flank = locus[0:y + 6]
    right_flank = locus[y + 6:len(locus)]

#R26D; note orientation of R26-TINS, sequence reads in csv input files are the strand containing the 5-CCN-3 PAM
elif locus_name == 'R26D':
    locus = 'TATTGCTTATCTTTAGTTCCCAACACTTGGGAGGCAGAGGCCAGCCAGGGCTATGTGACAAAAACCTTGTCTAGAGGAGAAACTTCATAGCTTATTTCCTATTCACGTAACCAGGTTAGCAAAATTTACCAGCCAGAGATGAAGCTAACAGTGTCCACTATATTTGTAGTGTTTTAAGTCAATTTTTTAAATATACTTAATAGAATTAAAGCTATGGTGAACCAAGTACAA'
    sgRNA = 'CCTATTCACGTAACCAGGTTAGC'
    y = locus.find (sgRNA)
    left_flank = locus[0:y + 6]
    right_flank = locus[y + 6:len(locus)]

#R26-MHD (R26A)
elif locus_name == 'R26A':
    locus = "GGGCCTGGGAGAATCCCTTCCCCCTCTTCCCTCGTGATCTGCAACTCCAGTCTTTCTAGAAGATGGGCGGGAGTCTTCTGGGCAGGCTTAAAGGCTAACCTGGTGTGTGGGCGTTGTCCT"
    sgRNA = "ACTCCAGTCTTTCTAGAAGATGG"
    cutsite = locus.find(sgRNA) + 17
    left_flank = locus[0:locus.find(sgRNA)+17]
    right_flank = locus[locus.find(sgRNA)+17:len(locus)]

#MEF204 (LBR locus)
elif locus_name == 'LBR':
    locus = "ACCCGAGGGGCACTTAGGGAGCAGGGGGGCAGTGAACACCTCTGCATGAGCAGGGGCATAAAAACGGAAGGTGACAACATCCAGCTTTGCCTATTTCAAC"
    sgRNA = "GAACACCTCTGCATGAGCAG"
    cutsite = locus.find(sgRNA) + 17
    left_flank = locus[0:locus.find(sgRNA)+17]
    right_flank = locus[locus.find(sgRNA)+17:len(locus)]

else:                                                       #manually sequence input, if using locus other than those predefined above
    left_flank = input("Enter left flank sequence:")        #manual entry of left and right flank
    right_flank = input ("Enter right flank sequence:")
    #locus = input ("Enter locus sequence:")                #manual entry of locus sequence and gRNA sequence, code will identify cutsite from gRNA location
    #sgRNA = input ("Enter sgRNA sequence (no PAM, but in NGG orientation):")
    #cutsite = locus.find(sgRNA) + 17
    #left_flank = locus[0:locus.find(sgRNA)+17]             #code automatically defines left and right flanks if locus and gRNA sequence are provided
    #right_flank = locus[locus.find(sgRNA)+17:len(locus)]

left_revcom = revcom(left_flank)            #establish left and right flank reverse complements for inverted TINS searches
right_revcom = revcom(right_flank)

#MH deletion lists, open empty lists for populating
ldel_size = []
rdel_size = []
tdel_size = []
junction_list = []
left_flank_list = []
right_flank_list = []
MH_size = []
MH_ID = []
rdel_size_MH = []

#TINS predictor lists, open empty lists for populating
seed_seq = []
seed_location = []
hairpins = []
seed_size = []
template_flank = []
TINS_seq = []
TINS_len = []
TINS_class = []

#repetitive resolutions filter, open empty lists for populating
left_repetitive_resolutions = []
right_repetitive_resolutions = []
junction_family = []
x_position = []
y_position = []
unique_family = []

#search for MH deletions
inv_MH_len = -1                                 #unless otherwise noted (for TINS) establish null insertion length (applied to MH.deletion products)
MH_search(left_flank,right_flank,'','')         #search for resolving MHs, no insertion (canonical TMEJ MH deletion products)

#search for inv.ss TINS
#left flank TINS identification
m = 6                                                           #initial MH search size of 6...iterating down to 2
while m >= 2:                                                   #only consider MHs 2 or larger, loop stops after searching for MHs of size 2 nt
    for i in range (0,16):                                      #search first 16 nts in each flank
        inv_MH_len = m                                          #iterate TINS priming MH from 6 down to 2

        left_seed_MH = left_revcom[i:i+m]                           #define MH sequence to search for resolving MH match
        left_seed_MH_revcom = revcom(left_seed_MH)                  #reverse complement of MH to search for inverted sequence match
        left_seed_pos = i                                           #position of the MH currently used for search
        left_search = left_revcom[(i + m + 3):]                     #restrain search to sequence beyond the MH in question, including a minimum hairpin of 3 nt
        left_ss_deletion = left_search.find(left_seed_MH_revcom)    #define deletion size for ss TINS, based on MH position, length, and hairpin size
        left_hairpin = 0                                            #initially define hairpin as 0, will be redefined if inv MH match exists
        if left_ss_deletion != -1:                                  #if revcom MH is located in flank search...
            left_ss_deletion = left_ss_deletion + i + m + 3         #redefine left ss TINS deletion based on MH size and position
            left_hairpin = left_ss_deletion - (left_seed_pos + m)   #redefine hairping size as distance between the MH matches in same flank

        #inverted TINS have two flavors, classical inverted (smaller deletion)...
        #and strand switching ('ss', relying on cleavage of the hairpin intermediate and resulting in a larger deletion)
        if ((left_ss_deletion != -1) and (12 >= left_hairpin >= 3)):                #limit left hairpin length
            ins_ID = "inv"                                                          #first iteration, searching for inverted TINS
            left_flank_inv_del_revcom = left_revcom[left_seed_pos:]                 #define new left flank, accounting for deletion up to the priming MH
            left_flank_inv_del = revcom(left_flank_inv_del_revcom)                  #reverse complement of deleted left flank to consider 'top' strand for resolving MH search
            left_inv_full = left_revcom[left_ss_deletion + inv_MH_len:left_ss_deletion + inv_MH_len + 15]       #take 15 nts of inv insertion
            MH_search(left_flank_inv_del,right_flank,'left',left_inv_full)          #begin MH search using inferred insertion according to classical inverted TINS mechanism

            ins_ID = "ss"                                                           #second iteration, searching for ss TINS
            left_flank_ss_del_revcom = left_revcom[left_ss_deletion:]               #define new left flank, accounting for deletion up to the priming MH + MH size + hairpin
            left_flank_ss_del = revcom(left_flank_ss_del_revcom)                    #reverse complement of deleted left flank to consider 'top' strand for resolving MH search
            left_ss_full = left_revcom[left_seed_pos + inv_MH_len:left_seed_pos + inv_MH_len + 15]              #take 15 nts of ss insertion
            MH_search(left_flank_ss_del,right_flank,'left',left_ss_full)            #begin MH search using inferred insertion according to strand switching TINS mechanism
    m -= 1                                                                          #iterate to smaller MH definition if match IS NOT found, loop stops when MH size = 2

#right flank TINS identification, same code block as left flank TINS identification with modifications outlined...
#swapped left and right flank positions and reverse complements to use same resolving MH finder code black (function = MH_search)
m = 6
while m >= 2:
    for i in range (0,16):
        inv_MH_len = m
        r_left_revcom = right_flank                             #swap left and right flanks to operate on bottom strand, thus MH_search function logic holds true
        r_right_flank = left_revcom

        left_seed_MH = r_left_revcom[i:i+m]
        left_seed_MH_revcom = revcom(left_seed_MH)
        left_seed_pos = i
        left_ss_deletion = r_left_revcom.find(left_seed_MH_revcom)
        left_hairpin = left_ss_deletion - (left_seed_pos + m)

        if ((left_ss_deletion != -1) and (12 >= left_hairpin >= 3)):               #limit right hairpin length
            ins_ID = "inv"
            left_flank_inv_del_revcom = r_left_revcom[left_seed_pos:]
            left_flank_inv_del = revcom(left_flank_inv_del_revcom)
            left_inv_full = r_left_revcom[left_ss_deletion + inv_MH_len:left_ss_deletion + inv_MH_len + 15]     #take 15 nts of inv insertion
            MH_search(left_flank_inv_del,r_right_flank,'right',left_inv_full)

            ins_ID = "ss"
            left_flank_ss_del_revcom = r_left_revcom[left_ss_deletion:]
            left_flank_ss_del = revcom(left_flank_inv_del_revcom)
            left_ss_full = r_left_revcom[left_seed_pos + inv_MH_len:left_seed_pos + inv_MH_len + 15]                      #take 15 nts of ss insertion
            MH_search(left_flank_ss_del,r_right_flank,'right',left_ss_full)
    m -= 1

#establish identifiers for each resolved TINS junction, will be used to filter out repetitive resolutiosn
for z in range(0,len(junction_list)):
    if TINS_class[z] == "MH.del":                                                           #for MH deletion junctions
        lrep = str(ldel_size[z]) + "." + seed_seq[z] + "." + MH_ID[z]                       #for MH deletions, filtering will be based on MH seq, retaining only the most proximal MH resolution to the breaksite
        rrep = str(rdel_size_MH[z]) + "." + seed_seq[z] + "." + MH_ID[z]
        family = "." + seed_seq[z] + "."                                                    #define family grouping for clustered plotting in R
    else:                                                                               #for inv and ss TINS, filtering will be based on both priming MH and resolving MH
        lrep = str(ldel_size[z]) + "." + seed_seq[z] + "." + TINS_seq[z]                    #define repetitive identifiers for left and right flank TINS
        rrep = str(rdel_size_MH[z]) + "." + seed_seq[z] + "." + TINS_seq[z]
        family = str(template_flank[z]) + "." + seed_seq[z] + "." + str(seed_location[z])   #define family grouping for clustered plotting in R
    left_repetitive_resolutions.append(lrep)                                                #establish left repetition list for filtering
    right_repetitive_resolutions.append(rrep)                                               #establish right repetition list for filtering
    junction_family.append(family)                                                          #establish family classifications for clustered plotting in R

xmax = max(ldel_size)                                                               #define maximum x axis position to place labels beyond data set on R plot
for family_ID in junction_family:
    if family_ID not in unique_family:                                              #populate list of unique family classifications
        unique_family.append(family_ID)

for y in range(0,len(unique_family)):                                               #set unique family list for identification in R plot
    for z in range(0,len(junction_list)):
        if junction_family[z] == unique_family[y]:
            if ((template_flank[z] == 'left') or (template_flank[z] == '')):        #keep left flank TINS and MH.dels in left label column
                x_position.append(xmax + 5)                                         #set label x.position 5 units past largest ldel
                y_position.append(y*2)                                              #set label y.position in ladder of 2s for MH.del and left flank TINS
            if template_flank[z] == 'right':                                        #move right flank TINS to right label column
                x_position.append(xmax + 12)                                        #set label x.position 7 units past prior x.position
                yr = y*2 - max(y_position)                                          #reset y.position in ladder of 2s, restarting at 2 for right flank TINS
                y_position.append(yr)

data = {'ldel': ldel_size, 'rdel': rdel_size_MH, 'MH.size': MH_size, 'resolve.MH.seq': MH_ID, 'tdel': tdel_size,        #establish data frame for output
         'junction': junction_list, 'left flank': left_flank_list, 'right flank': right_flank_list,
         'template.flank':template_flank, 'seed.MH.seq': seed_seq, 'seed.MH.location': seed_location,
         'hairpin.distance': hairpins, 'seed.MH.size': seed_size, 'TINS.seq': TINS_seq, 'TINS.len':TINS_len, 'TINS.class':TINS_class,
         'left.repetitive': left_repetitive_resolutions, 'right.repetitive': right_repetitive_resolutions,
         'junction.family': junction_family, 'x.position': x_position, 'y.position': y_position}

output_df = pd.DataFrame(data)
output_df.drop_duplicates(subset ="junction", keep = 'first', inplace = True)               # removes duplicates from unique_df based on junction
output_df.drop_duplicates(subset ="left.repetitive", keep = 'first', inplace = True)        # removes duplicates from unique_df based on repetitive resolutions
output_df.drop_duplicates(subset ="right.repetitive", keep = 'first', inplace = True)       # removes duplicates from unique_df based on repetitive resolutions
unique_junction_df = output_df.reset_index(drop = True)                                     # reset index after removing duplicates

unique_junction_df.to_csv('TMEJ_predictor_' + locus_name + '.csv')                          #output final data matrix to .csv