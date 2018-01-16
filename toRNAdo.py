import numpy
from copy import deepcopy
import os
from sys import argv
import csv

# filters a dictionary by length
def filter_by_length(unfiltered_list):
        copy1 = deepcopy(unfiltered_list)
        for nuc1 in copy1:
                count1 = 1
                position1 = nuc1 + 1
                while position1 in copy1:
                        count1 += 1
                        position1 += 1
                else:
                        if count1 < 50 and nuc1 - 1 not in unfiltered_list:
                                del unfiltered_list[nuc1]

# function to create text files
def write_to_file(sequences, filename):
        out_file = open(folder + "/toRNAdo_output/" + filename + '.wig','w')
        list2 = sorted(sequences.keys())
        lines1 = ["%s\t%f" % (g, sequences[g]) for g in list2]
        lines1.insert(0, "variableStep\tchrom=Rd")
        out_file.write('\n'.join(lines1))
        out_file.write('\n')
        out_file.close()    

# function to find ncRNAs to the "right" of UTR, and reclassify it
def find_RNA_right(UTR_dict):
        temp_dict = {}
        another_temp_dict = {}
        copy_UTR = deepcopy(UTR_dict)
        for nuc in copy_UTR:
                if nuc - 1 not in copy_UTR:
                        while nuc + 1 in copy_UTR and copy_UTR[nuc + 1] <= copy_UTR[nuc]:
                                nuc += 1
                        else:
                                if nuc + 1 not in copy_UTR:
                                        None
                                else:
                                        temp_thresh = copy_UTR[nuc]
                                        temp_dict[nuc] = copy_UTR[nuc]
                                        while nuc + 1 in copy_UTR:
                                                temp_dict[nuc + 1] = copy_UTR[nuc + 1]
                                                nuc += 1
                                        else:
                                                if any((x / 5) > temp_thresh for x in temp_dict.values()):
                                                        for nuc1 in temp_dict:
                                                                if UTR_dict == UTR5_minus:
                                                                        RNA_minus[nuc1] = temp_dict[nuc1]
                                                                        del UTR_dict[nuc1]
                                                                if UTR_dict == UTR3_plus:
                                                                        RNA_plus[nuc1] = temp_dict[nuc1]
                                                                        del UTR_dict[nuc1]
                                                else:
                                                        min_key = min(temp_dict, key=lambda k: temp_dict[k])
                                                        min_thresh = temp_dict[min_key]
                                                        pos = min_key + 1
                                                        another_temp_dict[min_key] = temp_dict[min_key]
                                                        while pos in temp_dict:
                                                                another_temp_dict[pos] = temp_dict[pos]
                                                                pos += 1
                                                        else:
                                                                if any((x1 / 5) > min_thresh for x1 in another_temp_dict.values()):
                                                                        for nuc2 in another_temp_dict:
                                                                                if UTR_dict == UTR5_minus:
                                                                                        RNA_minus[nuc2] = another_temp_dict[nuc2]
                                                                                        del UTR_dict[nuc2]
                                                                                if UTR_dict == UTR3_plus:
                                                                                        RNA_plus[nuc2] = another_temp_dict[nuc2]
                                                                                        del UTR_dict[nuc2]
                                                        another_temp_dict.clear()
                                        temp_dict.clear()

# function to find ncRNA to the "left" of UTR, and reclassify it
def find_RNA_left(UTR_dict1):
        temp_dict = {}
        another_temp_dict = {}
        copy_UTR = deepcopy(UTR_dict1)
        for nuc in copy_UTR:
                if nuc + 1 not in copy_UTR:
                        while nuc - 1 in copy_UTR and copy_UTR[nuc - 1] <= copy_UTR[nuc]:
                                nuc -= 1
                        else:
                                if nuc - 1 not in copy_UTR:
                                        None
                                else:
                                        temp_thresh = copy_UTR[nuc]
                                        temp_dict[nuc] = copy_UTR[nuc]
                                        while nuc - 1 in copy_UTR:
                                                temp_dict[nuc - 1] = copy_UTR[nuc - 1]
                                                nuc -= 1
                                        else:
                                                if any((x / 5) > temp_thresh for x in temp_dict.values()):
                                                        for nuc1 in temp_dict:
                                                                if UTR_dict1 == UTR3_minus:
                                                                        RNA_minus[nuc1] = temp_dict[nuc1]
                                                                        del UTR_dict1[nuc1]
                                                                if UTR_dict1 == UTR5_plus:
                                                                        RNA_plus[nuc1] = temp_dict[nuc1]
                                                                        del UTR_dict1[nuc1]
                                                else:
                                                        min_key = min(temp_dict, key=lambda k: temp_dict[k])
                                                        min_thresh = temp_dict[min_key]
                                                        pos = min_key - 1
                                                        another_temp_dict[min_key] = temp_dict[min_key]
                                                        while pos in temp_dict:
                                                                another_temp_dict[pos] = temp_dict[pos]
                                                                pos -= 1
                                                        else:   
                                                                if any((x1 / 5) > min_thresh for x1 in another_temp_dict.values()):
                                                                        for nuc2 in another_temp_dict:
                                                                                if UTR_dict1 == UTR3_minus:
                                                                                        RNA_minus[nuc2] = another_temp_dict[nuc2]
                                                                                        del UTR_dict1[nuc2]
                                                                                if UTR_dict1 == UTR5_plus:
                                                                                        RNA_plus[nuc2] = another_temp_dict[nuc2]
                                                                                        del UTR_dict1[nuc2]
                                                        another_temp_dict.clear()
                                        temp_dict.clear()

# function to find ncRNA in the intergenic region, which has been originally classified as part of an operon
def find_RNA_operon(operon_dict):
        temp_dict = {}
        temp_dict_right = {}
        temp_dict_left = {}
        new_left = {}
        new_right = {}
        copy_operon = deepcopy(operon_dict)
        for nuc in copy_operon:
                if nuc - 1 not in copy_operon:
                        while nuc + 1 in copy_operon:
                                temp_dict[nuc] = copy_operon[nuc]
                                nuc += 1
                        else:
                                temp_dict[nuc] = copy_operon[nuc]
                                min_key = min(temp_dict, key=lambda k: temp_dict[k])
                                min_thresh = temp_dict[min_key]
                                if any((x / 5) > min_thresh for x in temp_dict.values()):                                        
                                        while min_key + 1 in temp_dict:
                                                temp_dict_right[min_key + 1] = temp_dict[min_key + 1]
                                                min_key += 1
                                        else:
                                                min_key = min(temp_dict, key=lambda k: temp_dict[k])
                                                temp_dict_left[min_key] = temp_dict[min_key]
                                                while min_key - 1 in temp_dict:
                                                        temp_dict_left[min_key - 1] = temp_dict[min_key - 1]
                                                        min_key -= 1
                                                else:                                                       
                                                        if any((x1 / 5) > min_thresh for x1 in temp_dict_left.values()):
                                                                max_key_left = max(temp_dict_left, key=lambda k1: temp_dict_left[k1])
                                                                max_key_left_thresh = temp_dict_left[max_key_left]
                                                                while max_key_left - 1 in temp_dict_left:
                                                                        new_left[max_key_left - 1] = temp_dict_left[max_key_left - 1]
                                                                        max_key_left -= 1
                                                                else:
                                                                        if any((max_key_left_thresh / 5) > x2 for x2 in new_left.values()):
                                                                                min_key_left = min(new_left, key=lambda k2: new_left[k2])
                                                                                while min_key_left + 1 in temp_dict_left:
                                                                                        if operon_dict == operon_plus:
                                                                                                RNA_plus[min_key_left + 1] = temp_dict_left[min_key_left + 1]
                                                                                                del operon_dict[min_key_left + 1]
                                                                                        if operon_dict == operon_minus:
                                                                                                RNA_minus[min_key_left + 1] = temp_dict_left[min_key_left + 1]
                                                                                                del operon_dict[min_key_left + 1]
                                                                                        min_key_left += 1    
                                                        if any((x3 / 5) > min_thresh for x3 in temp_dict_right.values()):
                                                                max_key_right = max(temp_dict_right, key=lambda k3: temp_dict_right[k3])
                                                                max_key_right_thresh = temp_dict_right[max_key_right]
                                                                while max_key_right + 1 in temp_dict_right:
                                                                        new_right[max_key_right + 1] = temp_dict_right[max_key_right + 1]
                                                                        max_key_right += 1
                                                                else:
                                                                        if any((max_key_right_thresh / 5) > x4 for x4 in new_right.values()):
                                                                                min_key_right = min(new_right, key=lambda k4: new_right[k4])
                                                                                while min_key_right - 1 in temp_dict_right:
                                                                                        if operon_dict == operon_plus:
                                                                                                RNA_plus[min_key_right - 1] = temp_dict_right[min_key_right - 1]
                                                                                                del operon_dict[min_key_right - 1]
                                                                                        if operon_dict == operon_minus:
                                                                                                RNA_minus[min_key_right - 1] = temp_dict_right[min_key_right - 1]
                                                                                                del operon_dict[min_key_right - 1]
                                                                                        min_key_right -= 1                                                                                        
                                                        new_left.clear()
                                                        new_right.clear()
                                        temp_dict_right.clear()
                                        temp_dict_left.clear()                                                        
                        temp_dict.clear()

# function to reclassify former operon regions into UTRs and delete those regions from the operon dictionary
def operon_to_UTR(operon_dict1):
        temp_dict = {}
        operon_copy = deepcopy(operon_dict1)
        for nuc in operon_copy:
                if nuc - 1 not in operon_copy:
                        if operon_dict1 == operon_plus:
                                if nuc - 1 not in dict_strand:
                                        UTR5_plus[nuc] = operon_copy[nuc]
                                        del operon_dict1[nuc]
                                        while nuc + 1 in operon_copy:
                                                UTR5_plus[nuc + 1] = operon_copy[nuc + 1]
                                                del operon_dict1[nuc + 1]
                                                nuc += 1
                                else:
                                        while nuc + 1 in operon_copy:
                                                nuc += 1
                                        else:
                                                if nuc + 1 not in dict_strand:
                                                        UTR3_plus[nuc] = operon_copy[nuc]
                                                        del operon_dict1[nuc]
                                                        while nuc - 1 in operon_copy:
                                                                UTR3_plus[nuc - 1] = operon_copy[nuc - 1]
                                                                del operon_dict1[nuc - 1]
                                                                nuc -= 1
                        if operon_dict1 == operon_minus:
                                if nuc - 1 not in dict_strand:
                                        UTR3_minus[nuc] = operon_copy[nuc]
                                        del operon_dict1[nuc]
                                        while nuc + 1 in operon_copy:
                                                UTR3_minus[nuc + 1] = operon_copy[nuc + 1]
                                                del operon_dict1[nuc + 1]
                                                nuc += 1
                                else:
                                        while nuc + 1 in operon_copy:
                                                nuc += 1
                                        else:
                                                if nuc + 1 not in dict_strand:
                                                        UTR5_minus[nuc] = operon_copy[nuc]
                                                        del operon_dict1[nuc]
                                                        while nuc - 1 in operon_copy:
                                                                UTR5_minus[nuc - 1] = operon_copy[nuc - 1]
                                                                del operon_dict1[nuc - 1]
                                                                nuc -= 1
                                                                
# function to correct intergenic RNA into antisense of border RNAs
def RNA_to_antisense_border(rna_dict):
        RNA_copy = deepcopy(rna_dict)
        temp_dict = {}
        anti_dict = {}
        for rna in RNA_copy:
                if rna - 1 not in RNA_copy:
                        temp_dict[rna] = RNA_copy[rna]
                        while rna + 1 in RNA_copy:
                                temp_dict[rna + 1] = RNA_copy[rna + 1]
                                rna += 1
                        else:
                                if any(x in dict_strand for x in temp_dict.keys()):
                                        if any(y not in dict_strand for y in temp_dict.keys()):
                                                for rna1 in temp_dict:
                                                        if rna_dict == RNA_plus:
                                                                border_plus[rna1] = temp_dict[rna1]
                                                                del RNA_plus[rna1]
                                                        if rna_dict == RNA_minus:
                                                                border_minus[rna1] = temp_dict[rna1]
                                                                del RNA_minus[rna1]
                                        else:
                                                for rna2 in temp_dict:
                                                        if rna_dict == RNA_plus:
                                                                antisense_of_minus[rna2] = temp_dict[rna2]
                                                                del RNA_plus[rna2]
                                                        if rna_dict == RNA_minus:
                                                                antisense_of_plus[rna2] = temp_dict[rna2]
                                                                del RNA_minus[rna2]
                        temp_dict.clear()
                        
# function to look at all ncRNA dictionaries and find RNAs with a "5x" expression peak
def find_peaks(RNAdict, cooldict):
        temp_dict = {}
        copy_RNAdict = deepcopy(RNAdict)
        for abc in copy_RNAdict:
                if abc - 1 not in copy_RNAdict:
                        paul = abc
                        temp_dict[paul] = copy_RNAdict[paul]
                        while paul + 1 in copy_RNAdict:
                                temp_dict[paul + 1] = copy_RNAdict[paul]
                                paul += 1
                        else:
                                min_key = min(temp_dict, key=lambda k: temp_dict[k])
                                min_thresh = temp_dict[min_key]
                                if any((x / 5) > min_thresh for x in temp_dict.values()):
                                        for abcd in temp_dict:
                                                cooldict[abcd] = temp_dict[abcd]
                                                del RNAdict[abcd]
                        temp_dict.clear()

# file1 - nucleotide coverage for both strands; file2 - nucleotide coverage for minus strand; file 3 - nucleotide coverage for plus strand; file 4 - nucleotide coverage for all annotated genome features and their coordinates.
script, folder = argv

file1 = []
file2 = []
file3 = []
file4 = []

for filename in os.listdir(folder):
    if filename.endswith("nboth1.txt"):
        file1 = filename
    if filename.endswith("nminus1.txt"):
        file2 = filename
    if filename.endswith("nplus1.txt"):
        file3 = filename
    if filename.endswith("ncoverage1.txt"):
        file4 = filename

# open all files needed for analysis
data_both = open(folder + "/" + file1).readlines()
data_minus = open(folder + "/" + file2).readlines()
data_plus = open(folder + "/" + file3).readlines()
data_genes = open(folder + "/" + file4).readlines()

# empty lists that will have data appended to
realdata_minus = []
realdata_plus = []
realdata_genes = []
realdata_both = []
gene_ncoverage = []
both_ncoverage = []

# loops for appending into empty lists and changing lists of strings into lists of integers
for line_both in data_both:
        realdata_both.append(line_both.strip().split())
realdata_both = [[float(ab) for ab in bb] for bb in realdata_both]  # turns everything in a list of lists into integers
cov_both_list = [a[1] for a in realdata_both] # makes new list with just the coverage
sum_both = float(sum(cov_both_list)) # produces the sum of all coverage values, used for normalization later

for line_minus in data_minus:
        realdata_minus.append(line_minus.strip().split())
realdata_minus = [[float(ab) for ab in bb] for bb in realdata_minus] # turns everything in a list of lists into integers
for i in realdata_minus:                       # a loop to normalize each data based on the sum of coverage values calculated earlier
        i[1] = i[1] / sum_both * 10000000000

for line_plus in data_plus:
        realdata_plus.append(line_plus.strip().split())
realdata_plus = [[float(ab) for ab in bb] for bb in realdata_plus]
for e in realdata_plus:
        e[1] = e[1] / sum_both * 10000000000

for line_genes in data_genes:
        realdata_genes.append(line_genes.strip().split())
        
for value in realdata_genes:
        gene_ncoverage.append(value[3])   # makes a new list with just the gene coverage values. To calculate mean and standard deviation below. Not really used here yet!!
gene_ncoverage = [float(b) for b in gene_ncoverage]

threshold = 100.0 # expression threshold

# modifies a list of gene coordinates, removing the nucleotide column and
# turning the start position in the first column into nucleotide coordinate
for f in realdata_genes:
        f[0] = int(f[0])
        f[2] = int(f[2])
        f[3] = int(f[3])
        f[0] = f[0] + f[2] - 1
        f.pop(2)

# Dictionaries!! Easier to work with here than with lists...
dict_minus = {}
dict_plus = {}
dict_both = {}
dict_strand = {}

# Turning lists into dictionaries
for line1 in realdata_minus:
        dict_minus[line1[0]] = line1[1]
for line2 in realdata_plus:
        dict_plus[line2[0]] = line2[1]
for line3 in realdata_both:
        dict_both[line3[0]] = line3[1]
for line4 in realdata_genes:
        dict_strand[line4[0]] = line4[1]

# empty dictionaries
RNA_minus = {}
RNA_plus = {}
UTR3_minus = {}
UTR5_minus = {}
UTR3_plus = {}
UTR5_plus = {}
antisense_of_minus = {}
antisense_of_plus = {}
operon_minus = {}
operon_plus = {}

for key1 in dict_both:
        if key1 in dict_strand:
                if dict_strand[key1] == "-":
                        #UTR3_minus - puts all 3' UTR regions above the threshold on a minus strand into a new dictionary
                        if key1 - 1 not in dict_strand and key1 - 1 in dict_both:
                                while key1 != 1.0 and dict_minus[key1 - 1] > threshold and key1 - 1 not in dict_strand:
                                        UTR3_minus[key1 - 1] = dict_minus[key1 - 1]
                                        key1 -= 1
                elif dict_strand[key1] == "+":
                        #UTR5_plus - puts all 5' UTR regions above the threshold on a plus strand into a new dictionary
                        if key1 - 1 not in dict_strand and key1 - 1 in dict_both:
                                    while key1 != 1.0 and dict_plus[key1 - 1] > threshold and key1 - 1 not in dict_strand:
                                        UTR5_plus[key1 - 1] = dict_plus[key1 - 1]
                                        key1 -= 1
            
for key11 in dict_both:
        if key11 in dict_strand:
                if dict_strand[key11] == "-":  
                        #UTR5_minus - puts all 5' UTR regions above the threshold on a minus strand into a new dictionary
                        if key11 + 1 not in dict_strand and key11 + 1 in dict_both:
                                while dict_minus[key11 + 1] > threshold and key11 + 1 not in dict_strand:
                                        UTR5_minus[key11 + 1] = dict_minus[key11 + 1]
                                        key11 += 1
                elif dict_strand[key11] == "+":                     
                        #UTR3_plus - puts all 3' UTR regions above the threshold on a plus strand into a new dictionary
                        if key11 + 1 not in dict_strand and key11 + 1 in dict_both:
                                while dict_plus[key11 + 1] > threshold and key11 + 1 not in dict_strand:
                                        UTR3_plus[key11 + 1] = dict_plus[key11 + 1]
                                        key11 += 1

# finds any RNA that is antisense to a coding region and is above a threshold                     
for key21 in dict_both:
        if key21 in dict_strand:
                if dict_strand[key21] == "-":                       
                        #antisense_minus
                        if dict_plus[key21] > threshold:
                                antisense_of_minus[key21] = dict_plus[key21]                       
                elif dict_strand[key21] == "+":
                        #antisense_plus
                        if dict_minus[key21] > threshold:
                                antisense_of_plus[key21] = dict_minus[key21]                      

# finds any RNA that is in an intergenic region and above the threshold                                                                                              
for key2 in dict_both:
        if key2 not in dict_strand:                                
                if key2 not in UTR3_plus and key2 not in UTR5_plus:
                        if dict_plus[key2] > threshold:
                                RNA_plus[key2] = dict_plus[key2]                  
                if key2 not in UTR3_minus and key2 not in UTR5_minus :
                        if dict_minus[key2] > threshold:
                                RNA_minus[key2] = dict_minus[key2]
                                      

# creates a copy of a dictionary that will be used below. This is done so that I can modify the real dictionary while looping through the copy
copy_UTR3_minus = deepcopy(UTR3_minus)
# finds regions betwen genes that are likely to belong to an operon
for key3 in copy_UTR3_minus:
        if key3 in UTR5_minus:
                operon_minus[key3] = UTR3_minus[key3]
                del UTR3_minus[key3]                
                del UTR5_minus[key3]
                
copy_UTR3_plus = deepcopy(UTR3_plus)               
for key4 in copy_UTR3_plus:
        if key4 in UTR5_plus:
                operon_plus[key4] = UTR3_plus[key4]
                del UTR3_plus[key4]                
                del UTR5_plus[key4]

# the loop below correctly classifies UTRs that have been misclassified as antisense RNA due to overlapping coding regions
cop1_antisense_minus = deepcopy(antisense_of_minus)
cop1_antisense_plus = deepcopy(antisense_of_plus)
for pete in dict_strand:
        if dict_strand[pete] == "+":
                if pete + 1 in cop1_antisense_minus:
                        pos9 = pete + 1
                        while pos9 in cop1_antisense_minus:
                                UTR3_plus[pos9] = cop1_antisense_minus[pos9]
                                del antisense_of_minus[pos9]
                                pos9 += 1
                        cop1_antisense_minus = deepcopy(antisense_of_minus)
                elif pete - 1 in cop1_antisense_minus:
                        pos10 = pete - 1
                        while pos10 in cop1_antisense_minus:
                                UTR5_plus[pos10] = cop1_antisense_minus[pos10]
                                del antisense_of_minus[pos10]
                                pos10 -= 1
                        cop1_antisense_minus = deepcopy(antisense_of_minus)
        elif dict_strand[pete] == "-":
                if pete + 1 in cop1_antisense_plus:
                        pos11 = pete + 1
                        while pos11 in cop1_antisense_plus:
                                UTR5_minus[pos11] = cop1_antisense_plus[pos11]
                                del antisense_of_plus[pos11]
                                pos11 += 1
                        cop1_antisense_plus = deepcopy(antisense_of_plus)
                elif pete - 1 in cop1_antisense_plus:
                        pos12 = pete - 1
                        while pos12 in cop1_antisense_plus:
                                UTR3_minus[pos12] = cop1_antisense_plus[pos12]
                                del antisense_of_plus[pos12]
                                pos12 -= 1
                        cop1_antisense_plus = deepcopy(antisense_of_plus)

# joins the UTR region with antisense region if they are part of the same UTR 
temp_UTR5_minus = {}               
copy_antisense_of_plus = deepcopy(antisense_of_plus)
for key51 in UTR5_minus:
        if key51 + 1 in copy_antisense_of_plus:
                while key51 + 1 in copy_antisense_of_plus:
                        temp_UTR5_minus[key51 + 1] = antisense_of_plus[key51 + 1]
                        del antisense_of_plus[key51 + 1]
                        key51 += 1
                        
# merges two dictionaries together: will update the first dictionary with the second one, if any keys match
UTR5_minus.update(temp_UTR5_minus)

temp_UTR3_plus = {}
copy_antisense_of_minus = deepcopy(antisense_of_minus)
for key53 in UTR3_plus:
        if key53 + 1 in copy_antisense_of_minus:
                while key53 + 1 in copy_antisense_of_minus:
                        temp_UTR3_plus[key53 + 1] = antisense_of_minus[key53 + 1]
                        del antisense_of_minus[key53 + 1]
                        key53 += 1
                        
UTR3_plus.update(temp_UTR3_plus)

temp_UTR5_plus = {}
copy1_antisense_of_minus = deepcopy(antisense_of_minus)
for key54 in UTR5_plus:
        if key54 - 1 in copy1_antisense_of_minus:
                while key54 - 1 in copy1_antisense_of_minus:
                        temp_UTR5_plus[key54 - 1] = antisense_of_minus[key54 - 1]
                        del antisense_of_minus[key54 - 1]
                        key54 -= 1
                        
UTR5_plus.update(temp_UTR5_plus)

temp_UTR3_minus = {}
copy1_antisense_of_plus = deepcopy(antisense_of_plus)
for key52 in UTR3_minus:
        if key52 - 1 in copy1_antisense_of_plus:
                while key52 - 1 in copy1_antisense_of_plus:
                        temp_UTR3_minus[key52 - 1] = antisense_of_plus[key52 - 1]
                        del antisense_of_plus[key52 - 1]
                        key52 -= 1

UTR3_minus.update(temp_UTR3_minus)

# loops below merge any intergenic RNA with antisense RNA if they are part of the same transcript. Calls this new joined transcript either border_plus or border_minus ("mixed" ncRNAs)
border_minus = {}                           
for key5 in antisense_of_plus:
        if key5 - 1 in RNA_minus:
                border_minus[key5] = antisense_of_plus[key5]
                pos1 = key5
                while pos1 - 1 in RNA_minus:
                        border_minus[pos1 - 1] = RNA_minus[pos1 - 1]
                        pos1 -= 1
                pos2 = key5
                while pos2 + 1 in antisense_of_plus:
                        border_minus[pos2 + 1] = antisense_of_plus[pos2 + 1]
                        pos2 += 1
        elif key5 + 1 in RNA_minus:
                border_minus[key5] = antisense_of_plus[key5]
                pos3 = key5
                while pos3 + 1 in RNA_minus:
                        border_minus[pos3 + 1] = RNA_minus[pos3 + 1]
                        pos3 += 1
                pos4 = key5
                while pos4 - 1 in antisense_of_plus:
                        border_minus[pos4 - 1] = antisense_of_plus[pos4 - 1]
                        pos4 -= 1

border_plus = {}
for key7 in antisense_of_minus:
        if key7 - 1 in RNA_plus:
                border_plus[key7] = antisense_of_minus[key7]
                pos5 = key7
                while pos5 - 1 in RNA_plus:
                        border_plus[pos5 - 1] = RNA_plus[pos5 - 1]
                        pos5 -= 1
                pos6 = key7
                while pos6 + 1 in antisense_of_plus:
                        border_plus[pos6 + 1] = antisense_of_minus[pos6 + 1]
                        pos6 += 1
        elif key7 + 1 in RNA_plus:
                border_plus[key7] = antisense_of_minus[key7]
                pos7 = key7
                while pos7 + 1 in RNA_plus:
                        border_plus[pos7 + 1] = RNA_plus[pos7 + 1]
                        pos7 += 1
                pos8 = key7
                while pos8 - 1 in antisense_of_minus:
                        border_plus[pos8 - 1] = antisense_of_minus[pos8 - 1]
                        pos8 -= 1

# these loops just delete keys from antisense and RNA dictionaries that match the keys from border_plus or border_minus                        
for key9 in border_minus:
        if key9 in antisense_of_plus:
                del antisense_of_plus[key9]
        if key9 in RNA_minus:
                del RNA_minus[key9]

for key10 in border_plus:
        if key10 in antisense_of_minus:
                del antisense_of_minus[key10]
        if key10 in RNA_plus:
                del RNA_plus[key10]

# the following two loops joins mixed RNAs with antisense RNA if they are a part of the same transcript
border_plus_temp = {}
cop_antisense_minus = deepcopy(antisense_of_minus)
for a1 in border_plus:
        if a1 + 1 in cop_antisense_minus:
                apos1 = a1
                while apos1 + 1 in cop_antisense_minus:
                        border_plus_temp[apos1 + 1] = cop_antisense_minus[apos1 + 1]
                        del antisense_of_minus[apos1 + 1]
                        apos1 += 1
        elif a1 - 1 in cop_antisense_minus:
                apos2 = a1
                while apos2 - 1 in cop_antisense_minus:
                        border_plus_temp[apos2 - 1] = cop_antisense_minus[apos2 - 1]
                        del antisense_of_minus[apos2 - 1]
                        apos2 -= 1

border_plus.update(border_plus_temp)

border_minus_temp = {}
cop_antisense_plus = deepcopy(antisense_of_plus)
for a2 in border_minus:
        if a2 + 1 in cop_antisense_plus:
                apos3 = a2
                while apos3 + 1 in cop_antisense_plus:
                        border_minus_temp[apos3 + 1] = cop_antisense_plus[apos3 + 1]
                        del antisense_of_plus[apos3 + 1]
                        apos3 += 1
        elif a2 - 1 in cop_antisense_plus:
                apos4 = a2
                while apos4 - 1 in cop_antisense_plus:
                        border_minus_temp[apos4 - 1] = cop_antisense_plus[apos4 - 1]
                        del antisense_of_plus[apos4 - 1]
                        apos4 -= 1
                
border_minus.update(border_minus_temp)
                                
# here I call the functions to find ncRNAs that have been misclassified as operon or UTR
find_RNA_right(UTR5_minus)
find_RNA_right(UTR3_plus)
find_RNA_left(UTR5_plus)
find_RNA_left(UTR3_minus)
find_RNA_operon(operon_minus)
find_RNA_operon(operon_plus)
operon_to_UTR(operon_minus)
operon_to_UTR(operon_plus)

RNA_to_antisense_border(RNA_plus)
RNA_to_antisense_border(RNA_minus)

# here I call the function which filters dictionaries by length
filter_by_length(antisense_of_minus)
filter_by_length(antisense_of_plus)
filter_by_length(RNA_minus)
filter_by_length(RNA_plus)
filter_by_length(border_minus)
filter_by_length(border_plus)

proper_antisense_of_minus = {}
proper_antisense_of_plus = {}
proper_RNA_minus = {}
proper_RNA_plus = {}
proper_border_minus = {}
proper_border_plus = {}

find_peaks(antisense_of_minus, proper_antisense_of_minus)
find_peaks(antisense_of_plus, proper_antisense_of_plus)
find_peaks(RNA_minus, proper_RNA_minus)
find_peaks(RNA_plus, proper_RNA_plus)
find_peaks(border_minus, proper_border_minus)
find_peaks(border_plus, proper_border_plus)

os.system("mkdir " + folder + "/toRNAdo_output")

# here I call the function that writes a text file
write_to_file(antisense_of_minus, "antisense_of_minus") 
write_to_file(antisense_of_plus, "antisense_of_plus")
write_to_file(RNA_minus, "intergenic_minus")
write_to_file(RNA_plus, "intergenic_plus")
write_to_file(operon_minus, "operon_minus")
write_to_file(operon_plus, "operon_plus")
write_to_file(border_minus, "mixed_minus")
write_to_file(border_plus, "mixed_plus")
write_to_file(UTR3_plus, "UTR3_plus")
write_to_file(UTR5_plus, "UTR5_plus")
write_to_file(UTR3_minus, "UTR3_minus")
write_to_file(UTR5_minus, "UTR5_minus")
write_to_file(proper_antisense_of_minus, "filtered_antisense_of_minus")
write_to_file(proper_antisense_of_plus, "filtered_antisense_of_plus")
write_to_file(proper_border_minus, "filtered_mixed_minus")
write_to_file(proper_border_plus, "filtered_mixed_plus")
write_to_file(proper_RNA_minus, "filtered_intergenic_minus")
write_to_file(proper_RNA_plus, "filtered_intergenic_plus")

data_RNA_plus = open(folder + "/toRNAdo_output/filtered_intergenic_plus.wig").readlines()
data_RNA_minus = open(folder + "/toRNAdo_output/filtered_intergenic_minus.wig").readlines()
data_border_plus = open(folder + "/toRNAdo_output/filtered_mixed_plus.wig").readlines()
data_border_minus = open(folder + "/toRNAdo_output/filtered_mixed_minus.wig").readlines()
data_antisense_of_minus = open(folder + "/toRNAdo_output/filtered_antisense_of_minus.wig").readlines()
data_antisense_of_plus = open(folder + "/toRNAdo_output/filtered_antisense_of_plus.wig").readlines()

new_list = []

def sort_out_RNAs(openedfile, nameRNA, strand):
    temp_list = []
    temp_dict = {}
    temp_dict1 = {}
    for line in openedfile[1:]:
        temp_list.append(line.strip().split())
    temp_list = [[float(a) for a in b] for b in temp_list]
    for line1 in temp_list:
        temp_dict[line1[0]] = line1[1]
    for key in temp_dict:
        if key - 1 not in temp_dict:
            temp_dict1[key] = temp_dict[key]
            stop_pos = key
            while stop_pos + 1 in temp_dict:
                temp_dict1[stop_pos + 1] = temp_dict[stop_pos + 1]
                stop_pos += 1
            else:
                length = stop_pos - key + 1
                #value_array = numpy.array(temp_dict1.values())
                #mean_value = numpy.mean(value_array)
                peak_height_key = max(temp_dict1, key=lambda k: temp_dict1[k])
                peak_height = temp_dict1[peak_height_key]
                temp_row = []
                temp_row.extend([key, stop_pos, length, nameRNA, strand, peak_height])
                new_list.append(temp_row)
            temp_dict1.clear()
    temp_dict.clear()

new_list.insert(0, ["Start position", "Stop position", "Length", "Type of RNA", "Strand", "Peak height"])

sort_out_RNAs(data_antisense_of_minus, "antisense", "+")
sort_out_RNAs(data_antisense_of_plus, "antisense", "-")
sort_out_RNAs(data_border_plus, "mixed", "+")
sort_out_RNAs(data_border_minus, "mixed", "-")
sort_out_RNAs(data_RNA_plus, "intergenic", "+")
sort_out_RNAs(data_RNA_minus, "intergenic", "-")

with open(folder + "/toRNAdo_output/" + folder + "_ncRNAlist.csv", "wb") as savefile:
    writer = csv.writer(savefile)
    writer.writerows(new_list)