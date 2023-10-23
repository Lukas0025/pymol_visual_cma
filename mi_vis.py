###
# Visualize MI corelated protain parts in pymol
# @file mi_vis.py
# @autor Lukáš Plevač <xpleva07@vutbr.cz>
# @date 2023.10.21

import pymol

import csv
import numpy as np
import functools
import matplotlib.pyplot as plt
import os

####
#### Configuration begin 
####
FASTA_FILE    = "aligned_and_cated_74_180.fasta"  # file with aligned seq
# fisrt amino seq must by referece from PDB with header with PDB identifier only for example
# >5x56_B

TH_MI         = 0     # Treshhold for MI value for accept
TH_LEN        = 2     # Treshhold for min amino seq length
LEN_MAX       = 4     # max len of amino seq
SHOW_ID       = 3     # id of corelated seq pairs to show in pymol
SHOW_FIRST    = 10    # show first 10 corelated pairs in table

MAX_SAME      = 0.75   # how match same can be a aminoasid strands in percents

MI_F          = "mi"  # mi or mip

####
#### Configuration end
####

##
# Convert index of amino acid in MSA to PYMOL seq index (resi)
# @param start start of subpart in MSA
# @param size  size of subpart in MSA
# @param msa   seqvence in msa (refernce what is loaded in pymol)
# @return start and size in resi
#
def msa2protIdx(start, size, msa):

    # 0123456789
    # --AB---AAA
    #   AB   AAA
    #   12   345

    for char in msa[start:start + size]:
        if char == "-":
            size -= 1

    for char in msa[0:start]:
        if char == "-":
            start -= 1


    return start, size



# create tmp R program for find MI
f = open("tmp_rcode.R", "a")
f.write("""
# get Bios2cor lib
if(!require('Bios2cor')) {
  install.packages('Bios2cor')
  library('Bios2cor')
}

#Importing MSA file load MI
align <- import.fasta(\"""" + FASTA_FILE + """\")

#Creating correlation object with MI method for positions with gap ratio < 0.2 (Default)
mi <- """ + MI_F + """(align)

write.table(mi$score, file="tmp_mi.csv", row.names=FALSE, col.names=FALSE, sep=";")
""")
f.close()

# first run R for mi on file
os.system('Rscript tmp_rcode.R')

# Load corelation matrix from file csv
mi = []
with open('tmp_mi.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    for row in spamreader:
        mi.append([float(col) for col in row])

# array for fonded corelated seqvences
founded_parts = []

#
# get reference amino acid from file
#
fasta = open(FASTA_FILE, 'r')
Lines = fasta.readlines()

if (Lines[0][0] != ">"):
    print("Error not correct fasta file")
    exit(1)

protain = Lines[0].split(">")[1].split("_")[0]
chain   = Lines[0].split(">")[1].split("_")[1][0]



# read seqvence of reference protain
base    = ""
for i in range(1, len(Lines)):
    if (Lines[i].count(">") > 0):
        break

    base += Lines[i].replace("\n", "").replace("\r", "")

##
# Get square in corelation matrix
# @param row_amino row for square start
# @param col_amino col for square start
# @param seq       corelation matrix
# @param p_fail    error in % for stop square (how match worser can be element of this square that fisrt element)
# @return int square size and avg corelation
#
def get_square(row_amino, col_amino, seq, p_fail = 0.05):

    sq_size  = 1
    max_size = min(
        len(seq) - row_amino,
        len(seq) - col_amino,
        max(row_amino, col_amino) - min(row_amino, col_amino), # one seq to other
        LEN_MAX
    )
    aminos   = np.array([])

    for i in range(max_size):
        if (seq[row_amino + i][col_amino + i] > seq[row_amino][col_amino] * (1 + p_fail)):
            break

        aminos  = np.append(aminos, seq[row_amino + i][col_amino + i]) 
        sq_size = i + 1


    return sq_size, np.mean(aminos)

# loop over corelation matrix and find squares
for row_amino in range(len(mi)):
    for col_amino in range(row_amino): # loop only under main diagonal
        size, avg_score = get_square(row_amino, col_amino, mi)

        if col_amino != row_amino and avg_score <= TH_MI:
            _, sizeC = msa2protIdx(col_amino, size, base)
            _, sizeR = msa2protIdx(row_amino, size, base)

            if min(sizeC, sizeR) >= TH_LEN:
                founded_parts.append({
                    "col":  col_amino,
                    "row":  row_amino,
                    "size": size, 
                    "mi":   avg_score
                })

##
# Get how match same is to squares in corelation matrix
# @param x1 squre in cmatrix
# @param x2 squre in cmatrix
# @return hom match in % is two squares same
def same_part(x1, x2):
    # convert squares to array of numbers (position in amino acid strad)
    x1_a = []
    x2_a = []

    for i in range(x1["size"]):
        x1_a.append(x1["col"] + i)
        x1_a.append(x1["row"] + i)

    for i in range(x2["size"]):
        x2_a.append(x2["col"] + i)
        x2_a.append(x2["row"] + i)
    
    # compare with smaller
    smaller = x1_a
    bigger  = x2_a
    if len(x1_a) > len(x2_a):
        smaller = x2_a
        bigger  = x1_a

    # check how match same aminoacids it have
    same = 0
    for amino in smaller:
        if amino in bigger:
            same += 1

    return same / len(smaller)

# delete duplicity
clean_parts = []
for x1 in founded_parts:
    clean = True
    for x2 in founded_parts:
        if  not(
                (# same end
                    (x1["col"] + x1["size"] == x2["col"] + x2["size"]) and
                    (x1["row"] + x1["size"] == x2["row"] + x2["size"])
                )
                # same start
                or ((x1["col"] == x2["col"]) and (x1["row"] == x2["row"]))
                # have same part
                or (same_part(x1, x2) > MAX_SAME)
            ) or (x1["size"] > x2["size"]) or x1 == x2:
            continue

        clean = False

    if clean:
        clean_parts.append(x1)


founded_parts = clean_parts



# Plot coreleation matrix with fonded corelated seqvences
plt.matshow(mi)

for sq in founded_parts:
    plt.plot(
        [sq["col"], sq["col"] + sq["size"] - 1],
        [sq["row"], sq["row"] + sq["size"] - 1]
    )
plt.show()

#
# Now work with pymol
#


##
# Clear pymol
#
def clear():
    pymol.cmd.reinitialize()
    
    # fetch protain
    pymol.cmd.fetch(protain)

    # color whole protain to gray
    pymol.cmd.color("gray", "(all)")

    pymol.cmd.hide("cartoon",  "(all)")
    #pymol.cmd.show("surface",  "(all)")

    # show protain as sticks
    pymol.cmd.show("sticks",  "(all)")


##
# Measure destance between coraled parts in protain 3D strucure
# @param part square in corelation matrix
# @return distance mesured in pymol
def mesure(part):
    
    seqAS, seqAL = msa2protIdx(part['col'], part['size'], base)
    seqBS, seqBL = msa2protIdx(part['row'], part['size'], base)

    print(f"selecting seq range:  {seqAS + 1}-{seqAS + seqAL}")
    print(f"selecting seq range:  {seqBS + 1}-{seqBS + seqBL}")

    pymol.cmd.color("yellow", f"(resi {seqAS + 1}-{seqAS + seqAL} and chain {chain})")
    pymol.cmd.color("red",    f"(resi {seqBS + 1}-{seqBS + seqBL} and chain {chain})")
    dst = pymol.cmd.distance( f"(resi {seqAS + 1}-{seqAS + seqAL} and chain {chain})", 
                              f"(resi {seqBS + 1}-{seqBS + seqBL} and chain {chain})", mode=4)

    return dst, seqAL


##
# Coprate two seqvences by distance in seq
#
def compareDST(item1, item2):
    diff1 = np.abs(item1["col"] - item1["row"])
    diff2 = np.abs(item2["col"] - item2["row"])

    if diff1 < diff2:
        return -1
    elif diff1 > diff2:
        return 1
    else:
        return 0
    
##
# Coprate two seqvences by MI
#
def compareMI(item1, item2):
    if item1["mi"] < item2["mi"]:
        return -1
    elif item1["mi"] > item2["mi"]:
        return 1
    else:
        return 0

# sort by distance
# find seq with biggest distances
founded_parts = sorted(founded_parts,  key=functools.cmp_to_key(compareMI))

clear()

for i in range(min(SHOW_FIRST, len(founded_parts))):
    founded_parts[i]["dist"], founded_parts[i]["RealSize"] = mesure(founded_parts[i])
    pymol.cmd.png(f"MI_SEQ_{i}")
    clear()

if len(founded_parts) > SHOW_ID:
    mesure(founded_parts[SHOW_ID])

# Finaly print table
print("ID\tSTART 1\tSTART 2\tSIZE\tAVG MI\tDISTANCE")
for i in range(min(SHOW_FIRST, len(founded_parts))):
    print(f"{i}\t{founded_parts[i]['col']}\t{founded_parts[i]['row']}\t{founded_parts[i]['RealSize']}\t{founded_parts[i]['mi']}\t{founded_parts[i]['dist']}")

# delete TMP files
os.remove("tmp_rcode.R")
os.remove("tmp_mi.csv")