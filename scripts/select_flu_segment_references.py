import sys
import os
import argparse as ap

from Bio import SeqIO
import pandas as pd

fluAHAs = ["H1", "H2", "H3", "H4", "H5", "H6", "H7",
           "H8", "H9", "H10", "H11", "H12", "H13", "H14",
           "H15", "H16", "H17", "H18"]
fluANAs = ["N1", "N2", "N3", "N4", "N5", "N6", "N7",
           "N8", "N9", "N10", "N11"]

fluBHAs = ["Victoria", "Yamagata"]

def main():

    p = ap.ArgumentParser()

    p.add_argument("--alignmentStats", type=str, required=True, \
                   help="[Required] - A tab delimited file created using samtools idxstats containing reference genomes and the number of reads aligned to them.", \
                   action="store", dest="alignStats")
    p.add_argument("--referenceDB", type=str, required=True, \
                   help="[Required] - The database FASTA file containing sequences that reads were aligned to.", \
                    action='store', dest='refDB')
    p.add_argument("--outPref", type=str, required=True, \
                   help="[Required] - A prefix to use when naming output files.", \
                    action='store', dest="outPref")
    p.add_argument("--minReads", type=int, required=False, \
                   help='By default, the pipeline will only choose 1 reference for each segment to include in output. If supplied, the pipeline will include all sequences with greater than this number of reads aligned in the output, even if there are duplicate segments.', \
                    action='store', dest="minReads")
    
    args = p.parse_args()

    alignStats = ""
    if not os.path.exists(args.alignStats):
        sys.exit("ERROR: Alignment Stats File {0} does not exist. Please supply an existing file.".format(args.alignStats))
    else:
        alignStats = open(args.alignStats, "r")

    if not os.path.exists(args.refDB):
        sys.exit("ERROR: Reference Database File {0} does not exist. Please supply an existing file.".format(args.refDB))

    alignedDF = pd.DataFrame(columns= ["Virus", "Gene", "Classification", "Reads", "Reference"])
    for line in alignStats:
        split = line.strip("\n").split("\t")

        ref = line.strip("\n").split("\t")[0]
        readsAligned = int(line.strip("\n").split("\t")[1])

        if readsAligned > 0:

            gene = ref.split("_")[0]
            classification = ref.split("_")[1]

            if classification not in ["InfluenzaA", "InfluenzaB"]:
                if gene == "HA":
                    if gene in fluAHAs:
                        organism = "InfluenzaA"
                    elif gene in fluBHAs:
                        organism = "InfluenzaB"
                elif gene == "NA":
                    organism = "InfluenzaA"

            else:
                organism = classification

            alignedDF.loc[len(alignedDF)] = [organism, gene, classification, readsAligned, ref]

    print(alignedDF)

    alignedDF.to_csv("{0}-reference-alignment-stats.tab".format(args.outPref), sep="\t")

    referencesToUse = []
    if args.minReads:
        referencesToUse = alignedDF[alignedDF["Reads"] > args.minReads]["Reference"].to_list()
    else:
        referencesToUse = alignedDF.loc[alignedDF.groupby(["Virus", "Gene"])["Reads"].idxmax()]["Reference"].to_list()

    with open(args.refDB, "r") as dbHandle, open("{0}-references.fasta".format(args.outPref), "w+") as outHandle:
        for rec in SeqIO.parse(dbHandle, "fasta"):
            if rec.id in referencesToUse:
                SeqIO.write(rec, outHandle, "fasta")

if __name__ == "__main__":
    main()

