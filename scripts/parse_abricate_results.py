import sys
import os

import argparse as ap
import pandas as pd

fluAHAs = ["H1", "H2", "H3", "H4", "H5", "H6", "H7",
           "H8", "H9", "H10", "H11", "H12", "H13", "H14",
           "H15", "H16", "H17", "H18"]
fluANAs = ["N1", "N2", "N3", "N4", "N5", "N6", "N7",
           "N8", "N9", "N10", "N11"]

fluBHAs = ["Victoria", "Yamagata"]

def main():

    p = ap.ArgumentParser()

    p.add_argument("-i", "--input", type=str, \
                   help="[Required] - A file containing the tab delimited output of abricate for parsing.", \
                   action='store', dest='inFile')
    
    args = p.parse_args()

    data = ''
    if not os.path.exists(args.inFile):
        sys.exit("ERROR: Input File {0} does not exist. Please supply an existing file.".format(args.inFile))
    else:
        data = pd.read_csv(args.inFile, header=0, sep="\t")

    #print(data)
    #print(list(data.columns))

    TypesDetected = data[data["GENE"] == "M1"]["RESISTANCE"].str.replace("_", " ").unique().tolist()
    HAsDetected = data[data["GENE"] == "HA"]["RESISTANCE"].unique().tolist()
    NAsDetected = data[data["GENE"] == "NA"]["RESISTANCE"].unique().tolist()



    for ha in HAsDetected:
        if (ha in fluAHAs) and ("Influenza A" not in TypesDetected):
            TypesDetected.append("Influenza A")
        elif (ha in fluBHAs) and ("Influenza B" not in TypesDetected):
            TypesDetected.append("Influenza B")
    
    for na in NAsDetected:
        if (na in fluANAs) and ("Influenza A" not in TypesDetected):
            TypesDetected.append("Influenza A")

    results = []
    if len(TypesDetected) == 0:
        results.append("Unable to Type Sample")
    else:
        for type in TypesDetected:

            subtype = ""

            if type == "Influenza A":
                HAs = [ha for ha in HAsDetected if ha in fluAHAs]
                NAs = [na for na in NAsDetected if na in fluANAs]

                if len(HAs) == 1 and len(NAs) == 1:
                    subtype = HAs[0] + NAs[0]
                else:
                    if len(HAs) == 0:
                        HASubtypes = "Unable to Subtype HA"
                    else:
                        HASubtypes = ", ".join(HAs)
                    
                    if len(NAs) == 0:
                        NASubtypes = "Unable to Subtype NA"
                    else:
                        NASubtypes = ", ".join(NAs)

                    subtype = "HA: {0} NA: {1}".format(HASubtypes, NASubtypes)
            if type == "Influenza B":
                HAs = [ha for ha in HAsDetected if ha in fluBHAs]

                if len(HAs) == 1:
                    subtype = HAs[0]
                elif len(HAs) == 0:
                    subtype = "Unable to Subtype"
                else:
                    subtype = ", ".join(HAs)

            results.append("{0} {1}".format(type, subtype))

    if len(results) == 1:
        print(results[0])
    else:
        print("Multiple Detected: [{0}]".format("; ".join(results)))

if __name__ == "__main__":
    main()