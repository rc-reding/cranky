#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import random
import sys
from pathlib import Path

from Bio import SeqIO

AMBIG = {
    "A"    :"A", "C"    :"C", "G"    :"G", "N"    :"N", "T"    :"T",
    "*A"   :"a", "*C"   :"c", "*G"   :"g", "*N"   :"n", "*T"   :"t",
    "AC"   :"M", "AG"   :"R", "AN"   :"a", "AT"   :"W", "CG"   :"S",
    "CN"   :"c", "CT"   :"Y", "GN"   :"g", "GT"   :"K", "NT"   :"t",
    "*AC"  :"m", "*AG"  :"r", "*AN"  :"a", "*AT"  :"w", "*CG"  :"s",
    "*CN"  :"c", "*CT"  :"y", "*GN"  :"g", "*GT"  :"k", "*NT"  :"t",
    "ACG"  :"V", "ACN"  :"m", "ACT"  :"H", "AGN"  :"r", "AGT"  :"D",
    "ANT"  :"w", "CGN"  :"s", "CGT"  :"B", "CNT"  :"y", "GNT"  :"k",
    "*ACG" :"v", "*ACN" :"m", "*ACT" :"h", "*AGN" :"r", "*AGT" :"d",
    "*ANT" :"w", "*CGN" :"s", "*CGT" :"b", "*CNT" :"y", "*GNT" :"k",
    "ACGN" :"v", "ACGT" :"N", "ACNT" :"h", "AGNT" :"d", "CGNT" :"b",
    "*ACGN":"v", "*ACGT":"N", "*ACNT":"h", "*AGNT":"d", "*CGNT":"b",
    "*"    :"-", "*ACGNT":"N",
}

AMBIG_2 = {
        "A": "A", "a": "a", "T": "T", "t": "t", "G": "G", "g": "g",
        "C": "C", "c": "c", "N": "N", "n": "n"
}


# Expects a multivcf file
# In the multivcf, if a sample does not specify a difference it is assumed to be the same as reference
# Will produce MSA with reference as first sequence. Only uses SNPs from the vcf

def extract_sample_names(vcf_file):
    """
    Extract sample names from VCF file
    """
    if vcf_file.lower().endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    sample_names = []
    with opener(vcf_file, "rt") as vcf:
        for line in vcf:
            line = line.strip("\n")
            if line.startswith("#CHROM"):
                record = line.split("\t")
                sample_names = [record[i].replace("./", "") for i in range(9, len(record))]
                break
    return sample_names

def is_anomalous(record, num_samples):
    """
    Determine if the number of samples in current record corresponds to number of samples described
    in the line '#CHROM'
    """
    return bool(len(record) != num_samples + 9)


def is_snp(record):
    """
    Determine if current VCF record is a SNP (single nucleotide polymorphism) as opposed to MNP
    (multinucleotide polymorphism)
    """
    # <NON_REF> must be replaced by the REF in the ALT field for GVCFs from GATK
    alt = record[4].replace("<NON_REF>", record[3])
    return bool(len(record[3]) == 1 and len(alt) - alt.count(",") == alt.count(",") + 1)

def get_matrix_column(record, num_samples, resolve_IUPAC):
    """
    Transform a VCF record into a phylogenetic matrix column with nucleotides instead of numbers
    """
    ref_base = record[3].replace("-","*").upper()
    nt_dict = {str(0): ref_base, ".": ref_base}
    # <NON_REF> must be replaced by the REF in the ALT field for GVCFs from GATK
    alt = record[4].replace("-", "*").replace("<NON_REF>", nt_dict["0"])
    alt = alt.split(",")
    for n in range(len(alt)):
        nt_dict[str(n+1)] = alt[n]
    column = "" + ref_base
    for i in range(9, num_samples + 9):
        geno_num = record[i].split(":")[0].replace("/", "").replace("|", "")
        try:
            geno_nuc = "".join(sorted(set([nt_dict[j] for j in geno_num])))
        except KeyError:
            return "malformed"
        if resolve_IUPAC is False:
            if len(geno_nuc) > 1:
                column += str('N')
            else:
                column += AMBIG_2[geno_nuc]
        else:
            column += AMBIG[nt_dict[random.choice(geno_num)]]
    return column


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input",
        action = "store",
        dest = "filename",
        required = True,
        help = "Name of the input VCF file, can be gzipped")
    parser.add_argument("-f", "--fasta",
        action = "store",
        dest = "ref",
        help = "Reference fasta. If given than a full msa will be produced.")
    parser.add_argument("--output-folder",
        action = "store",
        dest = "folder",
        default = "./",
        help = "Output folder name, it will be created if it does not exist (same folder as input by "
               "default)")
    parser.add_argument("--output-prefix",
        action = "store",
        dest = "prefix",
        help = "Prefix for output filenames (same as the input VCF filename without the extension by "
               "default)")
    parser.add_argument("-r", "--resolve-IUPAC",
        action = "store_true",
        dest = "resolve_IUPAC",
        help = "Randomly resolve heterozygous genotypes to avoid IUPAC ambiguities in the matrices "
               "(disabled by default)")
    args = parser.parse_args()


    # Get samples names and number of samples in VCF
    if Path(args.filename).exists():
        sample_names = extract_sample_names(args.filename)
    else:
        print("\nInput VCF file not found, please verify the provided path")
        sys.exit()
    num_samples = len(sample_names)
    if num_samples == 0:
        print("\nSample names not found in VCF, your file may be corrupt or missing the header.\n")
        sys.exit()
    print("\nConverting file '{}':\n".format(args.filename))
    print("Number of samples in VCF: {:d}".format(num_samples))

    if args.ref != '':
        if Path(args.ref).exists():
            ref = [record for record in SeqIO.parse(args.ref, "fasta")]
            if len(ref) > 1:
                print("\nInput reference should be a single fasta")
                sys.exit()
            ref = ref[0].seq
            print(f"Reference has {len(ref)} bases")
        else:
            print("\nRef fasta file not found, please verify the provided path")
            sys.exit()

    # Output filename will be the same as input file, indicating the minimum of samples specified
    if not args.prefix:
        parts = Path(args.filename).name.split(".")
        args.prefix = []
        for p in parts:
            if p.lower() == "vcf":
                break
            else:
                args.prefix.append(p)
        args.prefix = ".".join(args.prefix)

    # Check if outfolder exists, create it if it doesn't
    if not Path(args.folder).exists():
        Path(args.folder).mkdir(parents=True)

    outfile = str(Path(args.folder, args.prefix))

    # We need to create an intermediate file to hold the sequence data vertically and then transpose
    # it to create the matrices
    temporal = open(outfile+".tmp", "w")


    ##########################
    # PROCESS GENOTYPES IN VCF

    if args.filename.lower().endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    with opener(args.filename, "rt") as vcf:
        # Initialize line counter
        snp_num = 0
        snp_accepted = 0
        snp_shallow = 0
        mnp_num = 0
        snp_biallelic = 0
        position = 1 # positions are index 1 based

        while 1:
            # Load large chunks of file into memory
            vcf_chunk = vcf.readlines(50000)
            if not vcf_chunk:
                break

            for line in vcf_chunk:
                line = line.strip()

                if not line or line.startswith("#"): # skip empty and commented lines
                    continue

                # Split line into columns
                record = line.split("\t")

                # Fill in with reference for inbetween region
                record_pos = int(record[1])
                while position < record_pos:
                    temporal.write(ref[position-1] * (num_samples + 1) +"\n")
                    position += 1

                # Keep track of number of genotypes processed
                snp_num += 1
                # Print progress every 500000 lines
                if snp_num % 500000 == 0:
                    print("{:d} genotypes processed.".format(snp_num))
                if is_anomalous(record, num_samples):
                    print("Skipping malformed line:\n{}".format(line))
                    continue
                else:
                    # Check that neither REF nor ALT contain MNPs
                    if is_snp(record):
                        # Uncomment for debugging
                        # print(record)
                        # Transform VCF record into an alignment column
                        site_tmp = get_matrix_column(record, num_samples,
                                                        args.resolve_IUPAC)
                        # Uncomment for debugging
                        # print(site_tmp)
                        # Write entire row of single nucleotide genotypes to temp file
                        if site_tmp == "malformed":
                            print("Skipping malformed line:\n{}".format(line))
                            continue
                        else:
                            # Add to running sum of accepted SNPs
                            snp_accepted += 1
                            temporal.write(site_tmp+"\n")
                            position += 1
                    else:
                        # Keep track of loci rejected due to multinucleotide genotypes
                        mnp_num += 1
        
        # Add remaining bases at the end
        while position < len(ref):
            temporal.write(ref[position-1] * (num_samples + 1) +"\n")
            position += 1

        # Print useful information about filtering of SNPs
        print("Total of genotypes processed: {:d}".format(snp_num))
        print("Genotypes excluded because they exceeded the amount "
              "of missing data allowed: {:d}".format(snp_shallow))
        print("Genotypes that passed missing data filter but were "
              "excluded for being MNPs: {:d}".format(mnp_num))
        print("SNPs that passed the filters: {:d}".format(snp_accepted))

    print("")

    temporal.close()


    #######################
    # WRITE OUTPUT MATRICES


    output_fas = open(outfile+".fasta", "w") 

    # write reference first, which will be index 0
    with open(outfile+".tmp") as temporal:
        seqout = ""

        for line in temporal:
            seqout += line[0]

        output_fas.write(">reference" + "\n"+seqout+"\n")
        print(f"Reference Printed with {len(seqout)} bases.")

    # Write sequences of the ingroup
    for s in range(0, len(sample_names)):
        with open(outfile+".tmp") as temporal:
            seqout = ""

            for line in temporal:
                seqout += line[s+1]

            output_fas.write(">" + sample_names[s] + "\n"+seqout+"\n")
        print(f"Sample {sample_names[s]} Printed with {len(seqout)} bases")

    Path(outfile+".tmp").unlink()

    print( "\nDone!\n")

if __name__ == "__main__":
    main()
