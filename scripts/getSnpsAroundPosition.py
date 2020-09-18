import sys
import argparse

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description="This script was " + 
        "designed to find any genes near SNPs found to be significant " + 
        "in a GWAS analysis performed in GAPIT.\n\n")

    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument(
        "--vcf", 
        type=str, 
        required=True, 
        help="VCF file.",
        action="store")

    requiredNamed.add_argument(
        "--midpoint",
        type=str,
        required=True,
        help="Position to pull SNPs around. Ex. S01_12250000.",
        action="store")

    parser.add_argument(
        "--distance",
        type=int,
        required=False,
        default=50000,
        help="Distance upstream and downstream to pull. Default: 50000 -- thus, 100000bp total.",
        action="store")

    return parser.parse_args()


def get_snps(vcf):
    snps = []
    with open(vcf) as f:
        for line in f:
            if not line.startswith("#"):
                split_line = line.split()
                snps.append(split_line[2])
    return snps


def pull_snps(snps, midpoint, distance):
    """Pull SNPs around 'midpoint' that are 'distance' upstream and down"""
    chrom = midpoint.split("_")[0]
    position = int(midpoint.split("_")[1])
    
    for snp in snps:
        split_snp = snp.split("_")
        temp_chrom = split_snp[0]
        temp_position = int(split_snp[1])
        if temp_chrom == chrom:
            if (max(position, temp_position) - min(position, temp_position)) <= distance:
                print(snp)


if __name__ == "__main__":
    args = parse_arguments()
    snps = get_snps(args.vcf)
    pull_snps(snps, args.midpoint, args.distance)
