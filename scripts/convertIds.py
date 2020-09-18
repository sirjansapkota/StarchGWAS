import os
import sys
import argparse

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description="This script was " + 
        "designed to convert between sorghum Sobic and Sb IDs. The " + 
        "sorghum ID is auto-recognized.")

    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument(
        "--input", 
        type=str, 
        required=True, 
        help="The name of the input file (single-column gene IDs).",
        action="store")

    parser.add_argument(
        "--output",
        type=str,
        required=False,
        help="The name of the output file.",
        action="store")

    parser.add_argument(
        "--synonyms",
        type=str,
        required=False,
        default="/zfs/tillers/Reference_Genomes/BTx623/v3.1.1/" \
                "annotation/Sbicolor_454_v3.1.1.synonym.txt",
        help="The synonym file containing key-value pairs.",
        action="store")

    parser.add_argument(
        "--alt_syn",
        type=str,
        required=False,
        default="/zfs/tillers/Reference_Genomes/BTx623/v3.1.1/" \
                "annotation/",
        help="Directory for '*annotation_info.txt' alternative " \
                "synonyms through A. thaliana.",
        action="store")

    return parser.parse_args()


def determine_orientation(genes):
    """Determine which direction to convert IDs"""
    if genes[0].startswith("Sobic"):
        orientation = "forward"
    elif genes[0].startswith("Sb"):
        orientation = "reverse"
    else:
        sys.stderr.write("Gene IDs not recognized\n")
        sys.exit(1)
    return orientation 


def generate_conversion_dictionary(orientation, synonyms):
    """Read in the IDs to convert between"""
    ids = {}

    with open(synonyms) as f:
        for line in f.read().splitlines():
            split_line = line.split()
            gene_id = ".".join(split_line[0].split('.')[:2])
            # Genes are Sobic
            if orientation == "forward":
                ids[gene_id] = split_line[1]
            # Genes are Sb
            elif orientation == "reverse":
                ids[split_line[1]] = gene_id
    return ids 


def annotations_to_dict(direction, file_name):
    """Read annotation file into dictionary"""
    annotation_dict = {}
    if "79" in file_name:
        sorghum_index = 0
        arab_index = 6
    elif "454" in file_name:
        sorghum_index = 1
        arab_index = 10
    with open(file_name) as f:
        for line in f.read().splitlines():
            split_line = line.split("\t")
            sorghum_gene = split_line[sorghum_index]
            arabidopsis_gene = split_line[arab_index]
            if direction == "to_sorghum":
                annotation_dict[arabidopsis_gene] = sorghum_gene
            elif direction == "to_ath":
                annotation_dict[sorghum_gene] = arabidopsis_gene
    return annotation_dict


def generate_alt_conversion_dictionaries(orientation, alt_syn_dir):
    """Read in the IDs to convert between"""
    # Sorghum [0], Arabidopsis [6]
    annotations = [i for i in os.listdir(alt_syn_dir) if "annotation_info" in i]

    # Set files
    # Genes are Sobic - 454_annotation
    if orientation == "forward":
        at_to_sb_file = alt_syn_dir + annotations[0]
        sb_to_at_file = alt_syn_dir + annotations[1]
    # Genes are Sb - 79_annotation
    elif orientation == "reverse":
        sb_to_at_file = alt_syn_dir + annotations[0]
        at_to_sb_file = alt_syn_dir + annotations[1]
    
    ath_to_sorghum = annotations_to_dict("to_sorghum", at_to_sb_file)
    sorghum_to_ath = annotations_to_dict("to_ath", sb_to_at_file)

    return ath_to_sorghum, sorghum_to_ath


def recover_missing_ids(gene, ath_to_sorghum, sorghum_to_ath):
    """For missing gene IDs, use Athaliana to recover"""
    try: 
        return ath_to_sorghum[sorghum_to_ath[gene]]
    except KeyError:
        try:
            gene_decimal = gene + ".1"
            return ath_to_sorghum[sorghum_to_ath[gene_decimal]]
        except KeyError:
            return "Unable to convert: {0}\n".format(gene)


def main(args):
    with open(args.input) as f:
        genes = f.read().splitlines()

    orientation = determine_orientation(genes)
    ids = generate_conversion_dictionary(orientation, args.synonyms)
    ath_to_sorghum, sorghum_to_ath = \
            generate_alt_conversion_dictionaries(orientation, args.alt_syn)
  
    if args.output:
        output = open(args.output, 'w')
        for gene in genes:
            try:
                output.write(ids[gene] + "\n")
            except KeyError:
                recovered = recover_missing_ids(gene, ath_to_sorghum, sorghum_to_ath)
                if "Unable" in recovered:
                    sys.stderr.write(recovered)
                else:
                    sys.stderr.write("Warning: {0} is ".format(recovered) + 
                            "recovered from Ath orthologs\n")
                    output.write(recovered + "\n")
        output.close()
    else:
        for gene in genes:
            try:
                print(ids[gene])
            except KeyError:
                recovered = recover_missing_ids(gene, ath_to_sorghum, sorghum_to_ath)
                if "Unable" in recovered:
                    sys.stderr.write(recovered)
                else:
                    sys.stderr.write("Warning: {0} is ".format(recovered) + 
                            "recovered from Ath orthologs\n")
                    print(recovered)
       

if __name__ == "__main__":
    args = parse_arguments()
    main(args)

