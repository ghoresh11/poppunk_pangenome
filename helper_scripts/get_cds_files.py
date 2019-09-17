import gffutils
import os
import argparse

def get_gene_sequences(out_file, gff_file, fasta_file):
    '''Convert a gff file with the appended FASTA to protein/all fasta file
    gff_file = input gff file
    out_files = output files to write'''
    out = open(out_file, "w")
    db = gffutils.create_db(gff_file, dbfn=gff_file + '.db', force=True)
    gene_ids = [f.id for f in db.features_of_type('CDS')]
    for member in gene_ids:
        out.write(">" + member + "\n" + db[member].sequence(fasta_file) + "\n")
    out.close()
    return

def run(args):
    curr_genome = args.fasta_file.split("/")[-1].split(".")[0]
    out_file = os.path.join(args.out_dir, curr_genome + ".cds")
    get_gene_sequences(out_file, args.gff_file, args.fasta_file)
    return

def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from the GFF files to have for a reference')
    # input options
    parser.add_argument('--fasta_file', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/genomes/enterobase/human/complete/esc_ga8906aa_as.fasta",
                        help='path to inpqut fasta file [%(default)s]')
    parser.add_argument('--gff_file', required=False,
                            type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/rerun_prokka/complete/esc_ga8906aa_as.gff",
                            help='path to input gff file [%(default)s]')
    parser.add_argument('--out_dir',
                        type=str, default = ".",
                        help='output path [%(default)s]')
    return parser.parse_args()

if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
