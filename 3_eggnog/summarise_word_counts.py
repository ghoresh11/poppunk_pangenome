import os
import re

def read_eggnog_annotation(cluster, eggnog_file, word_counts):
    ''' read an eggnot annotation file and count all the words for each
    COG category'''
    with open(eggnog_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if len(toks) < 12:
                continue
            cogs_cat = toks[-2].split(",")
            words = toks[-1].lower()
            words = re.sub('[^a-zA-Z0-9 \n\.]', '', words)
            words = words.split()
            for c in cogs_cat:
                if c not in word_counts:
                    word_counts[c] = {}
                for w in words:
                    if w not in word_counts[c]:
                        word_counts[c][w] = 0
                    word_counts[c][w] += 1


    return



def run():
    gene_types = ["rare", "inter", "core", "soft_core"]
    for gene_type in gene_types:
        word_counts = {}
        out = open(gene_type + "_words.csv", "w")
        out.write("COG_Cat, Word, Count\n")
        for cluster in range(1,40):
            read_eggnog_annotation(cluster, "eggnog/" + str(cluster) + "_" + gene_type + ".emapper.annotations", word_counts)

        for c in word_counts:
            for w in word_counts[c]:
                if word_counts[c][w] < 10:
                    continue
                out.write(",".join(map(str,[c, w, word_counts[c][w]])) + "\n")

        out.close()
    return

if __name__ == "__main__":
    run()
