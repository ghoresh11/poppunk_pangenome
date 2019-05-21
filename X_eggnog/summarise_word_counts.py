import os
import re
from itertools import combinations

def read_words_to_ignore():
    ignore = []
    with open("common_words.txt") as f:
        for line in f:
            toks = line.strip().split(",")
            ignore += toks
    return ignore

def read_eggnog_annotation(cluster, eggnog_file, word_counts, ignore):
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
            words_temp = toks[-1].lower()
            words_temp = re.sub('[^a-zA-Z0-9\- \n\.]', '', words_temp)
            words_temp = words_temp.split()
            words = []
            for w in words_temp: ## remove common words
                if w not in ignore:
                    words.append(w)
            for c in cogs_cat:
                if c not in word_counts:
                    word_counts[c] = {}

                ## take all word combinations
                for start, end in combinations(range(len(words)), 2):
                    phrase = " ".join(words[start:end+1])
                    if phrase not in word_counts[c]:
                        word_counts[c][phrase] = 0
                    word_counts[c][phrase] += 1

    return

def clean(word_counts):
    '''clean up the word word_counts
    alot of words are substrings of others, when the counts are the
    same only keep the longer one.
    Also, remove anything that only appears 10 times or less'''
    print("Cleaning duplicated phrases....")
    for c in word_counts:
        curr_words = word_counts[c].keys()
        for i in range(0, len(curr_words)-1):
            ## already been deleted
            if curr_words[i] not in word_counts[c]:
                continue
            ## ignore uncommon words
            if word_counts[c][curr_words[i]] < 10:
                del word_counts[c][curr_words[i]]
                continue
            for j in range(i+1, len(curr_words)):

                ## already been removed
                if curr_words[i] not in  word_counts[c]:
                    continue
                if curr_words[j] not in  word_counts[c]:
                    continue
                if word_counts[c][curr_words[j]] < 10:
                    del word_counts[c][curr_words[j]]
                    continue
                ## they have the same count almost exactly (+-10)
                if abs(word_counts[c][curr_words[i]] - word_counts[c][curr_words[j]]) < 10:
                    ## they're substrings of each other
                    if curr_words[i] in curr_words[j]:
                        del word_counts[c][curr_words[i]]
                    elif curr_words[j] in curr_words[i]:
                        del word_counts[c][curr_words[j]]
    return


def run():
    gene_types = ["rare", "inter", "core", "soft_core"]
    for gene_type in gene_types:
        print("Gene type: %s" %gene_type)
        word_counts = {}
        ignore = read_words_to_ignore()
        out = open(gene_type + "_words.csv", "w")
        out.write("COG_Cat, Word, Label, Count\n")
        for cluster in range(1,40):
            print("Reading cluster: %d" %cluster)
            read_eggnog_annotation(cluster, "eggnog/" + str(cluster) + "_" + gene_type + ".emapper.annotations", word_counts, ignore)

        clean(word_counts)
        for c in word_counts:
            for w in word_counts[c]:
                label = w
                if len(w) > 30:
                    label = w[:30]
                out.write(",".join(map(str,[c, w, label, word_counts[c][w]])) + "\n")

        out.close()
    return

if __name__ == "__main__":
    run()
