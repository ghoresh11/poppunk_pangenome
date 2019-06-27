import os

#

output_files = os.listdir(".")

print("Cluster1,Cluster2,Mean,SD,Min")

for f in files:
    toks = f.split("_")
    if len(toks) == 2:
        cluster1 = toks[0]
        cluster2 = toks[0]
    else:
        cluster1 = toks[0]
        cluster2 = toks[1]
    with open(f) as f_open:
        for line in f:
            if "Mean:" in line:
                mean = line.strip().split()[-1]
                continue
            if "SD:" in line:
                sd = line.strip().split()[-1]
                continue
            if "Min:" in line:
                min = line.strip().split()[-1]
                break
    print(",".join([cluster1, cluster2, mean, sd, min]))