import random

done = []
shapes = ["DIAMOND", "ELLIPSE", "RECTANGLE", "OCTAGON", "HEXAGON", "TRIANGLE", "VEE"]

with open("tree_labels.txt") as f:
	for line in f:
		toks = line.strip().split(",")
		cluster = toks[1]
		col = toks[3]
		if cluster in done:
			continue
		done.append(cluster)
	#	print("<discreteMappingEntry value=\"" + col + "\" attributeValue=\"" + cluster + "\"/>")	
		shape = random.choice(shapes)
		print('<discreteMappingEntry value="' + shape + '" attributeValue="' + cluster + '"/>')