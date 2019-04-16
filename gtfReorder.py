
import sys

infile = sys.argv[1]
outfile = infile[:-4] + '_reorder.gtf'

outlines = []

with open(infile) as a:
    for line in a.readlines():
        line = line.rstrip().split('\t')
        desc = line[-1].split('; ')
        if len(desc) == 3:
            newDesc = desc[1] +'; '+ desc[0] +'; '+ desc[2]
            newLine = '\t'.join(line[0:-1]) + '\t' + newDesc
            outlines.append(newLine)
        else:
            outlines.append('\t'.join(line))

with open(outfile,'w') as o:
    for line in outlines:
        print >> o, line
