import sys
import csv

infile = sys.argv[1]
outfile = infile + "_reduced.txt"

csv.field_size_limit(100000000)

with open(infile, "r") as source:
    rdr = csv.reader(source, delimiter='\t')
    with open(outfile, "w") as result:
        wtr = csv.writer(result, delimiter='\t')
        next(rdr)
        next(rdr)
        for r in rdr:
            wtr.writerow((r[0], r[6]))
