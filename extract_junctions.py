#!/usr/bin/env python3

import sys
from collections import defaultdict

#long read alignments jitter so must include a bin size to adjsut
BIN_SIZE=10

#Defining functions
def strand_from_flag(flag):
    return "-" if (flag & 16) else "+"

def bin_pos(pos, bin_size=BIN_SIZE):
    return (pos // bin_size) * bin_size

def normalize_junction(a, b):
    return tuple(sorted([a,b]))

def parse_sa_tag(sa_tag):  
    if not sa_tag.startswith("SA:Z:"):
        return []

    sa_body = sa_tag.replace("SA:Z:", "").rstrip(";")
    if not sa_body:
        return []

    sa_entries = sa_body.split(";")
    alignments = []

    for entry in sa_entries:
        fields = entry.split(",")

        # SA format:
        # rname,pos,strand,CIGAR,mapQ,NM
        if len(fields) < 3:
            continue

        rname = fields[0]

        try:
            pos = int(fields[1])
        except ValueError:
            continue

        strand = fields[2]

        alignments.append((rname, pos, strand))

    return alignments

def extract_junctions(sam_file):
    junction_counts = defaultdict(set)

    with open(sam_file) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            
            fields = line.rstrip().split("\t")
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])

            if rname =="*" or pos == 0:
                continue

            strand = strand_from_flag(flag)
            primary = (rname, bin_pos(pos), strand)

            sa_tags = [f for f in fields[11:] if f.startswith("SA:Z:")]
            if not sa_tags:
                continue

            secondary_alignments = parse_sa_tag(sa_tags[0])

            for rname2, pos2, strand2 in secondary_alignments:
                secondary = (rname2, bin_pos(pos2), strand2)

                junction = normalize_junction(primary, secondary)
                junction_counts[junction].add(qname)

    return junction_counts

def main():
    if len(sys.argv) !=2:
        sys.stderr.write("usage: extract_junctions.py sample.sam\n")
        sys.exit(1)

    sam_file = sys.argv[1]
    junctions = extract_junctions(sam_file)

    print("chrA\tposA\tstrandA\tchrB\tposB\tstrandB\tread_count")
    for (a, b), reads in junctions.items():
        chrA, posA, strandA = a
        chrB, posB, strandB = b
        print(
            f"{chrA}\t{posA}\t{strandA}\t"
            f"{chrB}\t{posB}\t{strandB}\t"
            f"{len(reads)}"
        )

if __name__ == "__main__":
    main()
