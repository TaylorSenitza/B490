import numpy
import matplotlib.pyplot as plt

def parse_gff3():
    posn = list()
    fh = open('s_meliloti.gff3', 'r')
    for line in fh:
        if '\tgene\t' in line:
            values = line.split('\t')
            name = values[0]
            start = int(values[3])
            end = int(values[4])
            posn.append((name, start, end))
    return posn


def parse_fasta(data):
    """
    Read DNA, RNA, or protein sequences in Fasta format.

    This generator function yields a tuple containing a defline and a sequence
    for each record in the Fasta data. Stolen shamelessly from
    http://stackoverflow.com/a/7655072/459780.
    """
    name, seq = None, []
    for line in data:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))

#main
positions = parse_gff3()
lengths = []
gc = []
for defline, sequence in parse_fasta(open('s_meliloti.fa', 'r')):
    for p in positions:
        if p[0] == defline[1:].split()[0]:
            length = p[2] - p[1] #Starting position - end position
            gc_value = ((float(sequence[p[1]:p[2]].count("G")) + float(sequence[p[1]:p[2]].count("C")))/float(length))
            lengths.append(length)
            gc.append(gc_value)

plt.scatter(lengths, gc)
plt.show()
