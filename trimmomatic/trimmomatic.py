from Bio import SeqIO
import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='My Trimmomatic')

    parser.add_argument('-c', '--headcrop', help='Cut c nucleotides from the begin', metavar='Int', type=int, default=0)
    parser.add_argument('-t', '--tailcrop', help='Cut t nucleotides from the end', metavar='Int', type=int, default=0)
    parser.add_argument('-q', '--quality', help='Quality threshold', metavar='Int', type=int, default=0)
    parser.add_argument('-s', '--slidewindow', help='Size of sliding window', metavar='Int', type=int, default=0)
    parser.add_argument('-i', '--input', help='Input filename (w/o extension)', metavar='FILE', required=True)
    parser.add_argument('-o', '--output', help='Results filename (w/o extension)', metavar='FILE', required=True)
    args = parser.parse_args()
    headcrop = args.headcrop
    tailcrop = args.tailcrop
    window = args.slidewindow
    treshold = args.quality
    inpname = args.input
    outname = args.output

handle = open(inpname + ".fastq")
rec = [record for record in SeqIO.parse(handle, "fastq")]
handle.close()
final = []
for key in rec:
    trimmed = key[headcrop:len(str(key.seq))-tailcrop]
    for k in range (len(str(trimmed.seq)) - window):
        if sum(trimmed.letter_annotations["phred_quality"][k:k+window])/window > treshold:
            k += 1
        else:
            trimmed = trimmed[:k]
            continue
    final.append(trimmed)

output_handle = open(outname + ".fastq", "w")
for seq in final:
    SeqIO.write(seq, output_handle, "fastq"),
output_handle.close()
