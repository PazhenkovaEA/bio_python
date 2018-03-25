from collections import defaultdict, Counter
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Kmer_spectrum:
    Kmer = defaultdict(int)
    c = Counter()
    def __init__(self, file, k, q):
        self.file = file
        self.k = k
        self.q = q
    def analyse(self):
        handle = open(self.file)
        self.rec = [record for record in SeqIO.parse(handle, "fastq")]
        handle.close()
        for key in self.rec:
            sequence = str(key.seq)
            for index in range(len(sequence) - self.k + 1):
                current_kmer = sequence[index:(index + self.k)]
                if all(t> self.q for t in key.letter_annotations["phred_quality"][index:(index + self.k)]):
                    self.Kmer[current_kmer] += 1
    def spec(self):
        for element in list(self.Kmer.values()):
            self.c[element] += 1
        self.spectrum = pd.DataFrame({'Number': list(self.c.keys()),
                                 'Frequency': list(self.c.values())}).sort_values("Number").reset_index(drop=True)
    def border(self): #граница отсечки шумовой части от нешумовой в первом локальном минимуме
        for n in range(10, len(self.spectrum["Frequency"]) - 10):
            if (self.spectrum["Frequency"][n - 10] > self.spectrum["Frequency"][n]) & (self.spectrum["Frequency"][n + 10] > \
                    self.spectrum["Frequency"][n]):
                self.locmin = n
                return self.locmin
            else:
                n += 1

    def spectrum_plot(self):
        self.xmax = len(self.spectrum['Number'])
        self.ymax = self.spectrum['Frequency'].max()
        plt.bar(np.arange(self.xmax), self.spectrum['Frequency'], color='darkgreen', width=1)
        plt.xlabel('Number of Kmers')
        plt.ylabel('Frequency')
        plt.title('Spectrum of Kmers')
        self.v = ([1, self.xmax, 1, self.ymax])
        plt.axis(self.v)
        plt.legend()
        plt.yscale('log')
        plt.axvline(x=self.border(), color='r')
        plt.show()
    def genome_size(self):
        self.maximum = self.spectrum['Frequency'][self.border():].idxmax()
        self.size = sum(map(lambda x, y:x*y, list(self.spectrum['Number'][self.border():]), list(self.spectrum['Frequency'][self.border():])))/self.spectrum['Number'][self.maximum]
        return(self.size)


a = Kmer_spectrum('test_kmer.fastq', 15, 20)
a.analyse()
a.spec()
a.spectrum_plot()
a.genome_size()


b = Kmer_spectrum('test_kmer.fastq', 15, 0)
b.analyse()
b.spec()
b.spectrum_plot()
b.genome_size()

