# Kmer distribution analysis 

The *kmer_spectrum* script contains Kmer_spectrum class included methods 
* ```analyse``` which opens fastq file, creates dictionary Kmer where keys are kmers with lenght k and quality more than q (in Phred scale) and values are number of kmers in reads;

* ```spec``` which transforms Kmer dictionary to dataframe spectum, contained number of kmers and frequency (how many kmers with this number was found);

* ```spectrum_plot``` which creates a barplot;

* ```genome_size``` which calculates the approximate size of the genome; 

* ```border``` which calculates the border between true spectrum and sequencing noise.

File **test_kmer.fastq** was analysed with kmer length = 15 with and without quality filter. 

Parameter| with filter|without
-------- | -----------|-------
Coverage of k-mers| |112
Genome size | |2157945


Plot with filter


![a](https://github.com/PazhenkovaEA/bio_python/blob/master/hw3/filter.png)

Plot without filter


![b](https://github.com/PazhenkovaEA/bio_python/blob/master/hw3/no%20filter.png)

The red line shows the border between true spectrum and noise.
