# Kmer distribution analysis 

The *kmer_spectrum* script contains Kmer_spectrum class included methods 
```analyse``` which opens fastq file, creates dictionary Kmer where keys are kmers with lenght k and quality more than q (in Phred scale) and values are number of kmers in reads;
```spec``` which transforms Kmer dictionary to dataframe spectum, contained number of kmers and frequency (how many kmers with this number was found);
```spectrum_plot``` which creates a barplot;
```genome_size``` which calculates the approximate size of the genome; 
```border``` which calculates the border between true spectrum and sequencing noise.
