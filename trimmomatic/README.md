# Mini trimmomatic

The trimming tool for NGS reads in fastq format.

## Running 

The script can be runned from terminal with the following command:
```
python3 trimmomatic.py [headcrop] [tailcrop] [quality] [slidewindow] [input] [output]

```
### Flags

* -c, --headcrop Number of the nucleotides trimmed from the 5' of the sequence
* -t, --tailcrop Number of the nucleotides trimmed from the end of sequence
* -q, --quality Quality threshold
* -s, --slidewindow Size of sliding window
* -i, --input Input filename (w/o extension)
* -o, --output Results filename (w/o extension)

## Test

612348 sequences, length 30-100bp

```

python3 trimmomatic.py -c 10 -t 10 -q 20 -s 10 -i test_classwork3 -o trimmed

```

**FastQC report before trimming**

Per base sequence quality


![before1](https://github.com/PazhenkovaEA/bio_python/blob/master/trimmomatic/fastqc/before1.png)

Per sequence quality score

![before2](https://github.com/PazhenkovaEA/bio_python/blob/master/trimmomatic/fastqc/before2.png)

**FastQC report after trimming**

Per base sequence quality

![after](https://github.com/PazhenkovaEA/bio_python/blob/master/trimmomatic/fastqc/after1.png)

Per sequence quality score

![after1](https://github.com/PazhenkovaEA/bio_python/blob/master/trimmomatic/fastqc/after.png)



