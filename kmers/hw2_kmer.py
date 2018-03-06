#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from Bio import SeqIO

class Kmer:
    counter = 1
    sequence = '' 

    def __init__(self,kmer_name,locus):
        self.sequence = kmer_name
        self.locus = [locus]
    
    def increase(self):
        self.counter += 1
        
    def increase_n(self, n):
        self.counter += n
        
    def add_locus(self, l):
        self.locus.append(l)
        
    def info(self):
        self.df = pd.DataFrame({'Kmer':[self.sequence],
                                'Locus':[self.locus],
                                'Amount': [self.counter]})
        print(self.df.to_string(index=False))

        
handle = open('./seq_y_pestis.fasta')
for record in SeqIO.parse(handle, "fasta"):
    seq = str(record.seq)
kmer_size = 23
seq_lng = len(seq)
kmer_dict = {}
index = 0
for index in range(seq_lng-kmer_size+1):
    current_kmer = seq[index:(index+kmer_size)]    
    if current_kmer in kmer_dict:
        kmer_dict[current_kmer].increase()
        kmer_dict[current_kmer].add_locus(index + 1)
    else:
        kmer_dict[current_kmer] = Kmer(current_kmer, index+1)

    
df = pd.DataFrame({'Kmer':[k for k in kmer_dict.keys()],
                   'Locus':[kmer_dict[key].locus for key in kmer_dict.keys()],
                   'Freq': [kmer_dict[key].counter for key in kmer_dict.keys()]})

    
df.loc[df['Freq'].idxmax()]  #Вывести инфу о самом распространенном к-мере
df.loc[df['Freq'].idxmax()]['Locus'] #Вывести только местоположения
#Если бластить самый распространенный к-мер, то можно определить, что со 100% идентичностью он совпадает только с Yersinia)
