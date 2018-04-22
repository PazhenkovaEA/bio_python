# -*- coding: utf-8 -*-

from Bio import SeqIO
import Bio
from collections import defaultdict, Counter
import pygraphviz as pgv
import pandas as pd

class Vertex:
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}
        
    def increase_coverage(self):
        self.coverage += 1

class Edge:
    
    def __init__(self, k1, k2):
        self.seq = k1 + k2[-1]
        self.coverage = 0
    
    def calc_coverage(self, c1,c2):
        self.coverage = (c1+c2)/2

class Graph:

    def __init__(self, k):
        self.vertices = {}
        self.k = k
        self.new_coverage = 0
        self.passed = []
        self.index = []
        self.locmin = 0

    def add_read(self, read):
        read_lng = len(read)
        if read_lng < self.k:
            return
            
        kmer = read[:self.k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)
        
        for next_kmer_indx in range(1, read_lng-self.k+1, 1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx+self.k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)
            
            new_edge = Edge(kmer, next_kmer)
            
            self.vertices[next_kmer].in_edges[kmer] = [new_edge]
            
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]
            kmer = next_kmer
    
    def calc_init_edge_coverage(self):
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0].calc_coverage(self.vertices[current_vertex].coverage,self.vertices[next_vertex].coverage)

    def visualize(self, short, name):
        self.G = pgv.AGraph(directed=True)
        if short == 1:
            for current_vertex in list(self.vertices.keys()):
                for next_vertex in self.vertices[current_vertex].out_edges.keys():
                    self.G.add_node(current_vertex, label = "C=" + str(self.vertices[current_vertex].coverage))
                    self.G.add_node(list(self.vertices[current_vertex].out_edges.keys())[0], label="C=" + str(self.vertices[list(self.vertices[current_vertex].out_edges.keys())[0]].coverage))
                    self.G.add_edge(current_vertex, list(self.vertices[current_vertex].out_edges.keys())[0],label= "C=" + str(self.vertices[current_vertex].out_edges[next_vertex][0].coverage) +",L=" + str(len(self.vertices[current_vertex].out_edges[next_vertex][0].seq)-self.k + 1))
                for prev_vertex in self.vertices[current_vertex].in_edges.keys():
                    self.G.add_node(current_vertex, label="C=" + str(self.vertices[current_vertex].coverage))
                    self.G.add_node(list(self.vertices[current_vertex].in_edges.keys())[0], label="C=" + str(self.vertices[list(self.vertices[current_vertex].in_edges.keys())[0]].coverage))
                    self.G.add_edge(list(self.vertices[current_vertex].in_edges.keys())[0], current_vertex,label="C=" + str(self.vertices[current_vertex].in_edges[prev_vertex][0].coverage) + ",L=" + str(len(self.vertices[current_vertex].in_edges[prev_vertex][0].seq)-self.k + 1))
        else:
            for current_vertex in list(self.vertices.keys()):
                for nextvertex in self.vertices[current_vertex].out_edges.keys():
                    self.G.add_edge(current_vertex, list(self.vertices[current_vertex].out_edges.keys())[0],label=self.vertices[current_vertex].out_edges[nextvertex][0].seq)
                for prevvertex in self.vertices[current_vertex].in_edges.keys():
                    self.G.add_edge(list(self.vertices[current_vertex].in_edges.keys())[0], current_vertex, label=self.vertices[current_vertex].in_edges[prevvertex][0].seq)
        self.G.layout()
        self.G.draw(name, format='png', prog='dot')

    def compress(self):
        self.passed = [current_vertex for current_vertex in list(self.vertices.keys()) if (len(self.vertices[current_vertex].out_edges.keys()) == 1) & (len(self.vertices[current_vertex].in_edges.keys()) == 1)]
        for vertex in self.passed:
            self.prev = list(self.vertices[vertex].in_edges.keys())[0]
            self.next = list(self.vertices[vertex].out_edges.keys())[0]
            self.vertices[self.prev].out_edges[self.next] = [Edge(self.prev, self.next)]
            self.vertices[self.next].in_edges[self.prev] = [Edge(self.prev, self.next)]
            self.vertices[self.next].in_edges[self.prev][0].seq = self.vertices[vertex].in_edges[self.prev][0].seq + self.vertices[self.next].in_edges[vertex][0].seq[self.k:]
            self.vertices[self.prev].out_edges[self.next][0].seq = self.vertices[vertex].in_edges[self.prev][0].seq + self.vertices[self.next].in_edges[vertex][0].seq[self.k:]
            self.new_coverage = (self.vertices[vertex].out_edges[self.next][0].coverage * (len(self.vertices[vertex].out_edges[self.next][0].seq) - self.k + 1) + self.vertices[vertex].in_edges[self.prev][0].coverage * (len(self.vertices[vertex].in_edges[self.prev][0].seq) - self.k + 1) - self.vertices[vertex].coverage) / (len(self.vertices[self.prev].out_edges[self.next][0].seq) - self.k + 1)
            self.vertices[self.prev].out_edges[self.next][0].coverage = self.new_coverage
            self.vertices[self.next].in_edges[self.prev][0].coverage = self.new_coverage
            del self.vertices[self.prev].out_edges[vertex]
            del self.vertices[self.next].in_edges[vertex]
            del self.vertices[vertex]

    def cut_tails(self):
        #Поиск локального минимума - минимальное покрытие
        self.c = Counter()
        for current_vertex in list(self.vertices.keys()):
            self.c[self.vertices[current_vertex].coverage] += 1

        self.spectrum = pd.DataFrame({'Number': list(self.c.keys()),
                                 'Frequency': list(self.c.values())}).sort_values("Number").reset_index(drop=True)
        for n in range(10, len(self.spectrum["Frequency"]) - 10):
            if (self.spectrum["Frequency"][n - 10] > self.spectrum["Frequency"][n]) & (self.spectrum["Frequency"][n + 10] > self.spectrum["Frequency"][n]):
                self.locmin = n
            else:
                n += 1
        #найти хвосты
        self.index = [current_vertex for current_vertex in list(self.vertices.keys()) if (self.vertices[current_vertex].coverage <= self.locmin) & ((len(self.vertices[current_vertex].in_edges.keys()) == 1) & (len(self.vertices[current_vertex].out_edges.keys()) == 0) | (len(self.vertices[current_vertex].out_edges.keys()) == 1) & (len(self.vertices[current_vertex].in_edges.keys()) == 0))]

        for current_vertex in self.index:
            if len(self.vertices[current_vertex].out_edges.keys()) == 1:
                self.vertices[list(self.vertices[current_vertex].out_edges.keys())[0]].in_edges.pop(current_vertex)
            elif len(self.vertices[current_vertex].in_edges.keys()) == 1:
                self.vertices[list(self.vertices[current_vertex].in_edges.keys())[0]].out_edges.pop(current_vertex)
            self.vertices.pop(current_vertex)



#test

dataset = 'out1.fas'
k = 55
test = Graph(k)

with open(dataset, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        read = str(record.seq)
        test.add_read(read)
        test.add_read( str(record.reverse_complement().seq) )

test.calc_init_edge_coverage()
test.compress()
test.visualize(1, "compressed.png") #Сохраняет сжатый граф с сокращенным выводом в файл "compressed.png"
test.cut_tails()
test.visualize(1, "untailed_sh.png") #Сохраняет сжатый граф без хвостов с сокращенным выводом в файл
test.visualize(10, "untailed_full.png" )#Сохраняет сжатый граф без хвостов с полным выводом в файл
