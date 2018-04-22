# Assembler, based on de Bruijn graph

The script contains the class "Graph"
with methods: 
* **add_read** to add sequences and to divide it on k-mers;
* **calc_init_edge_coverage** to calculate kmer coverage of nodes and edges of non-compressed graph;
* **compress** return the compressed graph (droppes vertices with indegree and outdegree = 1, coverage of edges is recalculated);
* **visualize** save the graph in png format. The method takes name of file and integer: 1 to build short graph (with coverages on nodes and coverages and lengths on edges) or any other to build the full graph with sequences on nodes and edges;
* **cut_tails** cut so-called tailes - vertices with low coverage and degree =1.

# Test
The script was tested on 1000 reads ([out1.fas](https://github.com/PazhenkovaEA/bio_python/blob/master/hw4_5/test/out1.fas)) with k=55.
1. After compression

![After compression](https://github.com/PazhenkovaEA/bio_python/blob/master/hw4_5/test/compressed.png)

2. After tail-cutting (full)

![tail-cutsh](https://github.com/PazhenkovaEA/bio_python/blob/master/hw4_5/test/untailed_full.png)

3. After tail-cutting (short)

![tail-cutsh](https://github.com/PazhenkovaEA/bio_python/blob/master/hw4_5/test/untailed_sh.png)



Then script was tested on reads obtained from genome of Hepatitis B virus (hw_4_5_dataset.fasta).
This is a graph after compression and tail-cutting.

![pain](https://github.com/PazhenkovaEA/bio_python/blob/master/hw4_5/test/hw4_5.png)

The assembled sequences was written in fasta with the following code and then blasted.
```python
contigs = []
for current_vertex in list(test.vertices.keys()):
    if len(test.vertices[current_vertex].out_edges.keys()) != 0:
        contigs.append(test.vertices[current_vertex].out_edges[list(test.vertices[current_vertex].out_edges.keys())[0]][0].seq)
    if len(test.vertices[current_vertex].in_edges.keys()) != 0:
        contigs.append(test.vertices[current_vertex].in_edges[list(test.vertices[current_vertex].in_edges.keys())[0]][0].seq)
con_1 = set(contigs)

fin = []
for seq1 in con_1:
    fin.append(SeqRecord(Seq(seq1), id = 'contig'))
output_handle = open("contigs.fasta", "w")
for seq in fin:
    SeqIO.write(seq, output_handle, "fasta"),
output_handle.close()
```
