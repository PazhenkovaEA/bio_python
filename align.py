#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import numpy as np

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Lets align it!')

    parser.add_argument('-s1', '--sequence1', help='The first sequence to align', metavar='Str', default='ALL')
    parser.add_argument('-s2', '--sequence2', help='The second sequence to align', metavar='Str', default='ALL')
    parser.add_argument('-f', '--flag', help='Print distance matrix', action='store_true')
    args = parser.parse_args()
    seq1 = args.sequence1
    seq2 = args.sequence2
    matr = args.flag
    

n = len(seq1)
m = len(seq2)
Wmin = np.zeros( (n+1,m+1), dtype = int )
trace = np.zeros( (n+1,m+1), dtype = str )
def comparison(a,b):
        if (a == b):
            return 0
        else:
            return 1
for i in range(n+1):
        Wmin[i][0] = 1 * i
        trace[i][0] = 'U'
        for j in range (m+1):
            Wmin [0][j] = 1 * j
            trace[0][j] = 'L'
            Wmin [0][0] = 0
for i in range(1, (n+1)):
        for j in range(1, (m+1)):
            a1 = Wmin[i,j-1] + 1
            a2 = Wmin[i-1,j] + 1
            a3 = Wmin[i-1,j-1] + comparison(seq1[i-1],seq2[j-1])
            Wmin[i][j] = min( [a1,a2,a3] ) 
i = n
j = m
r1 = ''
r2 = ''
while j > 0 or i > 0 :
        diag = Wmin[i-1][j-1]
        left = Wmin[i][j-1]
        up = Wmin [i-1][j]
        c1 = Wmin[i][j]
        if c1 == diag + comparison(seq1[i-1],seq2[j-1]):
            trace[i][j] = 'D'
            r1 += seq1[i-1]
            r2 += seq2[j-1]
            i -= 1
            j -= 1
        elif c1 == up + 1:
            trace[i][j] = 'U'
            r1 += seq1[i-1]
            r2 +='-' 
            i -= 1
        elif c1 == left + 1:
            trace[i][j] = 'L'
            r1 += '-'
            r2 += seq2[j-1]
            j -=1
if matr:
    print(Wmin)
print(r1[::-1])
print(r2[::-1])