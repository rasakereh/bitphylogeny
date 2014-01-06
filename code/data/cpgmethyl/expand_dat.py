import os
import string

f = open('input.dat', 'r')

f.readline()

out = open('output_ZZ_IRX2P_5_L.dat','w')
out.write('"V1","V2","V3","V4","V5","V6","V7","V8"\n')
for line in f:
    (binstring,count) = line.split(" ")
    for i in range(int(count)):
        out.write(string.join(list(binstring),",")+'\n')
    
f.close()