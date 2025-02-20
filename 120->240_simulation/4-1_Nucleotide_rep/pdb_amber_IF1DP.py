import numpy as np
import os

os.chdir('/Users/ryohei/desktop')

pdbfile = 'F1-IF1DP_syn120deg_stall_nucletide_rep.pdb'
filename ='F1-IF1DP_syn120deg_stall_4amber.pdb'

f = open(pdbfile,'r')
Head = []
Anum = []; Aname = []; Rnum = []; Rname = []; Cname = []; Coord = [];
Occu = []; Beta = []; Symb = [];
for line in f:
    if line[0:4]=="ATOM" or line[0:6]=="HETATM": #line[0:3]=="TER" or
        Head.append(line[0:6])
        Anum.append(int(line[6:11]))
        Aname.append(line[12:16])
        Rname.append(line[17:20])
        Cname.append(line[21:22])
        Rnum.append(int(line[22:26]))
        Coord.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        Occu.append(float(line[54:60]))
        Beta.append(float(line[60:66]))
        Symb.append(line[75:78])

na = len(Anum)
print(na)
out = open(filename,'w')
omit = 0
for i in range(na):
    if Rname[i]=='ATP' or Rname[i]=='ADP':
        if Aname[i][3]=="\'":
            Aname[i] = Aname[i][0:3] + "*"
    if Aname[i] == ' MG ':
        Aname[i] = 'MG  '
        Rname[i] = 'MG '
    if Rname[i]=='ATP':
        Rname[i] = 'atp'
    if Rname[i]=='ADP':
        Rname[i] = 'adp'
    if Rname[i]=='PO4':
        Rname[i] = 'h2p'
        if Aname[i]==' P  ':
            Aname[i] = ' P1 '

    out.write('%6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f         %3s  \n' % (Head[i],Anum[i],Aname[i],Rname[i],Cname[i],Rnum[i],Coord[i][0],Coord[i][1],Coord[i][2],Occu[i],Beta[i], Symb[i]))
   
    if Aname[i]==' OXT':
        out.write('TER\n')   
    if (i>=25890 and i <=26044) and Rnum[i] != Rnum[i+1]:
        out.write('TER\n')

out.write('TER\n')   
out.close()