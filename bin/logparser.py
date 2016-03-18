# Parse script for extracting data from log file of GECO solver. 


# Ellery Ames, 2015-09-21

import sys, math

# infile and outfile
ifilename = sys.argv[1]
ofilename = sys.argv[2]

ifile = open(ifilename, 'r')
ofile = open(ofilename, 'w')

# Data headers
ofile.write('E0, Eb, Zc\n')

# Write data as comma separated table.
for line in ifile:
    if "E0" in line:
        splitline = line.split()
        e0val = splitline[3]
    elif "Eb" in line:
        sline = line.split()
        ebval = sline[2]
    elif "Zc" in line:
        sline = line.split()
        zcval = sline[2]    
        ofile.write('%s, %s, %s\n' % (e0val,ebval,zcval))
        
ifile.close()
ofile.close()    
