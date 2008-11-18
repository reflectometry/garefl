#!/usr/bin/env python

f = open("par.dat","r")
for line in f:
    if line.startswith('#'):
        idx = line.find('parameter range:')
        if idx > 0:
            ranges = [pair.split(',') for pair in line[idx+16:].split()]
            ranges = [(float(x),float(y)) for (x,y) in ranges]
            [lo,hi] = zip(*ranges)
        print line.rstrip()
    else:
        values = [float(x) for x in line.split()]
        gen = int(values[0])
        chisq = float(values[1])
        pars = [float(x)*(b-a)+a for x,a,b in zip(values[2:],lo,hi)]
        print "  %15d %15g"%(gen,chisq)," ".join(["%15g"%(x) for x in pars])
