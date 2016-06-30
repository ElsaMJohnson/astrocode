#!/usr/bin/python
import pyfits
import numpy
from math import *

#There are subtle differences between dr8 and dr7 where
#variables are stored.

pi=3.14159265358
rr = pi/180.0
rd = open('radec.txt',"r")
ts = open('tsfiles.txt', "r")
of = open('radeccoords','w')

for file, line in zip(ts, rd):
#file = 'tsField-004264-5-40-0261.fit'
  fl = pyfits.open(file)
  tbdata = fl[1].data
  nrows = 5
  node = fl[0].header['NODE']
  incl = fl[0].header['INCL']
  line=line.strip('\n')
  cols = line.split(' ')
  ra = float(cols[0])
  dec = float(cols[1])
  row = [float(0.0)]*5
  col = [float(0.0)]*5
  for j in range(0,5):
    a = tbdata.field('a')[0,j]
    b = tbdata.field('b')[0,j]
    c = tbdata.field('c')[0,j]
    d = tbdata.field('d')[0,j]
    e = tbdata.field('e')[0,j]
    f = tbdata.field('f')[0,j]
    r0 = tbdata.field('dRow0')[0,j]
    r1 = tbdata.field('dRow1')[0,j]
    r2 = tbdata.field('dRow2')[0,j]
    r3 = tbdata.field('dRow3')[0,j]
    c0 = tbdata.field('dCol0')[0,j]
    c1 = tbdata.field('dCol1')[0,j]
    c2 = tbdata.field('dCol2')[0,j]
    c3 = tbdata.field('dCol3')[0,j]
    cra=cos((ra-node)*rr)
    sra=sin((ra-node)*rr)
    sd=sin(dec*rr)
    cd=cos(dec*rr)
    ci=cos(incl*rr)
    si=sin(incl*rr)
    yy=(ci*sd-si*sra*cd)
    nu=(asin(yy))/rr
    x=cra*cd/cos(nu*rr)
    mu = atan2(sra*cd*ci+sd*si,cd*cra)/rr + node 
    row[j]=(f*(mu-a)-c*(nu-d))/(-c*e+b*f)
    col[j]=(e*(mu-a)-b*(nu-d))/(c*e-b*f)
  of.write(str(round(row[0],2))+" "+str(round(col[0],2))+" "+str(round(row[1],2))+" "+str(round(col[1],2))+" "+str(round(row[2],2))+" "+str(round(col[2],2))+" "+str(round(row[3],2))+" "+str(round(col[3],2))+" "+str(round(row[4],2))+" "+str(round(col[4],2))+'\n')

