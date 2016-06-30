#!/usr/bin/python
import pyfits
import numpy
from math import *

#There are subtle differences between dr8 and dr7 where
#variables are stored.Node and Incl are in the table
# not in the header. -Note the header contains the the 
# words Node and incl, but these are WRONG! They are correct
# in the photoField file, but this requires another download.
# The frames file is a list of image files from dr8 that contain
# the astrometry info. We need only this header info and not the 
# image

# Be careful! The file order might put the frame files as giruz.
pi=3.14159265358
rr = pi/180.0
rd = open(radecfile,"r")
fr = open(framelistfile, "r")
of = open(radec_coordsfile,'w')
#file = 'frame-r-003366-2-0101.fits'
# expand this file out into 5 lines apiece:
tr = open('tempra.txt','w')
nrows=0
for lines in rd:
  tr.write(lines)
  tr.write(lines)
  tr.write(lines)
  tr.write(lines)
  tr.write(lines)
  nrows = nrows +1

# This is a test
tr.close()
rd.close()
tn = open('tempra.txt','r')
i=0
j=0
k=0
row = [[float(0.0)]*5 for m in range(nrows)]
col = [[float(0.0)]*5 for m in range(nrows)]
for file,line in zip(fr,tn): 
 fl = pyfits.open(file)
 tbdata = fl[3].data
 node = tbdata.field('node')[0]
 incl = tbdata.field('incl')[0]
 a = tbdata.field('a')[0]
 if a<0.0:
   a=a+360.0
 elif a>360.0:
   a=a-360.0
 b = tbdata.field('b')[0]
 c = tbdata.field('c')[0]
 d = tbdata.field('d')[0]
 if d<0.0:
   d=d+360.0
 elif d>360.0:
   d=d-360.0
 e = tbdata.field('e')[0]
 f = tbdata.field('f')[0]
 r0 = tbdata.field('dRow0')[0]
 r1 = tbdata.field('dRow1')[0]
 r2 = tbdata.field('dRow2')[0]
 r3 = tbdata.field('dRow3')[0]
 c0 = tbdata.field('dCol0')[0]
 c1 = tbdata.field('dCol1')[0]
 c2 = tbdata.field('dCol2')[0]
 c3 = tbdata.field('dCol3')[0]
 line=line.strip('\n')
 cols = line.split(' ')
 ra = float(cols[0])
 dec = float(cols[1])
 cra=cos((ra-node)*rr)
 sra=sin((ra-node)*rr)
 sd=sin(dec*rr)
 cd=cos(dec*rr)
 ci=cos(incl*rr)
 si=sin(incl*rr)
 yy=(ci*sd-si*sra*cd)
 nu=(asin(yy))/rr
 if nu<0.0:
   nu=nu+360.0
 elif nu>360.0:
   nu=nu-360.0
 x=cra*cd/cos(nu*rr)
 mu = atan2(sra*cd*ci+sd*si,cd*cra)/rr + node 
 if mu<0.0:
   mu=mu+360.0
 elif mu>360.0:
   mu=mu-360.0
 j=int(i%5)
 k=int(i/5)
 foo1=(f*(mu-a)-c*(nu-d))/(-c*e+b*f)
 foo2=(e*(mu-a)-b*(nu-d))/(c*e-b*f)
 row[k][j] = foo1
 col[k][j] = foo2
 print i,j,k,foo1,foo2
 i=i+1

for l in range(nrows):
  of.write(str(round(row[l][3],2))+" "+str(round(col[l][3],2))+" "+str(round(row[l][0],2))+" "+str(round(col[l][0],2))+" "+str(round(row[l][2],2))+" "+str(round(col[l][2],2))+" "+str(round(row[l][1],2))+" "+str(round(col[l][1],2))+" "+str(round(row[l][4],2))+" "+str(round(col[l][4],2))+'\n')

of.close()
