#!/usr/bin/python
import pyfits
import numpy
from math import *

#There are subtle differences between dr8 and dr7 where
#variables are stored.

pi=3.14159265358
rr = pi/180.0

#for file in list
file = 'tsField-004264-5-40-0261.fit'
fl = pyfits.open(file)
#for file in lt:
tbdata = fl[1].data
nrows = 5
node = fl[0].header['NODE']
incl = fl[0].header['INCL']
#xycoords = [[10,10],[10,1470],[2030,1470],[2030,10]]
xycoords = [[898.0135, 1724.6251],[905.0513, 1724.1894],[895.4520,1724.835],[899.405793548, 1726.60651127],[904.079365805, 1724.75080734]]
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
#   xycoords = [[4,4],[4,1495],[2045,1496],[2045,4]]
 for i in range(len(xycoords)):
      x = xycoords[i][0]
      y = xycoords[i][1]
      xp = x+r0+r1*y+r2*y**2+r3*y**3
      yp = y+c0+c1*y+c2*y**2+c3*y**3
      mu = a + b*xp + c*yp
      nu = d + e*xp + f*yp
      t1 = sin((mu-node)*rr)*cos(nu*rr)*cos(incl*rr)
      t2 = sin(nu*rr)*sin(incl*rr)
      t3 = cos((mu - node)*rr)*cos(nu*rr)
      ra = (atan2((t1-t2),t3) + node*rr)/rr
      t4 = sin((mu-node)*rr)*cos(nu*rr)*sin(incl*rr)
      t5 = sin(nu*rr)*cos(incl*rr)
      dec = asin(t4+t5)/rr
#      if j==2: 
      print ra,dec

