#!/usr/bin/python
import pyfits
import numpy
from math import *

#There are subtle differences between dr8 and dr7 where
#variables are stored. For dr8 you have only the appropriate
#constants (a-f drow-dcol) stored in 3rd fits extension of each
#frame so there isn't a loop.

pi=3.14159265358
rr = pi/180.0

#for file in list
#file = 'tsField-004264-5-40-0261.fit'
file = 'frame-r-003366-2-0101.fits'
fl = pyfits.open(file)
#for file in lt:
tbdata = fl[3].data
xycoords = [[10,10],[10,1470],[2030,1470],[2030,10]]
node = tbdata.field('node')[0]
incl = tbdata.field('incl')[0]
a = tbdata.field('a')[0]
b = tbdata.field('b')[0]
c = tbdata.field('c')[0]
d = tbdata.field('d')[0]
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
      dec = (asin(t4+t5))/rr
      if ra>360:
	ra = ra-360.0
      elif ra<0.0:
	ra = ra+360.0
      print ra,dec

