        PROGRAM RD2PIX

c	Converts ra dec to row and column values for SDSS images
c	(Sloan digital survey)
c	remember to specify with -o name when compiling.

c	This is the complete coordinate transformation, you will need
c	the following text files:
c	list containing ra and dec
c	list containing corresponding node and incllination (From ts files)
c	list containing a-f parameters (from ts files)
c	list containing drow[0-3] and dcols[0-3] from ts files
c	in ds9, the x direction has the most pixels (~2000)
c	note: I the amperstand, &, its only place is the 6th column 


        real*8 a(4000),b(4000),c(4000),d(4000),e(4000),f(4000)
	real*8 node(4000),iang(4000),mu,nu,rowp(4000),colp(4000)
        real*8 x,yy,cra,sra,sd,cd,ci,si,ra(4000),dec(4000)
	real*8 r0(4000),r1(4000),r2(4000),r3(4000)
	real*8 c0(4000),c1(4000),c2(4000),c3(4000),x2
        real*8 row(4000),col(4000),cll(40000),dy(40000),yp(40000)
        real*8 urw(4000),grw(4000),rrrw(4000),irw(4000),zrw(4000)
	real*8 ucl(4000),gcl(4000),rcl(4000),icl(4000),zcl(4000)
	CHARACTER*28 rdfile,inclfile,drcfile,atoffile


c	alpha=95.0	
	pi=3.14159265358
        open(unit=8,file='rd2pix.txt', status='unknown')
 3      format(A20)
	write(*,*) 'how many sn?'
	read(*,*) num
        write(*,*) 'Input radec file'
        read(*,3) rdfile
        write(*,*) 'Input incl file'
        read(*,3) inclfile
        write(*,*) 'Input drow dcol file'
        read(*,3) drcfile
        write(*,*) 'Input atof file'
        read(*,3) atoffile
	
        open(unit=11,file=rdfile,form='formatted',status='old')
	open(unit=12,file=inclfile,form='formatted',status='old')
 	open(unit=13,file=atoffile,form='formatted',status='old')
 	open(unit=14,file=drcfile,form='formatted',status='old')

	do 10 i=1,num
          read(11,*) ra(i),dec(i)
	  write(99,*) ra(i),dec(i)
          read(12,*) node(i),iang(i)
          write(99,*) node(i),iang(i),'node, incl'
	  sra=sin((ra(i)-node(i))*pi/180.0)
	  cra=cos((ra(i)-node(i))*pi/180.0)
c......all of the following are positive.?
	  sd=sin(dec(i)*pi/180.0)
	  cd=cos(dec(i)*pi/180.0)
          ci=cos(iang(i)*pi/180.0)
          si=sin(iang(i)*pi/180.0)
c      write(*,*) sra, sd,cd,ci,si,' sra sd cd ci si'
	  yy=(ci*sd-si*sra*cd)
c          yy=(ci*sd-si*sra)
          write(*,*) 'yy', yy
	  nu=(asin(yy))*180.0/pi
 	 write(*,*) nu, 'nu'        
c	  x=cra*cd/cos(nu*pi/180.0)
c        x=cra/(cos(nu*pi/180.0))
	x=(sra*cd*ci+sd*si)
        x2=cra*cd
	write(*,*) x	

        mu=node(i)+atan2(x,x2)*180.0/pi
c          if (ra(i).lt.90.0) then
c		mu=acos(x)*180.0/pi+node(i)
c         elseif (ra(i).gt.275.0) then
c              mu=acos(-x)*180.0/pi+node(i)
c         else if (node(i).gt.180.0) then
c	      mu=-acos(x)*180.0/pi+node(i)
c         else	
c              mu=acos(x)*180.0/pi+node(i)
c         endif

	 write(99,*) ra(i),dec(i),mu,nu, ' ra,dec,mu nu'
	 write(99,*) dec(i),nu,d(i),' dec nu d'
	write(*,*) ra(i),mu,nu,x
c        if (mu.gt.360) mu=mu-360.0
        k=(i-1)*5
	do 20 m=1,5
             l=k+m
         read(13,*) a(l),b(l),c(l),d(l),e(l),f(l)
         read(14,*) r0(l),r1(l),r2(l),r3(l),c0(l),c1(l),c2(l),c3(l)
         write(99,*) a(l),b(l),c(l),d(l),e(l),f(l)
         write(99,*) r0(l),r1(l),r2(l),r3(l),c0(l),c1(l),c2(l),c3(l)
	
         rowp(l)=(f(l)*(mu-a(l))-c(l)*(nu-d(l)))/(-c(l)*e(l)+b(l)*f(l))
	 colp(l)=(e(l)*(mu-a(l))-b(l)*(nu-d(l)))/(c(l)*e(l)-b(l)*f(l))
	 poop=0.5
	 ind=0
	 sht=colp(l)
	 do 30 j=1,160
c             sht=colp(l)
             cll(j)=sht-8.0+1.0*j/10.0

c       A few changes, calculate y first by trial and then compute x
	 
	yp(j)=c0(l)+cll(j)+c1(l)*cll(j)+c2(l)*(cll(j))**2+c3(l)*(cll(j))**3
        dy(j)=abs(yp(j)-sht)
	if (dy(j).lt.poop) then
		poop=dy(j)
		ind=j
		col(l)=cll(j)
c	else  poo=j
	end if
  30     continue
		 write(*,*) poop,ind
	
	 row(l)=rowp(l)-(r0(l)+r1(l)*col(l)+r2(l)*(col(l))**2+r3(l)*(col(l))**3)
	 if (m.eq.1) then
	 	urw(i)=row(l)
		ucl(i)=col(l)
	 else if (m.eq.2) then
		grw(i)=row(l)
                gcl(i)=col(l)
         else if (m.eq.3) then
                rrrw(i)=row(l)
                rcl(i)=col(l)
         else if (m.eq.4) then
                irw(i)=row(l)
                icl(i)=col(l)
         elseif (m.eq.5) then
                zrw(i)=row(l)
                zcl(i)=col(l)
 	 end if
  20	      continue
 104       format(2f16.3,10(1x,f11.3))
c      write(8,104) ra(i),dec(i),ucl(i),urw(i),gcl(i),
c     & grw(i),rcl(i),rrrw(i),icl(i),irw(i),zcl(i),zrw(i)  
c      write(*,*) ra(i),dec(i),ucl(i),urw(i),gcl(i),
c     & grw(i),rcl(i),rrw(i),icl(i),irw(i),zcl(i),zrw(i) 	

c      write(8,104) ra(i),dec(i),urw(i),grw(i),rrrw(i),irw(i),
c     & zrw(i),ucl(i),gcl(i),rcl(i),icl(i),zcl(i)
      write(8,104) ra(i),dec(i),gcl(i),grw(i),icl(i),
     & irw(i),rcl(i),rrrw(i),ucl(i),urw(i),zcl(i),zrw(i)
c       write(8,*) ra(i),dec(i),gcl(i),grw(i),icl(i),irw(i),
c     & rcl(i),rrrw(i),ucl(i),urw(i),zcl(i),zrw(i)
 

c	  write(8,*) row(i),col(i)	
c		ra(l)=(atan(nn))*180.0/pi+node(i)
c		if (mm.lt.0.and.ll.gt.0) then
c		   ra(l)=ra(l)+180.000000
c		else 
c		   ra(l)=ra(l)	
c		endif

  10		continue	
      END
