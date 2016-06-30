        PROGRAM RADEC2PIX

c	Converts ra dec to rowprime and colprime.

c	remember to specify with -o name when compiling.

        real*8 a(157),b(157),c(157),d(157),e(157),f(157)
	real*8 alpha,iang(157),mu,nu,row(157),col(157)
        real*8 x,yy,cra,sra,sd,cd,ci,si,ra(157),dec(157)
	
	alpha=95.0	
	pi=3.14159265358
        open(unit=8,file='rdpix.txt', status='unknown')
 3      format(A20)

        open(unit=11,file='radec.txt',form='formatted',status='old')
	do 10 i=1,157
          read(11,*) ra(i),dec(i),iang(i),a(i),b(i),c(i),d(i),e(i),f(i)

	  cra=cos((ra(i)-alpha)*pi/180.0)
	  sra=sin((ra(i)-alpha)*pi/180.0)
c......all of the following are positive.
	  sd=sin(dec(i)*pi/180.0)
	  cd=cos(dec(i)*pi/180.0)
          ci=cos(iang(i)*pi/180.0)
	  si=sin(iang(i)*pi/180.0)
c	  write(*,*) sra, sd,cd,ci,si,' sra sd cd ci si'
	  yy=(ci*sd-si*sra*cd)
	  nu=(asin(yy))*180.0/pi
c 	 write(*,*) nu, 'nu'        
	  x=cra*cd/cos(nu*pi/180.0)
        if (ra(i).lt.90.0) then
		mu=-acos(x)*180.0/pi+alpha + 360.0
        elseif (ra(i).gt.275.0) then
                mu=acos(-x)*180.0/pi+alpha+180.0
        else
                mu=acos(x)*180.0/pi+alpha
        endif

	 write(*,*) ra(i),mu,nu, ' ra,mu nu'
c	 write(*,*) dec(i),nu,d(i),' dec nu d'
	  row(i)=(f(i)*(mu-a(i))-c(i)*(nu-d(i)))/(-c(i)*e(i)+b(i)*f(i))
	  col(i)=(e(i)*(mu-a(i))-b(i)*(nu-d(i)))/(c(i)*e(i)-b(i)*f(i))

c	  write(*,*) row(i)
	  write(8,*) row(i),col(i)	
c		ra(i)=(atan(nn))*180.0/pi+alpha
c		if (mm.lt.0.and.ll.gt.0) then
c		   ra(i)=ra(i)+180.000000
c		else 
c		   ra(i)=ra(i)	
c		endif
  10		continue	
      END
