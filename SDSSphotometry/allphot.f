        PROGRAM ALLPHOT

c       This is like the gkphot program, but it now calculates
c       for all sdss filters (u to z)
c	Algorithm:
c	read in 2phot file
c	read in sdssphot file
c	read in zp file
	 
c        From 2phot take epadu, gain, background, flux, apertures, 
c 	names, 
c	sigma,number of pix (area), magnitude, merr
c        From sdssphot take  background, flux, apertures, names, 
c	sigma,
c        number of pix (area), magnitude, merr
c	From zpfiles take aa kk airmass

c	calculate error for 2phot
c	calculate zeroed magnitude for sdss
c	output filename, surface brightnesses, g-k color and error 
c       for both apertures

c	May have a problem with reading magnitudes,
c	some are undefined for low fluxes
c	can't perform arithmetic operations on characters

c	Okay, changed this around a bit. Since the error and mag
c	from the output files are so close together (no white spaces)
c	I decided not to read them into program.
c	This means I also got rid of checking for pier because it
c	lies right after it. With these small files it's easy to
c	check for this error via: tail -n 2 sphot* > crap. 
c	Just check for any nonzeroes after the mag and mag error
c	columns.

c	Interestingly enough, SDSS doesn't depend on the sky sigma
c	but 2MASS does. so I haven't been reading in the sky
c	sigma values for sdss data. But it really doesn't matter.


        CHARACTER*28 massfile(1000),sdssfile(1000)
        CHARACTER*28 zpfile(1000)
        CHARACTER*12 nmm,nms
        REAL skym,skys,apm(10),aps(10),arm(10),ars(10)
	REAL fxm(10),fxs(10),gm(10),ger(10),ress,smm(10)
	REAL prm(10),prs(10),ter(10),ker(10),tc(10),sms(10)
        REAL sigmam,sigmas,db1(10),db2(10),db3(10),km(10)
        REAL mgm(10),mgs(10),merm(10),mers(10)
        integer num(20000),amt,apn,filt
        REAL aa(5),kk(5),air(5),gain(5),dark(5),skyer(5)
        REAL cn,epadu,snr(10)
        CHARACTER*28 masslist,outfile,sdsslist
	

        write(*,*) 'how many files?'
        read(*,*) amt
	write(*,*) 'One or two apertures?'
	read(*,*) apn
        write(*,*) 'which sdss filter? g=2,r=3,i=4,z=5'
        read(*,*) filt	
 111    format(A23)
        write(*,*) 'enter name of 2mass list'
        read(*,111) masslist
        write(*,*) 'enter name of sdss list'
        read(*,111) sdsslist

        open(unit=7,file=masslist,status='old')
	open(unit=8,file=sdsslist,status='old')
	open(unit=9,file='zplist',status='old')
        if (filt.eq.1) outfile='ukcolor.txt'
        if (filt.eq.2) outfile='gkcolor.txt'
        if (filt.eq.3) outfile='rkcolor.txt'
        if (filt.eq.4) outfile='ikcolor.txt'
        if (filt.eq.5) outfile='zkcolor.txt'
	open(unit=10,file=outfile,status='unknown')
        do s=1,amt
         read(7,*) massfile(s)
	 read(8,*) sdssfile(s)
         read(9,*) zpfile(s)
         open(unit=11,file=massfile(s),status='old')
	 open(unit=12,file=sdssfile(s),status='old')
         open(unit=13,file=zpfile(s),status='old')
c        Exposure time (not factored into big grids and g-k photometry)

          et=2.5*log10(53.907456)
         read(12,100)
 100    format(//////////////////////////////////)
         read(12,101)
 101    format(/////////////////////////////////////////)
	 read(11,102)
 102    format(////////////////////)
	 read(11,103) epadu
 103    format(16x,f3.1)
	 write(*,*) epadu
         read(11,101)
	 read(11,104)
 104     format(////////////)
c	 read(11,105) nmm
c         read(12,105) nms
c 105    format(A23/)
	 read(11,106) skym,sigmam
         read(12,106) skys,sigmas
	 write(*,*) skym, skys
	 write(*,*) sigmam,sigmas 
 106     format(f10.2,f11.1/)
        f=filt
        read(13,107)
 107    format(///)
        do 30 k=1,5
         read(13,*) aa(k),kk(k),air(k),gain(k),dark(k),skyer(k)
 30      continue
         write(*,*) aa(f),kk(f),air(f),gain(f),dark(f),skyer(f)

c            read(13,*) aa,kk,air,gain,dark,skyer
c            write(*,*) aa,kk,air,gain,dark,skyer

	 do i=1,apn
c	 maybe store as 2 different variables?
	  read(11,*) apm(i),smm(i),arm(i),fxm(i),mgm(i)
          read(12,*) aps(i),sms(i),ars(i),fxs(i),mgs(i)
 	  write(*,*) apm(i),smm(i),arm(i),fxm(i),mgm(i)
	  write(*,*) aps(i),sms(i),ars(i),fxs(i),mgs(i)	
c 105      format(g18.7,g15.7,g15.7,A9,A7,i5)

          if (fxm(i).lt.0.7) then
            mgm(i)=100.00
          end if
          if (fxs(i).lt.0.7) then
            mgs(i)=100.00
          end if

c        if (prm(i).ne.0) write(*,*)'pier=',prm(i),s
c        if (prs(i).ne.0) write(*,*)'pier=',prs(i),s
            rs=0.396*0.396
c	Surface brightnesses:
         gm(i)=mgs(i)-(aa(f)+kk(f)*air(f))+
     * 2.5*log10(rs*ars(i))+et
                if(fxs(i).lt.0) gm(i)=200
	    km(i)=mgm(i)+2.5*log10(arm(i))
c	total color:
	 tc(i)=gm(i)-km(i)
c      g,r,i or z band error (from sdss)	
       ger(i)=1.0857*((gm(i)/gain(f))+ars(i)*(dark(f)+
     * skyer(f)))**.5/fxs(i)
c      K band error (from 2mass)
	 cn=sigmam/(16.656530250)
	 write(*,*) 'cn', cn
         db1(i)=(4*arm(i)*.024*cn)**2
	 db2(i)=4*arm(i)*(3.4*cn)**2
	 db3(i)=fxm(i)/(epadu*6)
	 snr(i)=fxm(i)/(db3(i)+db2(i)+db1(i))**.5	
	 ker(i)=1.0857/snr(i)
	write(*,*) 'ker',ker
c	total error
	ter(i)=(ger(i)**2+ker(i)**2)**0.5
	end do
 109    format(A15,8(1x,f9.3))
 110    format(A15,4(1x,f9.3))
	if (apn.gt.1) then
	write(10,109) massfile(s),tc(1),ter(1),tc(2),ter(2),gm(1),
     *  km(1),gm(2),km(2)
	else
	write(10,110) massfile(s),tc(1),ter(1),gm(1),km(1)
	end if
       end do
       END
