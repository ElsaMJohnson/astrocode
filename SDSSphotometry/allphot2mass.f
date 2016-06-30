        PROGRAM ALLPHOT2MASS

c	Algorithm:
c	read in k,h,j band files and calculate colors (for 2.5 and 3.0
c       apertures)

	 

        CHARACTER*28 massfile(1000),ufile(1000),gfile(1000)
        CHARACTER*12 nmm,nms
        REAL*8 apmk,armk,fxmk,mgmk
        REAL*8 apmj,armj,fxmj,mgmj
        REAL*8 apmh,armh,fxmh,mgmh

        REAL*8 skyk,skyh,skyj,ebv
        REAL*8 sigmak,sigmaj,sigmah,ker,jer,her
        REAL*8 jh,jk,hk,hker,jher,jker
        
        REAL*8 db1k,db2k,db3k,km
        REAL*8 db1j,db2j,db3j,jm
        REAL*8 db1h,db2h,db3h,hm

        integer num(20000),amt,apn,filt
        REAL*8 cnk,epaduk,snrk,zptk
        REAL*8 cnj,epaduj,snrj,zptj
        REAL*8 cnh,epaduh,snrh,zpth
        real*8 ksc,kskb,jsc,jskb,hsc,hskb

        real*8 kpass,jpass,hpass,crapj,crapk,craph

        CHARACTER*12 snfile(1000),kfile(1000),efile(1000)
        CHARACTER*28 masslist,outfile,sdsslist
        CHARACTER*12 tfile(1000),hfile(1000),jfile(1000)
	

        write(*,*) 'how many files?'
        read(*,*) amt
 111    format(A23)
        apn=1

        open(unit=7,file='masslist',status='old')
	open(unit=10,file='hjkphot.txt',status='unknown')
        

        do s=1,amt
          read(7,*) tfile(s)
          kfile(s)='ksm'//tfile(s)
          jfile(s)='jsm'//tfile(s)
          hfile(s)='hsm'//tfile(s)
          efile(s)='e'//tfile(s)

         write(*,*) efile(s)
         open(unit=11,file=kfile(s),status='old')
         open(unit=12,file=jfile(s),status='old')
         open(unit=13,file=hfile(s),status='old')
         open(unit=18,file=efile(s),status='old')
c        Exposure time (not factored into big grids and g-k photometry)

         read(18,*) ebv

 100    format(//////////////////////////////////)

 101    format(/////////////////////////////////////////)
 120    format(//////////////////////////////)
	 read(11,102)
         read(12,102)
         read(13,102)
 102    format(////////////////////)
	 read(11,103) epaduk
         read(12,103) epaduj
         read(13,103) epaduh

 103    format(16x,f3.1)
 121    format(16x,f7.4) 
 122    format(/)
 	 write(*,*) epaduk, epaduj,epaduh
         read(11,120)
         read(12,120)
         read(13,120)
         read(11,121) zptk
         read(12,121) zptj
         read(13,121) zpth
         write(*,*) zpt
	 read(11,102)
         read(11,122)
         read(12,102)
         read(12,122)
         read(13,102)
         read(13,122)

 104     format(////////////)
	 read(11,108) skyk,sigmak
         read(12,108) skyj,sigmaj
         read(13,108) skyh,sigmah

 108    format(f10.2,f11.1//)
 106    format(f10.2,f11.1////) 

        f=filt
 107    format(///)
      
c	 maybe store as 2 different variables?
	  read(11,*) apmk,smmk,armk,fxmk
          read(12,*) apmj,smmj,armj,fxmj
          read(13,*) apmh,smmh,armh,fxmh
        write(*,*) apmk,apmj,apmh

c 105      format(g18.7,g15.7,g15.7,A9,A7,i5)



c       Using UKIRT value for extinction (based on 2.2um lambda
c       where as 2mass is 2.17um in K band) But dust value is so
c       small it wont matter. (saw in McIntosh et al 2006)

c       1.84 is the conversion from Vega to AB system (SDSS)
c       for K band
c	km=mgm-0.367*ebv+2.5*log10(arm)+1.84
c       km=mgm-0.367*ebv+1.84
c      The conditional statements are for crappy weak fluxes.

        if (fxmk.gt.0) then
        mgmk=-2.5*log10(fxmk)+zptk   
        km=mgmk-0.367*ebv
        else if (fxmk.lt.0) then
        km=99.
        end if

c       for J band
c        jm=mgm-0.902*ebv+.9
        if (fxmj.gt.0) then
        mgmj=-2.5*log10(fxmj)+zptj
        jm=mgmj-0.902*ebv
        else if (fxmj.lt.0) then
        jm=99.
        end if

c        for h band
        if (fxmh.gt.0) then
        mgmh=-2.5*log10(fxmh)+zpth
        hm=mgmh-0.576*ebv
        else if (fxmh.lt.0) then
        hm=99.
        end if


        kskb=zptk-2.5*log10(skyk*armk)+2.5*log10(armk)
        ksc=zptk-2.5*log10(smmk)+2.5*log10(armk)
        jskb=zptj-2.5*log10(skyj*armj)+2.5*log10(armj)
        jsc=zptj-2.5*log10(smmj)+2.5*log10(armj)
        hskb=zpth-2.5*log10(skyh*armh)+2.5*log10(armh)
        hsc=zpth-2.5*log10(smmh)+2.5*log10(armh)

c      K band error (from 2mass)
	 cn=1/(16.656530250)
c        Got rid of the above because it diminishes the errors by too 
c       much and its not used in the example provided on the page
        cnk=sigmak*cn
        cnj=sigmaj*cn
        cnh=sigmah*cn

	 write(*,*) 'cn', cnk
c         db1k=(armk*.024*cnk)**2
        
	 db2k=armk*(3.4*cnk)**2
	 db3k=abs(fxmk)/(epaduk*6)
	 snrk=abs(fxmk)/(db3k+db2k)**.5	
	 ker=1.0857/snrk
	write(*,*) 'ker',ker

         db2j=armj*(3.4*cnj)**2
         db3j=abs(fxmj)/(epaduj*6)
         snrj=abs(fxmj)/(db3j+db2j)**.5  
         jer=1.0857/snrj
        write(*,*) 'ker',jer

         db2h=armh*(3.4*cnh)**2
         db3h=abs(fxmh)/(epaduh*6)
         snrh=abs(fxmh)/(db3h+db2h)**.5  
         her=1.0857/snrh
        write(*,*) 'her',her

        crapj=jskb-3*jer
        craph=hskb-3*her
        crapk=kskb-3*ker
        
        if (snrj.ge.7) then
           jpass=1.
        elseif (snrj.lt.7) then
           jpass=0.
        end if

        if (snrh.gt.7) then
           hpass=1.
        else if (snrh.lt.7) then
           hpass=0.
        end if
        if (snrk.ge.7) then
           kpass=1.
        else if (snrk.lt.7) then
           kpass=0.
        end if

        jh=jm-hm
        jher=((jer)**2.+(her)**2.)**0.5
        hk=hm-km
        hker=((her)**2.+(ker)**2.)**0.5
        jk=jm-km
        jker=((jer)**2.+(ker)**2.)**0.5
        
       
        
 109    format(A15,8(1x,f9.3))
 110    format(A15,12(1x,f9.3))
	write(10,110) tfile(s),jm,jer,hm,her,km,ker,jh,jher,
     *  hk,hker,jk,jker
 112   format(A15,3f5.2)
        write(99,112)tfile(s),jpass,hpass,kpass

       end do
       END
