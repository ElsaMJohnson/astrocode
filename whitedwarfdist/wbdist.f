        PROGRAM WBDIST

        REAL*8 psmal,bb,nn,tt1a,tt2a,tt3a
        REAL*8 tt1b,tt2b,tt3b,xx,yy
        REAL*8 g,nrm,pbig,ttta,tttb
        REAL*8 gmma,eeta,bta,norma,hist,num
        REAL*8 bin(10000000),c2s
        REAL*8 tchi,bbest,nbest,gambest,normabest
        REAL*8 tt1,tt2,tt3,ttt,frac,cnt,cnt2
        CHARACTER*11 winput,cp(100000)
        integer seed
C       this program takes WD mass histogram data 
c       and tries to fit it to a weilbull dist.
c       chi squares it until it works. Needs
c       to have a start in order to guess a better fit.

       write(*,*) 'enter in numbr of mass bin fraction (e.g. .02)'
        read(*,*) frac
 150    format(A11)
        write(*,*) 'enter number of samples'
        read(*,*) wdnum
        write(*,*) 'enter seed for rand function'
        read(*,*) seed

c        seed=398
        
        cnt=0
        cnt2=0

        bb=1.71
        nn=0.21
        nrm=frac
        g=.46
        pbig=0
        psmal=0
c        num=3418
        
        
 100    format(f5.3,1x,f3.0)
c
         c2s=0
         tt1a=0
         tt2a=0
         tt3a=0
         tt1b=0
         tt2b=0
         tt3b=0

        nbins=(1.4-g)/frac
        write(*,*) nbins
        do i=1,nbins
          
          bin(i)=frac*(i)+g
c          write(*,*) bin(i)
        end do
        do i=1,wdnum
        xx=int(nbins*(ran3(seed)))
        yy=int(nbins*(ran3(seed)))
c       write(*,*) xx,yy

c       Weilbull distribution         
        tt1a=nrm*(bb/nn)
        tt2a=(((bin(xx)-g)/nn)**(bb-1))
        tt3a=exp(-((bin(xx)-g)/nn)**bb)
c        write(*,*) nrm,tt1,tt2,tt3
        ttta=tt1a*tt2a*tt3a
        tt1b=nrm*(bb/nn)
        tt2b=(((bin(yy)-g)/nn)**(bb-1))
        tt3b=exp(-((bin(yy)-g)/nn)**bb)

        tttb=tt1b*tt2b*tt3b
c        if (g.le.bin(xx)) then
c        write(*,*) bin(xx)
          sum=bin(xx)+bin(yy)
          pt=ttta*tttb
        if (sum.lt.1.4) then
        write(*,*) sum,pt
        end if
c        end if


        if (sum.ge.1.4) then
        pbig=pbig+pt
        cnt=cnt+1
        else  if (sum.lt.1.4) then
         psmal=psmal+pt
        cnt2=cnt2+1
        end if
        end do
        psmal=psmal/(wdnum)
        pbig=pbig/(wdnum)
        write(*,*) "prob less than 1.4, prob greater than 1.4"
        write(*,*) psmal,pbig,cnt,cnt2     
        
         END

      Real*8 FUNCTION RAN3(IDUM)
      Implicit None
      Integer iff,idum,i,inext,inextp,k,ii
      Real*8 mbig, mseed, mz, fac, ma(55), mj, mk
      save inext,INEXTP, ma

      PARAMETER (MBIG=4000000.D0,MSEED=1618033.D0,MZ=0.D0,FAC=2.5D-7)
CC      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DATA IFF /0/
        IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
         IFF=1
         MJ=MSEED-dble(IABS(IDUM))
         MJ=MOD(MJ,MBIG)
         MA(55)=MJ
         MK=1
         DO I=1,54,1
            II=MOD(21*I,55)
            MA(II)=MK
            MK=MJ-MK
            IF(MK.LT.MZ)MK=MK+MBIG
            MJ=MA(II)
         END DO
         DO K=1,4,1
            DO I=1,55,1
               MA(I)=MA(I)-MA(1+MOD(I+30,55))
               IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
            END DO
         END DO
         INEXT=0
         INEXTP=31
         IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.ge.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.ge.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END
 
