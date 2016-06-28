        PROGRAM WDMASSPROB2
c       This pulls from the actual distribution
        REAL*8 l(20),m(10000),mass1,mass2,sum
        REAL*8 ptbig,ptsmal,pt,num,ar2(10000),wdnum
        real*8 cc(20),cc2(20),ss(50),jj,nn,pairs
        REAL*8 ar1(10000)
        CHARACTER*28 massfile
        INTEGER seed


       wdnum=10.**6
        write(*,*) 'enter seed for rand function'
        read(*,*) seed
        write(*,*) 'enter mass file - must be one column list'
        read(*,150) massfile

 150    format(A28)
        open(unit=9,file=massfile,form='formatted',status='old')


 101    format(f4.2)
        num=0
 4      continue              
              read(9,101,end=5000) m(num)
              num=num+1
              goto 4
 5000         continue
              rewind(9)       
            write(*,*) 'total events', num 
        do i=1,num
        ar1(i)=0
        ar2(i)=0
        end do
        
c       seed=398
c      get rid of first point
       xx=ran3(seed)
        pt=0
        ptbig=0
        ptsmal=0
        pairs=num/2
        if (mod(num,2).ne.0) then
         pairs=(num-1)/2
        end if
c        do t=1,wdnum
        do t=1,num
 10       xx=int(num*ran3(seed))
c          write(*,*) xx
          if (ar1(xx).ne.0) goto 10
          mass1=m(xx)
          ar1(xx)=mass1
 20        yy=int(num*ran3(seed))
c           write(*,*) yy
        if (ar2(yy).ne.0) goto 20
          mass2=m(yy)
          ar2(yy)=mass2
          pt=(1/num)**2  
          sum=mass1+mass2
         if (sum.ge.1.4) then   
          ptbig=ptbig+pt
         else if (sum.lt.1.4) then
          ptsmal=ptsmal+pt
         end if
        end do        
        write(*,*) "# less than 1.4, # greater than 1.4"
        write(*,*) ptsmal,ptbig        

        END

c    ******************************************
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

