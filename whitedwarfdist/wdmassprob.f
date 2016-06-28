        PROGRAM WDMASSPROB
        REAL*8 l(25),m(25),mass1,mass2,sum
        REAL*8 bb,p1,p2,pt,psmal,pbig,ps,lm1,lm2
        real*8 cc(20),cc2(20),ss(50),jj,nn,wdnum
        DATA l(1),l(2),l(3),l(4),l(5),l(6),
     * l(7),l(8),l(9),l(10),l(11),l(12),
     * l(13),l(14),l(15),l(16),l(17),l(18)
c     * l(19),l(20)
c     *,l(21),l(22),
c     * l(23),l(24),l(25)
c    SDSS+others
c     */0.000292569,0.000292569,0.001462844
c     *,0.001462844,0.018139263,0.043300176,
c     *0.053540082,0.131948508,0.189877121,
c     *0.148332358,0.098888239,0.086307782,
c     *0.071679345,0.050029257,0.038326507,
c     *0.023698069,0.016091281,0.007606788,
c     *0.00702165,0.004388531,0.003510825,
c     *0.001755413,0.000585138,0.001170275,0/
c
c     Koester 2001 values
c     * /0,0.006134969,0.030674847,0.036809816,
c     * 0.09202454,0.26993865,0.503067485,
c     * 0.668711656,0.773006135,0.846625767,
c     * 0.877300613,0.901840491,0.932515337,
c     * 0.944785276,0.957055215,0.963190184,
c     * 0.981595092,0.993865031/
     */0,0.005763689,0.034582133,0.086455331,
     *0.121037464,0.30259366,0.55907781,
     *0.74351585,0.809798271,0.858789625
     *,0.899135447,0.925072046,0.948126801
     *,0.962536023,0.971181556,0.976945245
     *,0.985590778,0.9913544676/
c     Holmberg et al 2008 limits
c     * l(19),l(20)/0,0.008264463,0.016528926,
c     *0.024793388,0.033057851,
c     *0.049586777,0.082644628,0.173553719,
c     *0.396694215,0.561983471,0.702479339,
c     *0.785123967,0.842975207,0.900826446,
c     *0.909090909,0.917355372,0.94214876,
c     *0.966942149,0.983471074,0.991735537/
       DATA m(1),m(2),m(3),m(4),m(5),m(6),
     * m(7),m(8),m(9),m(10),m(11),m(12),
     * m(13),m(14),m(15),m(16),m(17),m(18)
c     *, m(19)/
c     *m(20),m(21),m(22),
c     * m(23),m(24),m(25)/
c     *.175,.225,.275,.325,.375,.425,.475,.525,.575,
     */.325,.375,.425,.475,.525,.575,
     *.625,.675,.725,.775,.825,.875,.925,.975,1.025,
     *1.075,1.125,1.175/
c     * 1.075,1.125,1.175,1.225,1.275,1.325,1.375/

c     * m(20)/
c       Holmberg et al masses 2008
c     * .2,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,
c     * .85,.9,.95,1.0,1.05,1.1,1.2,1.25/
        INTEGER seed
        
c       wdnum=1.
        write(*,*) 'enter seed for rand function'
        read(*,*) seed
        write(*,*) 'enter number of samples'
        read(*,*) wdnum
c       seed=398
c      get rid of first point
       xx=ran3(seed)
        write(*,*) xx
        cnt1=0
        cnt2=0
        do j=1,18
         cc(j)=0
         cc2(j)=0
        end do
        do j=1,40
         ss(j)=0
        end do
        pbig=0
        psmal=0
        sht=0
        lim=18
        do t=1,wdnum
 10          xx=1.2*ran3(seed)
 20         yy=1.2*ran3(seed)
        if (xx.lt.0.25) goto 10
        if (yy.lt.0.25) goto 20
 
        write(*,*) xx,yy
        
        do i=1,lim
           lm1=m(i) -.025
           lm2=m(i) +.025
          
          if (i.lt.lim) then
 
           if (xx.gt.lm1.and.xx.le.lm2) then
            mass1=m(i)
            p1=l(i+1)-l(i)
c            write(*,*) xx,l(i),p1,lm1,lm2

            cc(i)=cc(i)+1
           end if
          elseif (i.eq.lim) then
           if (xx.gt.lm1) then
            mass1 = m(i)
            p1=.99711819-l(i)
            cc(i)=cc(i)+1
          end if
          end if
        end do
        lm1=0
        lm2=0
        do i=1,lim
           lm1=m(i) -.025
           lm2=m(i) +.025

          if (i.lt.lim) then
           if (yy.gt.lm1.and.yy.le.lm2) then
            mass2=m(i)
            p2=l(i+1)-l(i)
            cc2(i)=cc2(i)+1
           end if
          elseif (i.eq.lim) then
           if (yy.gt.lm1) then
            mass2 = m(i)
            p2=.99711819-l(i)
            cc2(i)=cc2(i)+1
c          write(*,*) 'p2',p2
           end if
          end if
        end do
c         write(*,*) 'p2,p1',p2,p1
          sum=mass1+mass2
          pt=p1*p2
          write(*,*) "sum, pt",sum,pt
   
          if(sum.ge.1.4) then
          pbig=pbig+pt
          sht=sht+1
          write(*,*) pt,pbig
          end if
          if(sum.lt.1.4) psmal=psmal+pt
          write(*,*) 'psmal and pbig',psmal,pbig 
        end do
c       This normalizes the distribution otherwise you get the amplitude
c       of the distribution
        psmal=psmal/wdnum
        pbig=pbig/wdnum
        write(*,*) "prob less than 1.4, prob greater than 1.4"
        write(*,*) psmal,pbig     

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

