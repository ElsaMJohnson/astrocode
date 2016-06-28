        PROGRAM WEILBULL

        REAL*8 chi(100),bb(10),nn(10)
        REAL*8 g(10),nrm(10),wbb(100)
        REAL*8 gmma,eeta,bta,norma,hist,num
        REAL*8 bin(10000),amp(10000),c2s
        REAL*8 tchi,bbest,nbest,gambest,normabest
        REAL*8 tt1,tt2,tt3,ttt,namp(10000)
        CHARACTER*11 winput,cp(100000)

C       this program takes WD mass histogram data 
c       and tries to fit it to a weilbull dist.
c       chi squares it until it works. Needs
c       to have a start in order to guess a better fit.

c        write(*,*) 'enter in wd info file'
c        read(*,150) winput
 150    format(A11)
        winput='wdinput.txt'

        bta=1.76
        eeta=0.22
        norma=.05
        gmma=.42
        num=3418
        
        normabest=norma
        gambest=gmma
        nbest=eeta
        bbest=bta
        tchi=1000.
        
        open(unit=11,file=winput,form='formatted',status='old')
        hist=0
 4      continue
              read(11,150,end=5000) cp(hist)
              hist=hist+1
              goto 4
 5000         continue
              rewind(11)       

        c2s=0
 100    format(f5.3,1x,f3.0)
c
        do t=1,1
        do m=1,10
        do j=1,10
        do x=1,10
         c2s=0
         tt1=0
         tt2=0
         tt3=0

        do i=1,hist
        read(11,100) bin(i),amp(i) 
        namp(i)=amp(i)/num
          nrm(t)= norma 
          bb(m)=bta -0.05 + 0.01*(m-1)
          nn(j)=eeta -.1 + (j-1)*.01
          g(x) = gmma+( -5.+(x-1))*.01
c           nrm(t)= norma
c           bb(m)=bta
c           g(x) = gmma
c           nn(j)=eeta

c       write(*,*) nrm(t),bb(m),g(x),nn(j)
   
c        write(*,*) nrm(t),bb(m),g(x),nn(j)
        tt1=nrm(t)*(bb(m)/nn(j))
        tt2=(((bin(i)-g(x))/nn(j))**(bb(m)-1))
        tt3=exp(-((bin(i)-g(x))/nn(j))**bb(m))
c        write(*,*) nrm(t),tt1,tt2,tt3
        ttt=tt1*tt2*tt3
        write(*,*) ttt,namp(i)
        
c        write(*,*) 
            chi(i)=((ttt-namp(i))**2)/ttt
        
        if (g(x).lt.bin(i)) then
c        write(*,*) chi(i) 
        c2s=c2s+chi(i)
        end if
        end do
        write(*,*) c2s
        rewind(11)
       if (c2s.lt.tchi) then
       write(*,*) c2s   
        tchi=c2s     
        nbest=nn(j)
        bbest=bb(m)
        normabest= nrm(t)
        gambest=g(x)
        end if
c        write(*,*) c2s
        end do
        end do
        end do
        end do
        write(*,*) 'best parameters: c2s,bta,eeta,norma,gmma'
        write(*,*) tchi,bbest,nbest,normabest,gambest
        
         END 
