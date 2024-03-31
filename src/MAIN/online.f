      subroutine online(tim2,gamma,del,no,kkman,ktype,ks)
c      subroutine online(tim2,gamma,del,no,l,kkman,ktype,ks)

c..   analysis of new time step for on-line output

      implicit none

      include 'dimenfile'

      include 'compu'
      include 'comod'
      include 'ctmstp'
      include 'cqdv.h'
      include 'cruntm.h'
      include 'cbug'
      include 'czone'
      include 'cnabla'
      include 'cphot'
      include 'cgtr.h'
      include 'cconst'
      include 'conline'
      include 'cenv'
      include 'cburn'
      include 'cenrchk'
      include 'caeps'

      real*8 gamma, del, xlum, area, te4, te2, te, tim2

      integer*4 k
      integer*4 i

      integer*4 ioflag, n
      integer*4 no,kkman,ktype,ks

      data dtcflag/'  ','cv'/

      data cnit/'R','T','V','L',' '/

      data ioflag/0/
c--------------------------------------------------------------------

c      write(*,*)'ONLINE: ncyc'
c      write(*,'(50i2)')(ncyc(k),k=1,kk)
cccccccccccccccccccccccccc

c..   worst zones (last to converge)
      ditmax = 0.0d0
      kitmax = 0
      nitest = 5
      do n = 1, 4
         if( abs( wtest(n) ) .gt. abs(ditmax) )then
            ditmax = wtest(n)
            kitmax = nwtest(n)
            nitest = n
         endif
      enddo
      if( nitest .lt. 1 .or. nitest .gt. 5 )then
         write(*,*)'outline error: nitest=',nitest
         stop'outline'
      endif

      if( modes .eq. 2 )then
         call xcheck(kk+1,x)
      else
         call xcheck(kk,x)
      endif

c..   enrchk only computes energy as a check value

      call enrchk(no)

c..   header..........
      if( l .eq. 1 )then

         write(*,20)"l","model","kk ","no","it","dth","tim",
     1        "Te","xlol", "Tc","L","Ts","cit","kit","cdt","kmn"
     2        ,"jmaxz","echk","X(he4)","ncytot"
         write(3,20)"l","model","kk ","no","it","dth","tim",
     1        "Te","xlol", "Tc","L","Ts","cit","kit","cdt","kmn"
     2        ,"jmaxz","echk","X(he4)","ncytot"
      endif

 20   format(a5,a6,a5,a2,a3,2a10,2a6,3a10,4a4,a6,a10,a10,a6)

c..   save initial energies and mass
      if( l .eq. 1 )then
         bindz = - eint + pote
         dmgc2z = dmgc2 + esou - elum
         echk1 = 0
         xmintl = xm(kk) + dmh(kk+1)
      endif

      echk1o   = echk1
      echk1    = - dmgc2 - esou + elum - echkz
      xmfinl   = xm(kk) + dmh(kk+1)

      if( ioflag .eq. 0 )then
c..   output to screen (standard io = unit 6)

         call stabil(2,gamma,del)

c..   calculate position in h-r diagram, uses ks , rphot, uphot
c..   from hstat,dynam

         xlum = tl(2,ks)
         if( xlum .ne. 0.0d0 )then
            xlol = dlog10( dabs( xlum/sollum ) )
            area = 4.0d0*pi*rphot**2
            if( area .gt. 0.0d0 )then
               te4   = dabs( xlum )/(sigma*area)
               if( te4 .ne. 0.0d0 )then
                  te2   = dsqrt(dabs(te4))
                  te    = dsqrt(te2)
                  telog = dlog10(dabs(te))
               else
                  telog = -10.0d0
               endif
c               xmbol = -2.5d0*xlol + 4.72d0
            endif
         endif

         write(*,10) l,model,kk,no,it,dth(2),tim2,
     1        telog,xlol,t(2,2),tl(2,kk)
     1        ,t(2,kk),cnit(nitest),kitmax, dthflag(ktype+1),kkman
     1        ,jmaxz,echeck-echkz,x(nnuc,2),ncytot

         write(3,10) l,model,kk,no,it,dth(2),tim2,
     1        telog,xlol,t(2,2),tl(2,kk)
     1        ,t(2,kk),cnit(nitest),kitmax, dthflag(ktype+1),kkman
     1        ,jmaxz,echeck-echkz,x(nnuc,2),ncytot

c..old monitoring of convection (slow)
         if( mode .ne. 0 .and. modec .ne. -1 )then
            write(*,'(10(10i1,1x))')(ic(k),k=1,kk)

cc            do i = 1, kk/100
cc               write(*,'(10(10i1,1x),i5)')(ic(k+(i-1)*100),k=1,100),
cc     1              i*100
cc            enddo
cc            write(*,'(10(10i1,1x))')(ic(k),k=(kk/100)*100+1,kk)
cc            write(*,'(10i5)')(kk/100)*100+1,kk

c            write(3,'(10(10i1,1x))')(ic(k),k=1,kk)
c            if( kkman .gt. 0 .and. kkman .le. kk )then
c               write(*,'(/a20,i5,1p8e12.4)')'online',
c     1              kk,t(2,kk),tl(1,kk),tl(2,kk),tl(1,kk+1),tl(2,kk+1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            endif
         endif

      endif

 10   format(i5,i6,i5,i2,i3,1p2e10.2,0pf6.3,0pf7.3,1p3e10.2,
     1     2(a3,i5),i5,1pe10.2,1pe10.3,i6)


      return
      end


