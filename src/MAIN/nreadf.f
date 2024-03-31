      subroutine nreadf(note,j)

c..   READS MODEL USING NEW FORMAT...Thu Apr 15 11:31:09 MST 1999
c..   complement to nritef

      implicit none

      include 'dimenfile'
      include 'comod'

      include 'czone'
      include 'cconst'
      include 'cburn'

      integer*4  j, k, n, kkp2, netsize0
      character*72 note
      character*30 dumb

      integer*4 io
      parameter( io = 8 )
c---------------------------------------------------------------------
c..   reads all updated values for j=2,
c..   reads initial values of r,u,t,p,tl,v (updated x,h,dmh) for j=1

      read(io,'(a72)')note
      write(*,'(a72)')note

      read(io,10)dumb,model
      write(*,10)dumb,model

      read(io,10)dumb,kk
      write(*,10)dumb,kk

      kkp2 = kk+2

      read(io,10)dumb,mode

      read(io,10)dumb
c..line having l1 entry, now dummy
      read(io,10)dumb,l2
      read(io,10)dumb,l3
      read(io,10)dumb,ncond
      read(io,10)dumb,iter
      read(io,10)dumb,modec
      read(io,10)dumb,modes
      read(io,10)dumb,modex
      read(io,10)dumb,modez

 10   format(a30,i10)

      read(io,11)dumb,time
      read(io,11)dumb,dth(2)
c..assume last time step was the same size (to avoid 0)
      dth(1) = dth(2)

      read(io,11)dumb,qcons
      read(io,11)dumb,xm(1)
      read(io,11)dumb,resid

      read(io,11)dumb,temin
      read(io,11)dumb,cdelt
      read(io,11)dumb,cdelv
      read(io,11)dumb,cdeln

 11   format(a30,1pe25.16)

      read(io,12)dumb
      read(io,13)(r(j,k),k=1,kkp2)

      read(io,12)dumb
      read(io,13)(u(j,k),k=1,kkp2)

      read(io,12)dumb
      read(io,13)(omeg(k),k=1,kkp2)

      read(io,12)dumb
      read(io,13)(ajay(k),k=1,kkp2)

      read(io,12)dumb
      read(io,13)(t(j,k),k=1,kkp2)

      read(io,12)dumb
      read(io,13)(v(j,k),k=1,kkp2)

      read(io,12)dumb
      read(io,13)(p(j,k),k=1,kkp2)

      read(io,12)dumb
      read(io,13)(dmh(k),k=1,kkp2)

      read(io,12)dumb
      read(io,13)(h(k),k=1,kkp2)

      read(io,12)dumb
      read(io,13)(tl(j,k),k=1,kkp2)

 12   format(a17)
 13   format(1p5e17.8)

c      cnuc(ndim) = '   Ye'

      if( newnet .eq. 0 )then
c..   new definitions will replace these................
         read(io,'(17i5)',err=1000)(  lz(n),n=1,ndim)

c..   check network consistency
c..   only neutrons have Z=0
         netsize0 = 1
         do n = 1, ndim
            if( lz(n) .gt. 0 )then
               netsize0 = netsize0 + 1
            endif
         enddo
         if( netsize .ne. netsize0 )then
c..inconsistent networks
            write(*,'(/a40)')'NREADF: network sizes are inconsistent'
            write(*,*)' make new net.rc using newnet .ne. 0'
            write(*,'(2(a10,i5))')'code',netsize,'imodel',netsize0
            stop'nreadf'
         endif

         read(io,'(17i5)',err=1000)(  ln(n),n=1,ndim)
         read(io,'(17a5)',err=1000)(cnuc(n),n=1,ndim)

         do n = 1, netsize+1
            read(io,12)dumb
            read(io,13)(x(n,k),k=1,kkp2)
         enddo

c-----------------------------------------------------
c     read(8,10)model,iflag
c     read(8,13)ll,runtmx
c     read(8,10)kmax,kk,ifm,ieos,nbug
c     read(8,10)mode,l1,l2,l3,ncond,iter,nfall,nshel,
c     1           modec, modes, modex, modez
c     read(8,11)time , tenvelop
c     read(8,23)dth(1),dth(2),dti(1),idt
c
c     read(8,11)qcons,xm(1),etar,resid
c     read(8,11)temin,cdelt,cdelv,cdeln
c
c     read(8,12)(k,r(j,k),u(j,k),t(j,k),p(j,k),k=1,kmax)
c     read(8,12)(k,dmh(k),h(k),tl(j,k),v(j,k),k=1,kmax)
c
c     read(8,12)(k,(x(n,k),n=1,4),k=1,kmax)
c     read(8,12)(k,(x(n,k),n=5,8),k=1,kmax)
c     read(8,12)(k,(x(n,k),n=9,12),k=1,kmax)
c
c     read(8,12)(k,(x(n,k),n=13,16),k=1,kmax)
c     read(8,12)(k,(x(n,k),n=17,20),k=1,kmax)
c     read(8,12)(k,(x(n,k),n=21,24),k=1,kmax)
c     read(8,12)(k,(x(n,k),n=25,28),k=1,kmax)
c     read(8,12)(k,(x(n,k),n=29,32),k=1,kmax)
c     read(8,12)(k, x(33,k), x(34,k),
c     1     omeg(k),ajay(k),k=1,kmax)
c--------------------------------------------------------------------
      else
         write(*,*)'nreadf: newnet ',newnet
         write(*,*)
     1'*** changed network or reset abundances, will try to fix'
      endif

      return

 1000 continue
c..   error in abundance read versus network
      write(*,*)'nreadf: read error, inconsistent abundances'
      stop'nreadf'

      end








