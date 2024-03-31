      subroutine readf(note,j)

c..   READS MODEL USING NEW FORMAT..Thu Jan 19 13:27:39 PST 2006
c..   complement to ritef

      implicit none

      include 'dimenfile'
      include 'comod'

      include 'czone'
      include 'cconst'
      include 'cburn'
      include 'crate'

      integer*4  j, k, n, kkp2, netsize0,lonuc0
      character*72 note
      character*30 dumb

      integer*4 io
      parameter( io = 8 )
c---------------------------------------------------------------------
c..   reads all updated values for j=2,
c..   reads initial values of r,u,t,p,tl,v (updated x,h,dmh) for j=1

      read(io,'(a72)')note
c      write(*,*)note

      read(io,10)dumb,model
c      write(*,10)dumb,model

      read(io,10)dumb,kk
c      write(*,10)dumb,kk

      kkp2 = kk+2

       read(io,10)dumb,mode
c       write(*,10)dumb,mode

       read(io,10)dumb,modec
c       write(*,10)dumb,modec

       read(io,10)dumb,modes
c       write(*,10)dumb,modes

       read(io,10)dumb,modex
c       write(*,10)dumb,modex

       read(io,10)dumb,modez
c       write(*,10)dumb, modez

       read(io,10)dumb,mixmode
c       write(*,10)dumb,mixmode

       read(io,10)dumb,netsize0
c       write(6,10)dumb,netsize0

       read(io,10)dumb,lonuc0
c       write(6,10)dumb,lonuc0

       read(io,10)dumb,nopac
c       write(*,*)dumb,nopac

       read(io,10)dumb,nkscale
       write(*,*)dumb,nkscale

       read(io,10)dumb,noion
       write(*,*)dumb,noion

       read(io,10)dumb,norad
c       write(*,*)dumb,norad

       read(io,10)dumb,nomole
c       write(*,*)dumb,nomole

       read(io,10)dumb,nocoul
c       write(*,*)dumb,nocoul

       read(io,10)dumb,nburn
c       write(*,*)dumb,nburn

       read(io,10)dumb,neutro
c       write(*,*)dumb,neutro

       read(io,10)dumb,jnb
c       write(*,*)dumb,jnb

       read(io,10)dumb,nopaleos
c       write(*,*)dumb,nopaleos

       read(io,10)dumb,mrot
c       write(*,*)dumb,mrot

       read(io,10)dumb,nsweep
c       write(*,*)dumb,nsweep

       read(io,10)dumb,newt
c       write(6,10)dumb,newt

       read(io,10)dumb,mloss
c       write(*,*)dumb,mloss

       read(io,10)dumb,altloss
c       write(*,*)dumb,altloss

c      read(io,10)dumb
c..line having l1 entry, now dummy
c      read(io,10)dumb,l2
c      read(io,10)dumb,l3
c      read(io,10)dumb,ncond
c      read(io,10)dumb,iter
c      read(io,10)dumb,modec
c      read(io,10)dumb,modes
c      read(io,10)dumb,modex
c      read(io,10)dumb,modez

 10   format(a30,i10)

      read(io,11)dumb,time
c      write(*,*)dumb,time

      read(io,11)dumb,dth(2)
c      write(*,*)dumb,dth(2)
c..assume last time step was the same size (to avoid 0)
      dth(1) = dth(2)

      read(io,11)dumb,qcons
c      write(*,*)dumb,qcons

      read(io,11)dumb,xm(1)
c      write(*,*)dumb,xm(1)

      read(io,11)dumb,xm(kk)
c      write(*,*)dumb,xm(kk)

      read(io,11)dumb,xm(kk+1)
c      write(*,*)dumb,xm(kk+1)

      read(io,11)dumb,alphaml
c      write(*,*)dumb,alphaml

      read(io,11)dumb,uuml
c      write(*,*)dumb,uuml

      
      read(io,11)dumb,fthoul
c      write(*,*)dumb,fthoul


       read(io,11)dumb,zhyd0
c       write(*,*)dumb,zhyd0

       read(io,11)dumb,zpop0
c       write(*,*)dumb,zpop0

      read(io,11)dumb,peryear
c      write(*,*)dumb,peryear



c..   space for new entries: 5 lines

       read(io,'(a30,a15)')dumb,copal
c       write(*,*)dumb,copal

       do n = 1, 4
          read(io,'(a30)')dumb
c          write(*,*)dumb

       enddo

c..done in gen.f via read from params.d
c      read(io,11)dumb,resid
c      read(io,11)dumb,temin
c      read(io,11)dumb,cdelt
c      read(io,11)dumb,cdelv
c      read(io,11)dumb,cdeln

 11   format(a30,1pe25.16)

      read(io,12)dumb
c      write(*,12)dumb

c      write(*,*)io,kkp2,j

      read(io,13)(r(j,k),k=1,kkp2)
c      write(*,13)(r(j,k),k=1,kkp2)

      read(io,12)dumb
c      write(*,*)dumb
      read(io,13)(u(j,k),k=1,kkp2)

      read(io,12)dumb
c      write(*,*)dumb
      read(io,13)(omeg(k),k=1,kkp2)

      read(io,12)dumb
c      write(*,*)dumb
      read(io,13)(ajay(k),k=1,kkp2)

      read(io,12)dumb
c      write(*,*)dumb
      read(io,13)(t(j,k),k=1,kkp2)

      read(io,12)dumb
c      write(*,*)dumb
      read(io,13)(v(j,k),k=1,kkp2)

      read(io,12)dumb
c      write(*,*)dumb
      read(io,13)(p(j,k),k=1,kkp2)

      read(io,12)dumb
c      write(*,*)dumb
      read(io,13)(dmh(k),k=1,kkp2)

      read(io,12)dumb
c      write(*,*)dumb
      read(io,13)(h(k),k=1,kkp2)

      read(io,12)dumb
c      write(*,*)dumb
      read(io,13)(tl(j,k),k=1,kkp2)



 12   format(a17)
 13   format(1p5e17.8)

c      cnuc(ndim) = '   Ye'

      if( newnet .eq. 0 )then
c..   new definitions will replace these................
         read(io,'(17i5)')(  lz(n),n=1,ndim)
c         write(*,'(17i5)')(  lz(n),n=1,ndim)


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
            write(*,'(/a40)')'READF: network sizes are inconsistent'
            write(*,*)' make new net.rc using newnet .ne. 0'
            write(*,'(2(a10,i5))')'code',netsize,'imodel',netsize0
            stop'readf'
         endif

         read(io,'(17i5)',err=1000)(  ln(n),n=1,ndim)
c         write(*,'(17i5)')(  ln(n),n=1,ndim)

         read(io,'(17a5)',err=1000)(cnuc(n),n=1,ndim)
c         write(*,'(17a5)')(cnuc(n),n=1,ndim)

         read(io,'(1p6e14.5)')(xxin(n),n=1,ndim)
c         write(*,'(1p6e14.5)')(xxin(n),n=1,ndim)

         do n = 1, netsize+1
            read(io,12)dumb
c            write(*,*)dumb
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
         write(*,*)'readf: newnet ',newnet
         write(*,*)
     1'*** changed network or reset abundances, will try to fix'
      endif

      return

 1000 continue
c..   error in abundance read versus network
      write(*,*)'readf: read error, inconsistent abundances'
      stop'readf'

      end








