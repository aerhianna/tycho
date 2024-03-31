c     
c     
c
      subroutine readb(note,j)
      
c     
c..   READS MODEL FILE WITH NEW FORMAT
c     Thu Jan 19 13:27:39 PST 2006
c..   (complement to ritef)
c     
c     
      
      implicit none

      include 'dimenfile'
      include 'comod'

      include 'czone'
      include 'cconst'
      include 'cburn'
      include 'crate'

      integer*4  j, k, n, kkp2, netsize0,lonuc0
      character*72 note,cdumb72
      character*30 dumb

      integer*4 io,ierr
      parameter( io = 8 )



c     -------------------------------------
c     reads all updated values for j=2,
c     reads initial values of r,u,t,p,tl,v 
c     (updated x,h,dmh) for j=1
c     -------------------------------------

      read(io,iostat=ierr)note
      if( ierr .ne. 0 )then
         write(*,*)'readb has io error on first read'
         write(*,*)'file may not be binary'
         stop'readb.f'
      endif
c     write(*,*)note
      
      read(io)dumb,model
c     write(*,10)dumb,model
      
      read(io)dumb,kk
c     write(*,10)dumb,kk
      
      kkp2 = kk+2
      
      read(io)dumb,mode
c     write(*,10)dumb,mode
      
      read(io)dumb,modec
c     write(*,10)dumb,modec
      
      read(io)dumb,modes
c     write(*,10)dumb,modes
      
      read(io)dumb,modex
c     write(*,10)dumb,modex
      
      read(io)dumb,modez
c     write(*,10)dumb, modez
      
      read(io)dumb,mixmode
c     write(*,10)dumb,mixmode
      
      read(io)dumb,netsize0
c     write(6,10)dumb,netsize0
      
      read(io)dumb,lonuc0
c     write(6,10)dumb,lonuc0
      
      read(io)dumb,nopac
c     write(*,*)dumb,nopac
      
      read(io)dumb,nkscale
c     write(*,*)dumb,nkscale
      
      read(io)dumb,noion
c     write(*,*)dumb,noion
      
      read(io)dumb,norad
c     write(*,*)dumb,norad
      
      read(io)dumb,nomole
c     write(*,*)dumb,nomole
      
      read(io)dumb,nocoul
c     write(*,*)dumb,nocoul
      
      read(io)dumb,nburn
c     write(*,*)dumb,nburn
      
      read(io)dumb,neutro
c     write(*,*)dumb,neutro
      
      read(io)dumb,jnb
c     write(*,*)dumb,jnb
      
      read(io)dumb,nopaleos
c     write(*,*)dumb,nopaleos
      
      read(io)dumb,mrot
c     write(*,*)dumb,mrot
      
      read(io)dumb,nsweep
c     write(*,*)dumb,nsweep
      
      read(io)dumb,newt
c     write(6,10)dumb,newt
      
      read(io)dumb,mloss
c     write(*,*)dumb,mloss
      
      read(io)dumb,altloss
c     write(*,*)dumb,altloss
      
c     return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     read(io,10)dumb
c..   line having l1 entry, now dummy
c     read(io,10)dumb,l2
c     read(io,10)dumb,l3
c     read(io,10)dumb,ncond
c     read(io,10)dumb,iter
c     read(io,10)dumb,modec
c     read(io,10)dumb,modes
c     read(io,10)dumb,modex
c     read(io,10)dumb,modez
c     
c     
 10   format(a30,i10)
      
      read(io)dumb,time
c     write(*,*)dumb,time
      
      read(io)dumb,dth(2)
c     write(*,*)dumb,dth(2)
c..   assume last time step was the same size (to avoid 0)
      dth(1) = dth(2)
      
      read(io)dumb,qcons
c     write(*,*)dumb,qcons
      
      read(io)dumb,xm(1)
c     write(*,*)dumb,xm(1)
      
      read(io)dumb,xm(kk)
c     write(*,*)dumb,xm(kk)
      
      read(io)dumb,xm(kk+1)
c     write(*,*)dumb,xm(kk+1)
      
      read(io)dumb,alphaml
c     write(*,*)dumb,alphaml
      
      read(io)dumb,uuml
c     write(*,*)dumb,uuml
      
      read(io)dumb,fthoul
c     write(*,*)dumb,fthoul
      
      read(io)dumb,zhyd0
c     write(*,*)dumb,zhyd0
      
      read(io)dumb,zpop0
c     write(*,*)dumb,zpop0

      read(io)dumb,peryear
c     write(*,*)dumb,peryear

      read(io)dumb,copal
c     write(*,*)dumb,copal
      
c..   space for new entries: 4 lines
      do n = 1, 4
         read(io)dumb
c     write(*,*)dumb
      enddo


      
c..   done in gen.f via read from params.d
c     read(io)dumb,resid
c     read(io)dumb,temin
c     read(io)dumb,cdelt
c     read(io)dumb,cdelv
c     read(io)dumb,cdeln
      
 11   format(a30,1pe25.16)
      
      read(io)dumb
c     write(*,*)dumb
      
      read(io)(r(j,k),k=1,kkp2)
c     write(*,13)(r(j,k),k=1,kkp2)
      
      read(io)dumb
c     write(*,*)dumb
      read(io)(u(j,k),k=1,kkp2)
      
      read(io)dumb
c     write(*,*)dumb
      read(io)(omeg(k),k=1,kkp2)
      
      read(io)dumb
c     write(*,*)dumb
      read(io)(ajay(k),k=1,kkp2)
      
      read(io)dumb
c     write(*,*)dumb
      read(io)(t(j,k),k=1,kkp2)
      
      read(io)dumb
c     write(*,*)dumb
      read(io)(v(j,k),k=1,kkp2)
      
      read(io)dumb
c     write(*,*)dumb
      read(io)(p(j,k),k=1,kkp2)
      
      read(io)dumb
c     write(*,*)dumb
      read(io)(dmh(k),k=1,kkp2)
      
      read(io)dumb
c     write(*,*)dumb
      read(io)(h(k),k=1,kkp2)
      
      read(io)dumb
c     write(*,*)dumb
      read(io)(tl(j,k),k=1,kkp2)




c     ------------------------
c     NUCLEAR REACTION NETWORK
c     ------------------------
      
      if( newnet .eq. 0 )then
c..   new definitions will replace these ...
         read(io)( lz(n),n=1,ndim)
c     write(*,*)( lz(n),n=1,ndim)
         
c..   check network consistency
c..   only neutrons have Z=0
         netsize0 = 1
         do n = 1, ndim
            if( lz(n) .gt. 0 )then
               netsize0 = netsize0 + 1
            endif
         enddo
         if( netsize .ne. netsize0 )then
c..   inconsistent networks
            write(*,'(/a40)')'READF: network sizes are inconsistent'
            write(*,*)' make new net.rc using newnet .ne. 0'
            write(*,'(2(a10,i5))')'code',netsize,'imodel',netsize0
            stop'readf'
         endif
         
         read(io,err=1000)(  ln(n),n=1,ndim)
c     write(*,'(17i5)')(  ln(n),n=1,ndim)
         
         read(io,err=1000)(cnuc(n),n=1,ndim)
c     write(*,'(17a5)')(cnuc(n),n=1,ndim)
         
         read(io)(xxin(n),n=1,ndim)
c     write(*,'(1p6e14.5)')(xxin(n),n=1,ndim)

         do n = 1, netsize+1
            read(io)cnuc(n)
c     write(*,*)cnuc(n)
            read(io)(x(n,k),k=1,kkp2)
         enddo
      else
         write(*,*)'readf: newnet ',newnet
         write(*,*)
     .        '*** changed network or reset abundances, will try to fix'
      endif
      

c     
      return
      
      
 1000 continue
c     error in abundance read versus network
      write(*,*)'ERR(readb): read error, inconsistent abundances'
      stop'ERR(readb)'

c     
      end









c     --------------------------------------------------
c     NOTE: moved the following commented code here for
c     readability.
c     --------------------------------------------------
c     
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





