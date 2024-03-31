c
c
      subroutine ritef(note,j)
c     
c     
      
      implicit none

      include 'dimenfile'
      include 'comod'
      include 'czone'
      include 'cconst'
      include 'cburn'

c     uses unix intrinsic function fdate
      
      integer*4  j, k, n
      character*72 note
      character*25 fdate,cdate
      
      integer*4 io
      parameter( io = 8 )


c     ---------------------------------------
c..   writes all updated values for j=2,
c..   writes initial values of r,u,t,p,tl,v 
c     (updated x,h,dmh) for j=1
c     
c..   add present date/time at end of header
c     ---------------------------------------
      
      
      cdate = fdate()
      note(47:72)=cdate
      
      write(io,'(a72)')note
      write(io,10)'model',model
      write(io,10)'kk',kk
      write(io,10)'mode',mode
      
c     write(io,10)'l1'
c     write(io,10)'l2',l2
c     write(io,10)'l3',l3
c     write(io,10)'ncond',ncond
c     write(io,10)'iter',iter
      write(io,10)'modec',modec
      write(io,10)'modes',modes
      write(io,10)'modex',modex
      write(io,10)'modez',modez
      write(io,10)'mixmode',mixmode
      write(io,10)'netsize',netsize
      write(io,10)'lonuc',lonuc
      write(io,10)'nopac',nopac
      write(io,10)'nkscale',nkscale
      write(io,10)'noion',noion
      write(io,10)'norad',norad
      write(io,10)'nomole',nomole
      write(io,10)'nocoul',nocoul
      write(io,10)'nburn',nburn
      write(io,10)'neutro',neutro
      write(io,10)'jnb',jnb
      write(io,10)'nopaleos',nopaleos
      write(io,10)'mrot',mrot
      write(io,10)'nsweep',nsweep
      write(io,10)'newt',newt
      write(io,10)'mloss',mloss
      write(io,10)'altloss',altloss
      
      
 10   format(a30,i10)
      
      write(io,11)'time',time
      write(io,11)'time step',dth(2)
      
      write(io,11)'qcons',qcons
      write(io,11)'xm(1)',xm(1)
      
      write(io,11)'xm(kk)',xm(kk)
      write(io,11)'xm(kk+1)',xm(kk+1)
      write(io,11)'alphaml',alphaml
      write(io,11)'uuml',uuml
      write(io,11)'fthoul',fthoul
      write(io,11)'zhyd0',zhyd0
      write(io,11)'zpop0',zpop0
      write(io,11)'peryear',peryear
      
      write(io,'(a30,a15)')'opacity tables used',copal
      write(*,*)copal
      
      do n = 1, 4
         write(io,'(a72)')' '
      enddo
      
c     write(io,11)'resid',resid
c     write(io,11)'temin',temin
c     write(io,11)'cdelt',cdelt
c     write(io,11)'cdelv',cdelv
c     write(io,11)'cdeln',cdeln
      
 11   format(a30,1pe25.16)
      
      write(io,12)'radius'
      write(io,13)(r(j,k),k=1,kk+2)
      
      write(io,12)'velocity'
      write(io,13)(u(j,k),k=1,kk+2)
      
      write(io,12)'rotation'
      write(io,13)(omeg(k),k=1,kk+2)
      
      write(io,12)'angular momentum'
      write(io,13)(ajay(k),k=1,kk+2)
      
      write(io,12)'temperature'
      write(io,13)(t(j,k),k=1,kk+2)
      
      write(io,12)'specific volume'
      write(io,13)(v(j,k),k=1,kk+2)
      
      write(io,12)'pressure'
      write(io,13)(p(j,k),k=1,kk+2)
      
      write(io,12)'zone mass'
      write(io,13)(dmh(k),k=1,kk+2)
      
      write(io,12)'convection speed'
      write(io,13)(h(k),k=1,kk+2)
      
      write(io,12)'luminosity'
      write(io,13)(tl(j,k),k=1,kk+2)
      
 12   format(a17)
 13   format(1p5e17.8)
      
      cnuc(netsize+1) = '   Ye'
      
      write(io,'(17i5)')(  lz(n),n=1,ndim)
      write(io,'(17i5)')(  ln(n),n=1,ndim)
      write(io,'(17a5)')(cnuc(n),n=1,ndim)
cccccccccccccccccccccc
      write(io,'(1p6e14.5)')(xxin(n),n=1,ndim)
      
      do n = 1, netsize+1
         write(io,12)cnuc(n)
         write(io,13)(x(n,k),k=1,kk+2)
      enddo
      

      
      
c     SUCCESS
      return
      
      end
