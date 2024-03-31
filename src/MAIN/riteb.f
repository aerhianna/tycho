      subroutine riteb(note,j)

      implicit none

      include 'dimenfile'
      include 'comod'
      include 'czone'
      include 'cconst'
      include 'cburn'

c..   uses unix intrinsic function fdate

      integer*4  j, k, n
      character*72 note,cdumb72
      character*25 fdate,cdate
      character*30 dumb
      character*5  cmodel
c..add a character $ to denote binary
      character*8  filout
      character*1  prefix(2),bfix

      logical tobe

      integer*4 io
      parameter( io = 8 )
      data bfix/'0'/
c------------------------------------------------------------------
c..   writes all updated values for j=2,
c..   writes initial values of r,u,t,p,tl,v (updated x,h,dmh) for j=1
c..   binary format (unformatted)

c..   add present date/time at end of header
      cdate = fdate()
      note(47:72)=cdate

      prefix(1) = note(12:12)
      prefix(2) = note(13:13)

c..   save imodel in labeled file to insure a record of it
      call modflgo(model,cmodel)
c..   add bfix as a prefix to denote binary
      filout = bfix//prefix(1)//prefix(2)//cmodel

      inquire(file=filout,exist=tobe)
      if( .not. tobe )then
         open(8,file=filout,form='unformatted',status='new')
         rewind 8
c         write(*,*)filout,' does not exist, opening it'
      else
         write(*,*)filout,' already exists, skipping overwrite'
         rewind 8
         return
      endif

      write(io)note

      dumb = 'model'
      write(io)dumb,model
      dumb = 'kk'
      write(io)dumb,kk
      dumb = 'mode'
      write(io)dumb,mode
      dumb = 'modec'
      write(io)dumb,modec
      dumb = 'modes'
      write(io)dumb,modes
      dumb = 'modex'
      write(io)dumb,modex

      dumb = 'modez'
      write(io)dumb,modez
      dumb = 'mixmode'
      write(io)dumb,mixmode
      dumb= 'netsize'
      write(io)dumb,netsize
      dumb = 'lonuc'
      write(io)dumb,lonuc
      dumb = 'nopac'
      write(io)dumb,nopac
      dumb = 'nkscale'
      write(io)dumb,nkscale
      dumb = 'noion'
      write(io)dumb,noion
      dumb = 'norad'
      write(io)dumb,norad
      dumb = 'nomole'
      write(io)dumb,nomole
      dumb = 'nocoul'
      write(io)dumb,nocoul
      dumb = 'nburn'
      write(io)dumb,nburn
      dumb = 'neutro'
      write(io)dumb,neutro
      dumb = 'jnb'
      write(io)dumb,jnb
      dumb = 'nopaleos'
      write(io)dumb,nopaleos
      dumb = 'mrot'
      write(io)dumb,mrot
      dumb = 'nsweep'
      write(io)dumb,nsweep
      dumb = 'newt'
      write(io)dumb,newt
      dumb = 'mloss'
      write(io)dumb,mloss
      dumb = 'altloss'
      write(io)dumb,altloss

      dumb = 'time'
      write(io)dumb,time
      dumb = 'time step'
      write(io)dumb,dth(2)
      dumb = 'qcons'
      write(io)dumb,qcons
      dumb = 'xm(1)'
      write(io)dumb,xm(1)

      dumb = 'xm(kk)'
      write(io)dumb,xm(kk)
      dumb = 'xm(kk+1)'
      write(io)dumb,xm(kk+1)
      dumb = 'alphaml'
      write(io)dumb,alphaml
      dumb = 'uuml'
      write(io)dumb,uuml
      dumb = 'fthoul'
      write(io)dumb,fthoul
      dumb = 'zhyd0'
      write(io)dumb,zhyd0
      dumb = 'zpop0'
      write(io)dumb,zpop0
      dumb = 'peryear'
      write(io)dumb,peryear

      dumb = 'opacity tables used'
      write(io)dumb,copal

      cdumb72 = ' '
      do n = 1, 4
         write(io)cdumb72
      enddo

      dumb = 'radius'
      write(io)dumb
      write(io)(r(j,k),k=1,kk+2)

      dumb = 'velocity'
      write(io)dumb
      write(io)(u(j,k),k=1,kk+2)
      dumb = 'rotation'
      write(io)dumb
      write(io)(omeg(k),k=1,kk+2)
      dumb = 'angular momentum'
      write(io)dumb
      write(io)(ajay(k),k=1,kk+2)
      dumb = 'temperature'
      write(io)dumb
      write(io)(t(j,k),k=1,kk+2)
      dumb ='specific volume'
      write(io)dumb
      write(io)(v(j,k),k=1,kk+2)

      dumb = 'pressure'
      write(io)dumb
      write(io)(p(j,k),k=1,kk+2)
      dumb = 'zone mass'
      write(io)dumb
      write(io)(dmh(k),k=1,kk+2)
      dumb = 'convection speed'
      write(io)dumb
      write(io)(h(k),k=1,kk+2)
      dumb = 'luminosity'
      write(io)dumb
      write(io)(tl(j,k),k=1,kk+2)

      cnuc(netsize+1) = '   Ye'

      write(io)(  lz(n),n=1,ndim)
      write(io)(  ln(n),n=1,ndim)
      write(io)(cnuc(n),n=1,ndim)
      write(io)(xxin(n),n=1,ndim)

      do n = 1, netsize+1
         write(io)cnuc(n)
         write(io)(x(n,k),k=1,kk+2)
      enddo


      write(*,*)'leaving riteb, ',filout,' written,',
     1 ' to restart: cp ',filout,' imodel'
cccccccccccccccccccccccccccccccccc

      return
      end

