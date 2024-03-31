      subroutine nritef(note,j)

      implicit none

      include 'dimenfile'
      include 'comod'
      include 'czone'
      include 'cconst'
      include 'cburn'

c..uses unix intrinsic function fdate

      integer*4  j, k, n
      character*72 note
      character*25 fdate,cdate

      integer*4 io
      parameter( io = 8 )
c------------------------------------------------------------------
c..writes all updated values for j=2,
c..writes initial values of r,u,t,p,tl,v (updated x,h,dmh) for j=1

c..add present date/time at end of header
      cdate = fdate()
      note(47:72)=cdate

       write(io,'(a72)')note
       write(io,10)'model',model
       write(io,10)'kk',kk
       write(io,10)'mode',mode
       write(io,10)'l1'
       write(io,10)'l2',l2
       write(io,10)'l3',l3
       write(io,10)'ncond',ncond
       write(io,10)'iter',iter
       write(io,10)'modec',modec
       write(io,10)'modes',modes
       write(io,10)'modex',modex
       write(io,10)'modez',modez

 10    format(a30,i10)

       write(io,11)'time',time
       write(io,11)'time step',dth(2)

       write(io,11)'qcons',qcons
       write(io,11)'xm(1)',xm(1)
       write(io,11)'resid',resid

       write(io,11)'temin',temin
       write(io,11)'cdelt',cdelt
       write(io,11)'cdelv',cdelv
       write(io,11)'cdeln',cdeln

 11    format(a30,1pe25.16)

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

 12    format(a17)
 13    format(1p5e17.8)

       cnuc(netsize+1) = '   Ye'

       write(io,'(17i5)')(  lz(n),n=1,ndim)
       write(io,'(17i5)')(  ln(n),n=1,ndim)
       write(io,'(17a5)')(cnuc(n),n=1,ndim)

       do n = 1, netsize+1
         write(io,12)cnuc(n)
         write(io,13)(x(n,k),k=1,kk+2)
       enddo

      return
      end

