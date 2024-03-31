      subroutine opacgn93 (z,xh,t6,r,kzone,lstep)

c..   wda 4-1-08
c..   added lstep in call sequence above
c..   added lstep and kzone in call sequence for opac1 below

      implicit none

c.....The purpose of this subroutine is to interpolate the data along Z
      save
      integer*4 mx,mz,nrm,nrb,nre,nr,ntm,ntb,nt,kzone,lstep
      integer*4 imd,is,iw,izi,ilo,m1,izz,iz,itime,m,mfm,n,mzz,m2,m3,m4
      integer*4 mk,ihi,mf,i

      real*8 z,xh,t6,r,za,opact,zval,xzf,xx,xa,opk2,zza,zzl,opk
      real*8 xz,t6list,alr,dopactd,dopact,dopacr,t6listf,opl,theta
      real*8 quad,dkapz1,dkapz2,dkapz3,dkapz4,dix,dfsx,dfsr,dfsz,dfs
      real*8 dkapdtr,b,dkapdrt,alrf,alt,a,ak

      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/a1/ mzz, xz(mx,mz,nt,nr),  
     . t6list(nt),alr(nr),n(mx),alt(nt),opk(nt,nr),opk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),dfsz(mz),a(3,mx),b(3),m,mf,xa(mx)
     . ,alrf(nrm),xzf(nt,nr),t6listf(ntm),za(mz)
c..... OPACT- opacity obtained from a quadraric interpolation at
c      fixed log T6 at three values of log R; followed by quadratic
c      interpolation along log T6. Results smoothed bt mixing
c      overlapping quadratics.
c..... DOPACT- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics
c              at fixed R  
c..... DOPACR- is  Dlog(k)/Dlog(R) smoothed by mixing quadratics.
c..... DOPACTD- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics
c               at fixed rho
      common/e/ opact,dopact,dopacr,dopactd
      common/ee1/ opl(mx,nt,nr),xx(mx),zza(mz)
      common/eee1/m1,zval
      real*8 kapz,kapz1,kapz2
      dimension kapz(mz),dkapdtr(mz),dkapdrt(mz)
c--------------------------------------------------------------------

      zval=z
      zzl=z   ! use zzl=log10(.0001+z) for log interpolation
      if(itime .ne. 12345678) then
         do i=1,mz
            zza(i)=za(i)        
c.. use zza=log10(0.0001+za(i)) for log interpolation
         enddo
      endif
      if(itime .ne. 12345678) then
        itime=12345678
      endif

c..   look for close match in table
      do i=1,mz
        if(abs(z-za(i)) .lt. 1.e-7 ) then 
          izz=i
          call opac1(0,izz,xh,t6,r,lstep,kzone)

c          if (opact .gt. 9.0) write (*,'(" logK > 9.0, X=",f7.5," Z=",
c     x    f7.5," T6=",f10.5," R=",e12.4)') xh,z,t6,r
c..   keep opacity above klein-nishina limit
c          theta = 2.2d-3*t6
c          ak = log10( 0.39d0*(xh + 1.0d0 )*0.5d0
c     1         /(1.0d0 + theta) )
c          opact = dmax1( opact, ak )
c.....................................................wda 4/1/06

          return
        endif
      enddo

      ilo=2
      ihi=mz
    8 if(ihi-ilo .gt. 1) then
      imd=(ihi+ilo)/2
      if(z .le. za(imd)+1.e-7) then
        ihi=imd
      else
        ilo=imd
        endif
        go to 8
      endif
      i=ihi
      m1=i-2
      m2=i-1
      m3=i
      m4=i+1
      mfm=m4
c.....check whether Z is near a table limit 
      if((z .le. za(2)+1.e-7) .or. (z .ge. za(mz-1))) mfm=m3
c.....  Check if Z+X interpolation sums exceed unity at needed indices.
c       If so, backup to lower Z indices to perform interpolation.
c       This should work OK, due to density of Z-grid points and the 
c       slow Z variation(except at very small Z)
      if(xh+za(mfm) .gt. 1.0d0) mfm=m3 
        if(xh+za(mfm) .gt. 1.0d0) then
            if(m1 .le. 1) then
              write(*,'("special case: X,Z location not covered by"
     x                 ," logic")')
              stop
            endif 
          m1=m1-1
          m2=m2-1
          m3=m3-1
        endif
c
      izi=0
      do iz=m1,mfm
        izz=iz
        call opac1(izi,izz,xh,t6,r,lstep,kzone)

c        if (opact .gt. 9.0) write (*,'(" logK > 9.0, X=",f7.5," Z=",
c     x       f7.5," T6=",f10.5," R=",e12.4)') xh,z,t6,r
c..   keep opacity above klein-nishina limit
c        theta = 2.2d-3*t6
c        ak = log10( 0.39d0*(xh + 1.0d0 )*0.5d0
c     1       /(1.0d0 + theta) )
c        opact = dmax1( opact, ak )
c.....................................................wda 4/1/06
        izi=1
        kapz(iz)=10.0d0**opact ! converts logK to K
        dkapdtr(iz)=dopact
        dkapdrt(iz)=dopacr
      enddo
      is=0
      iw=1
      kapz1=quad(is,iw,zzl,kapz(m1),kapz(m2),kapz(m3)
     x ,zza(m1),zza(m2),zza(m3))
c..   keep opacity above klein-nishina limit
      theta = 2.2d-3*t6
      ak = ( 0.39d0*(xh + 1.0d0 )*0.5d0
     1     /(1.0d0 + theta) )
      kapz1 = dmax1( kapz1, ak )
c.....................................................wda 4/1/06
      is=1
      dkapz1=quad(is,iw,zzl,dkapdtr(m1),dkapdtr(m2),dkapdtr(m3)
     x ,zza(m1),zza(m2),zza(m3))
      dkapz3=quad(is,iw,zzl,dkapdrt(m1),dkapdrt(m2),dkapdrt(m3)
     x ,zza(m1),zza(m2),zza(m3))
      if (mfm .eq. m3) then
        opact=log10(kapz1)   ! converts K to logK
        dopact=dkapz1
        dopacr=dkapz3
        dopactd=-3.*dopacr+dopact
        is=0

        return

      endif
      is=0
      iw=2
      kapz2=quad(is,iw,zzl,kapz(m2),kapz(m3),kapz(m4)
     x ,zza(m2),zza(m3),zza(m4))

c..   keep opacity above klein-nishina limit
      theta = 2.2d-3*t6
      ak = ( 0.39d0*(xh + 1.0d0 )*0.5d0
     1     /(1.0d0 + theta) )
      kapz2 = dmax1( kapz2, ak )
c.....................................................wda 4/1/06

      is=1
      dkapz2=quad(is,iw,zzl,dkapdtr(m2),dkapdtr(m3),dkapdtr(m4)
     x ,zza(m2),zza(m3),zza(m4))
      dkapz4=quad(is,iw,zzl,dkapdrt(m2),dkapdrt(m3),dkapdrt(m4)
     x ,zza(m2),zza(m3),zza(m4))
      dix=(zza(m3)-zzl)*dfsz(m3)
      opact=log10(kapz1*dix+kapz2*(1.-dix))   ! converts K to logK
      dopact=dkapz1*dix+dkapz2*(1.-dix)
      dopacr=dkapz3*dix+dkapz4*(1.-dix)
      dopactd=-3.*dopacr+dopact
      is=0

      return 
      end
