      subroutine opac (z,xh,xxc,xxo,t6,r)

c.....The purpose of this subroutine is to interpolate log kappa
c     and obtain smooth derivatives.
c     in C/O abundance and in T6, R,i.e. (xc,xo,T6,R)
c     xxc=Xc=carbon mass fraction
c     xxo=Xo=oxygen mass abundance
c     t6=T6=temperature in millions of degrees kelvin
c     r=R=density(g/cm**3)/T6**3
c.....to use opac insert common/e/ in the calling routine.
c     This common contains interpolated values for kappa and its
c     first derivities.
c..   made double precision, passes through if conduction dominant
c..   (wda 4-2-01)
c..   t6list = list of temperatures in K/1e6
c..   alt = log 10 of these (vector)
c..   alr = log 10 of Ropal list (vector)
c..   nta = maximum nonzero element in opacity array in temperature
c..   slr = log10 of r (Ropal) variable
c..   slt = log10 of T6 variable

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

      save
c     integer w
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb
     .     ,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)

c      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),
c     .     opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
cccccccccccccccccccccc
      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),
     .     opl(mx,nt,nr),cxdf(mx,mc),oxdf(mx,mc),itime(mx)

      common/aa/ q(4), h(4), xcd(mc),xod(mc), xc(mc),xo(mo)
     .     ,xcs(mc),xos(mo), cxd(mc),oxd(mo),cx(mc),ox(mo),
     .     zzz,xxh,xx(mx),nc,no
c      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr), index(101),
c     .     t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),
c     .     opk(nt,nr),dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),
c     .     alrf(nrm),cof(ntm,nrm),t6listf(ntm),opk2(nt,nr),dfsx(mx)

      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr),
     . t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr)
     . ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) 
     .,t6listf(ntm),opk2(nt,nr),dfsx(mx), index(101)

c      common/b/ itab(mx,ntabs),nta(nrm),x(mx,ntabs),y(mx,ntabs),
c     .     zz(mx,ntabs),xca(mx,ntabs),xoa(mx,ntabs)
ccccccccccccccccccc
      common/b/ x(mx,ntabs),y(mx,ntabs),zz(mx,ntabs),
     1     xca(mx,ntabs),xoa(mx,ntabs),itab(mx,ntabs),nta(nrm)

      common/d/dkap
      common/bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq,xodp,xcdp,xxco,cxx,oxx
c.....OPACT- opacity obtained from a quadraric interpolation at
c     fixed log T6 at three values of log R; followed by quadratic
c     interpolation along log T6. Results smoothed by mixing
c     overlapping quadratics.
c.....DOPACT- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics.
c.....DOPACR- is  Dlog(k)/Dlog(R) smoothed by mixing quadratics.
      common/e/ opact,dopact,dopacr,dopactd
      common/recoin/ itimeco,mxzero
c-------------------------------------------------------------
c.....INDEX refers to the index,i, of xc(i) or xo(i); xc(i) and xo(i)
c     are the abundance grid points.
      iop=1
c..   provides smoothed interpolations; iop=0 gives no smoothing
      if(nr .lt. 6) go to 65
      if((xh .gt. 1.d-6) .and. (mx .lt.4)) go to 64

c.....set-up C/O axis points
      xxco=xxc+xxo
      if(z+xh+xxco-1.d-6 .gt. 1.0d0 ) go to 61
      zzz=z+0.001d0
      xxh=xh
      xxci=xxc
      xxoi=xxo
      xxi=xh
      t6i=t6
      ri=r

c.....convert xxc and xxo to logarithmic shifted by Z
      cxx=log10(zzz+xxc)
      oxx=log10(zzz+xxo)
      xxx=log10(0.005d0+xh)
      slt=log10(t6)
      slr=log10(r)

c.....set X indices
      ilo=2
      ihi=mx
 8    if(ihi-ilo .gt. 1) then
         imd=(ihi+ilo)/2
         if(xh .le. xa(imd)+1.d-7) then
            ihi=imd
         else
            ilo=imd
         endif
         go to 8
      endif
      i=ihi
      mf=i-2
      mg=i-1
      mh=i
      mi=i+1
      mf2=mi
      istep1=1
      if (mx .gt. 1) then
         istep1=mx-1
         if((xh .le. xa(2)+1.d-7) .or. (xh .ge. xa(istep1)-1.d-7))
     1        mf2=mh
      endif

      if ((mx .eq. 1) .or. (xh .lt. 1.d-6)) then
         mf=1
         mg=1
         mh=1
         mi=2
         mf2=1
      endif

      if (itime(1) .ne. 12345678) then
         alr(1)=-8.d0+(nrb-1)*0.5d0
         do i=2,nr
            alr(i)=alr(i-1)+0.5d0
         enddo
         alt(1)=-2.25d0+(ntb-1)*0.05d0
         do i=ntb+1,46
            alt(i)=alt(i-1)+.05d0
         enddo
         ntd=47
         if ((ntb +1) .gt. 47) ntd=ntb+1
         do i=ntd,68
            alt(i)=alt(i-1)+.1d0
         enddo
         do i=68,70
            alt(i)=alt(i-1)+.2d0
         enddo
         do i=1,nt
            t6list(i)=10.d0**alt(i)
         enddo
      endif
      ilo=2
      ihi=nr
 12   if(ihi-ilo .gt. 1) then
         imd=(ihi+ilo)/2
         if(slr .le. alr(imd)+1.d-7) then
            ihi=imd
         else
            ilo=imd
         endif
         go to 12
      endif
      i=ihi
      l1=i-2
      l2=i-1
      l3=i
      l4=l3+1

      ilo=2
      ihi=nt
 11   if(ihi-ilo .gt. 1) then
         imd=(ihi+ilo)/2
         if(t6 .le. t6list(imd)+1.d-7) then
            ihi=imd
         else
            ilo=imd
         endif
         go to 11
      endif
      i=ihi
      k1=i-2
      k2=i-1
      k3=i
      k4=k3+1
      k3s=k3+ntb-1
      l3s=l3+nrb-1

c-----set-up indices to skip when low T&R data missing for X=0.
      kmin=0
      k1in=k1
      iadvance=0
      mfin=mf
      if ((mfin .eq. 1) .and. (co(1,1,1,k1,l1) .gt. 9.d0)) then
c..   data missing
         do i=1,6
            if (co(1,1,1,i,l1) .gt. 9.d0)  then
               if (xh .lt. .1) then
                  kmin=i+1
               else

                  if (iadvance .eq. 0) then
c..   shift X index to avoid X=0.
                     iadvance=iadvance+1
                     mf=mf+1
                     mg=mg+1
                     mh=mh+1
                     mi=mi+1
                     mf2=mf2+1
                  endif
               endif
            endif
         enddo
         if ((iadvance .eq. 0) .and. (k1 .le. kmin) .and.
     x        (slt .le. alt(kmin))) then
            k1=kmin
            if ((co(1,1,1,kmin,l1+1) .lt. 9.d0) .and.
     x           ((slr+.01d0) .gt. alr(l1+1))) then
               l1=l1+1
               kmin=0
               k1=k1in
               do i=1,6
                  if (co(1,1,1,i,l1) .gt. 9.d0) kmin=i+1
               enddo
               if ((kmin .ne. 0) .and. (k1in .lt. kmin)) k1=kmin
            endif
         endif
         if ((slt+.001d0) .lt. alt(k1)) then
            write (*,'("OPAL data not available for X=", f7.5,
     1           " logT6=", f7.3," logR=",f7.3)') xh,slt,slr
            opact=30.d0
            dopact=99.d0
            dopacr=99.d0
            dopactd=99.d0
            return
         endif

         l2=l1+1
         l3=l2+1
         l4=l3+1
         l3s=l3+nrb-1
         k2=k1+1
         k3=k2+1
         k4=k3+1
         k3s=k3+ntb-1
      endif

c-----end of check for missing data

      do 123 m=mf,mf2
         if(mx .ge. 4) then
c.....C and O  fractions determined by the ray through the origin that
c     also passes through the point (Xc,Xo). Specific interpolation
c     values determined by tabulated X values;i.e. xa(m).  Inter-
c     polation along the ray gives log (kappa(Xc,Xo)).  (Advantage
c     of method: keeps indices within table boundaries)
c     Subtract Z to prevent out-of-range C+O values for small X
            if(1.d0-xh-z.gt.1.d-6)then
               cmod=(1.d0-xa(m)-z)/(1.d0-xh-z)
            else
               cmod=0.d0
            endif
            xxc=cmod*xxci
            xxo=cmod*xxoi
            cxx=log10(zzz+xxc)
            oxx=log10(zzz+xxo)
         endif
c..   ythism=z+xa(m)+xxc+xxo

         do i=1,mc
            xhe=1.d0-xa(m)-z
            nc=i
            no=i
            xc(i)=xcs(i)
            xo(i)=xos(i)
c..   if(xcs(i) .ge. xhe-1.e-6) then
            if(xcs(i) .gt. xhe) then
               xc(i)=xhe
               xo(i)=xhe
               go to 3
            endif
         enddo
    3    continue

         if(itime(m) .ne. 12345678) then
            itime(m)=12345678
            mxzero=0
            do i=1,mx
               xx(i)=log10(0.005d0+xa(i))
               if(xa(i) .eq. 0.0d0) mxzero=i
            enddo

c.....this is the first time through this m. Calculate the decadic
c     log of the perimeter points shifted by Z+0.001(to avoid
c     divergence
c     at origin); m refers to xa(m); the hydrogen table value.

c     note that the nc-th elements are sometimes needed!
            do i=1,nc
               oxf(m,i)=log10(zzz+xo(i))
               cxf(m,i)=log10(zzz+xc(i))
               xcdf(m,i)=-xo(i)+xo(no)
               xodf(m,i)=-xc(i)+xc(nc)
               cxdf(m,i)=log10(zzz+xcdf(m,i))
               oxdf(m,i)=log10(zzz+xodf(m,i))
            enddo

c     note that the nc-th elements are sometimes needed!
            do i=1,nc
               ox(i)=oxf(m,i)
               cx(i)=cxf(m,i)
               xcd(i)=xcdf(m,i)
               xod(i)=xodf(m,i)
               cxd(i)=cxdf(m,i)
               oxd(i)=oxdf(m,i)
            enddo

            if( m .le. 0 .or. m .gt. 5 )then
c..test for sane values
               write(*,*)'readco error in opac.f, m=',m
               stop'opac: readco'
            endif
c.....read the data files
            call readco

         endif

         do i=1,nc
            ox(i)=oxf(m,i)
            cx(i)=cxf(m,i)
            xcd(i)=xcdf(m,i)
            xod(i)=xodf(m,i)
            cxd(i)=cxdf(m,i)
            oxd(i)=oxdf(m,i)
         enddo

c.....Determine log R and log T6 grid points to use in the
c     interpolation.
         if((slt .lt. alt(1)).or.(slt .gt. alt(nt))) go to 62
         if((slr .lt. alr(1)).or.(slr .gt. alr(nr))) go to 62

         if (m .eq. mf) then
c..   calculate table indices
            if((mf2 .ne. mxzero) .and. (k3s .gt. ntm)) go to 62
            do i=14,18
               if((l3s .gt. i) .and. (k3s .gt. nta(i+1))) go to 67
            enddo
            ip=3
            iq=3
            ntlimit=nta(l3s)
            if((k3s .eq. ntlimit) .or. (iop .eq. 0)) then
               ip=2
               iq=2
            endif
            if(t6 .le. t6list(2)+1.d-7 .or. iop .eq. 0) ip=2

            if((l3 .eq. nr) .or. (iop .eq. 0)) then
c..   right edge of full table
               iq=2
               ip=2
            endif
            if(slr .le. alr(2)+1.d-7 .or. iop .eq. 0) iq=2
         endif

         xodp=max(-xxc+xc(nc),0.0d0)
         xcdp=max(-xxo+xo(no),0.0d0)
         is=0

         call cointerp(xxc,xxo)
 123  continue

      if((zz(mg,1) .ne. zz(mf,1)) .or. (zz(mh,1) .ne. zz(mf,1))) then 
         write(*,'("Z does not match Z in codata* files you are"
     x        ," using, line 444, xcotrin21.f")')
         stop
      endif
      if(z .ne. zz(mf,1)) go to 66
      xxc=xxci
c..   restores input value; necessary if stop replaced with return
      xxo=xxoi
c..   restores input value
      is=0
      iw=1
      do 45 ir=l1,l1+iq
         do 46 it=k1,k1+ip
            if((mx .eq. 1) .or. (mf2 .eq. 1)) then
               opk(it,ir)=opl(mf,it,ir)
               go to 46
            endif
            opk(it,ir)=quad(is,iw,xxx,opl(mf,it,ir),opl(mg,it,ir)
     x           ,opl(mh,it,ir),xx(mf),xx(mg),xx(mh))
            is=1
 46      continue
 45   continue

      if (mi .eq. mf2) then     ! interpolate between quadratics
         is=0
         iw=1
         dixr=(xx(mh)-xxx)*dfsx(mh)
         do 47 ir=l1,l1+iq
            do it=k1,k1+ip
               opk2(it,ir)=quad(is,iw,xxx,opl(mg,it,ir),opl(mh,it,ir)
     x              ,opl(mi,it,ir),xx(mg),xx(mh),xx(mi))
               opk(it,ir)=opk(it,ir)*dixr+opk2(it,ir)*(1.d0-dixr)
               is=1
            enddo
 47      continue
c..   interpolate X between two overlapping quadratics
      endif

      is=0

c.....completed H,C,O interpolation. Now interpolate T6 and log R on a 
c     4x4 grid. (log(T6(i)),i=i1,i1+3),log(R(j)),j=j1,j1+3)).Procedure 
c     mixes overlapping quadratics to obtain smoothed derivatives.

      call t6rinterp(slr,slt)

      return

c..error and exception exits....................................
 61   write(*,'("opac: Mass fractions exceed unity")')
      write(*,'(8a12)')'z','xh','xxco','T(K)','Ropal'
      write(*,'(1p8e12.3)')z,xh,xxco,t6*1.0d6,r
      xxc=xxci
c..   restores input value; required if stop replaced
c..   with a return
      xxo=xxoi
c..   restores imput value
      stop'opac mass error'
 62   write(*,'("opac: T6/LogR outside of table range",1p3e11.3)')
     1     t6,r,dlog10(t6)+6.0d0
      write(*,'(5(a10,1pe12.3))')'T6',t6,'ropal',r,'log T',
     1     dlog10(t6)+6.0d0,'log rho',dlog10(r)+3.0d0*dlog10(t6),
     2     'logR',dlog10(r)
      write(*,'(5(a10,1pe12.3))')'z',z,'xh',xh,'xxc',xxc,'xxo',xxo
      xxc=xxci
c..   restores input value; required if stop replaced with a return
      xxo=xxoi
c..   restores input value
      stop'opac range error'

 64   write(*,'("opac: X not equal to zero: To run this case it
     .     is necessary"/ "to recompile with parameter (mx=5)")')
      stop'opac X error'
 65   write(*,'("opac: Too few R values; NRE+1-NRB < 6")')
      stop'opac R values error'
 66   write(*,'("opac: Z does not match Z in codata* files you are",
     .     " using, 66 xcotrin21.f")')
      stop'opac Z match error'

 67   continue
c..   conduction region (wda 4-2-01)
c..   the entries (log T6, log ropal) are in the table range 
c..   but the tabular entries for log opacity are blank, so we
c..   approximate the opacity near thomson (klein-nishina) limit

      call cofill(z,xh,xxc,xxo,t6,r,slt,slr,i,l3s,k3s,mf2)

      xxc=xxci
c..   restores input value; required if stop replaced with a return
      xxo=xxoi
c..   restores input value
      return

      end


      subroutine bile(xx,yy,x,y,z,zz,zx,zy)
      implicit none

c..   does bilinear interpolation in x,y square on variable z
c..   xx and yy are coordinates, and must lie inside x,y square
c..   zz is the bilinearly interpolated function
c..   zx and zy are the d(zz)/dx and d(zz)/dy derivatives

      real*8 x(2),y(2),z(2,2),xx,yy,zz,zx,zy
      real*8 denom,dx1,dx2,dy1,dy2
      real*8 dxmin,dxmax,dymin,dymax

c---------------------------------------------------------------------
c..test ranges
      dxmin =  dmin1(x(1),x(2))
      dxmax =  dmax1(x(1),x(2))
      dymin =  dmin1(y(1),y(2))
      dymax =  dmax1(y(1),y(2))

      if( xx .lt. dxmin -1.0d-6 .or. xx .gt. dxmax+1.0d-6 )then
         write(*,'(a30,1p3e15.7)')'bile.f: x out of range: x1,x,x2 ',
     1        x(1),xx,x(2)
         stop'bile.f: x out of range'
      endif


      if( yy .lt. dymin -1.0d-6 .or. yy .gt. dymax +1.0d-6)then
         write(*,'(a30,1p3e15.7)')'bile.f: y out of range: y1,y,y2 ',
     1        y(1),yy,y(2)
         stop'bile.f: y out of range'
      endif


      denom = ( x(2)-x(1) )*( y(2)-y(1) )
      if( denom .eq. 0.0d0 )then
         write(*,*)' bile.f: demoninator = ',denom
         write(*,*)' dx = ', x(2)-x(1),'    dy = ', y(2)-y(1)
         stop'bile error: singular'
      endif

      dx1 = xx - x(1)
      dx2 = xx - x(2)
      dy1 = yy - y(1)
      dy2 = yy - y(2)
      zz  = ( z(1,1)*dx2*dy2 - z(1,2)*dx2*dy1 - z(2,1)*dx1*dy2
     1     + z(2,2)*dx1*dy1 ) / denom
      zx  = ( z(1,1)*dy2 - z(1,2)*dy1 - z(2,1)*dy2
     1     + z(2,2)*dy1 ) / denom
      zy  = ( z(1,1)*dx2 - z(1,2)*dx2 - z(2,1)*dx1
     1     + z(2,2)*dx1 ) / denom
      return
      end



      subroutine cofill(z,xh,xxc,xxo,t6,r,slt,slr,i,l3s,k3s,mf2)

c..darnett 6-6-02
c..fills in blank elements in co array for extrapolation into
c..conduction-dominated region

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

      save
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb
     .     ,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),
     .     opl(mx,nt,nr),cxdf(mx,mc),oxdf(mx,mc),itime(mx)
      common/aa/ q(4), h(4), xcd(mc),xod(mc), xc(mc),xo(mo)
     .     ,xcs(mc),xos(mo), cxd(mc),oxd(mo),cx(mc),ox(mo),
     .     zzz,xxh,xx(mx),nc,no
      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr),
     . t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr)
     . ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) 
     .,t6listf(ntm),opk2(nt,nr),dfsx(mx), index(101)
      common/b/ x(mx,ntabs),y(mx,ntabs),zz(mx,ntabs),
     1     xca(mx,ntabs),xoa(mx,ntabs),itab(mx,ntabs),nta(nrm)
      common/d/dkap
      common/bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq,xodp,xcdp,xxco,cxx,oxx
c.....OPACT- opacity
c.....DOPACT- is Dlog(k)/Dlog(T6)
c.....DOPACR- is  Dlog(k)/Dlog(R)
      common/e/ opact,dopact,dopacr,dopactd
      common/recoin/ itimeco,mxzero

      real*8 xint(2),yint(2),zint(2,2)
c-------------------------------------------------------------
c      if( xxc+xxo .lt. 1.0d-7 )then
c..non-alpha enhanced case
         do ivar = 1, nr
            if( nta(ivar) .lt. nt )then
               ntlow = ivar-1
               go to 100
            endif
         enddo
 100     continue

         do m = mf,mf2
            if(  abs( co(m,1,1,nt,ntlow+1)) .lt. 1.0d-6 .and.
     1           abs( co(m,1,1,nt,nr)) .lt. 1.0d-6 )then
c..   fill in co array
               do j = ntlow+1,nr
                  co(m,1,1,nt,j) =  co(m,1,1,nt,ntlow)
               enddo
c..   work out to larger log ropal
               do l = ntlow+1,nr 
                  if( nta(l) .lt. nt-1 )then
c..   empty values, fill them by linear interpolation in logs
                     op1 = co(m,1,1,nta(l),l)
                     op2 = co(m,1,1,nt,l)
                     dlgt = alt(nt) - alt(nta(l))
                     dopdlgt = (op2-op1)/dlgt
                     do j = nta(l),nt
                        fact =  co(m,1,1,nta(l),l) 
     1                       + dopdlgt*(alt(j)-alt(nta(l)))
                        co(m,1,1,j,l) = fact
                     enddo
                  endif
               enddo
            endif
         enddo

         m = mf2
c..   11
         op11 = co(m,1,1,k3s,l3s)
c..   10
         op01 = co(m,1,1,k3s,l3s-1)
c..   01
         op10 = co(m,1,1,k3s-1,l3s)
c..   00
         op00 = co(m,1,1,k3s-1,l3s-1)

         yint(1) = alr(l3s-1)
         yint(2) = alr(l3s)
         xint(1) = alt(k3s-1)
         xint(2) = alt(k3s)
         zint(1,1) = op00
         zint(1,2) = op10
         zint(2,1) = op01
         zint(2,2) = op11

         call bile(slt,slr,xint,yint,zint,opact,dopact,dopacr)
         
         return

c      else
c         write(*,*)'OPAC: error: alpha enhanced, xxc, xxo ',xxc,xxo
c         write(*,*)'must add algorithm for missing opal table points'
c         stop'cofill from 67 in opac '

c      endif

c      stop'cofill'
      end
