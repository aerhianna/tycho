
      subroutine opac1(izi,mzin,xh,t6,r,lstep,kzone)

c..   modified wda 4-1-08
c..   added lstep and kzone for debugging
c..   prettyprint indenting
c..   go to 67 option revised to use Klein-Nishina out of table

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
c      implicit none

c      real*8 ak,akt,theta,xxo,xxc,dixr,quad,dfsxmx1,xxmx1,xamx1
c      real*8 slt,slr,xxx,r,ri,t6i,t6,xxi,xh,z,dopactd,dopacr
c      real*8 dopact,opact,zval,dkap,b,xxh,h,q,xx,a1,a,dfsx,zza
c      real*8 zz,y,x,alrf,dfsz,za,xa,t6listf,xzf,dfsr,dfs

c      integer*4 iw,it,ir,is,ntlimit,mfin,iadvance,k1in,kmin,k3s,l3s
c      integer*4 mf2,mi,mh,mg,imd,ihi,ilo,mxend,itime,izi,mzin,iop
c      integer*4 m1,iq,ip,k4,k3,k2,k1,l4,l3,l2,l1,mf,m,mzz,i,n
c      integer*4 itab

      real*8 tlow1,thigh1,tlow2,thigh2,ropcut,
     1       opact2, dlgt, fact
      parameter( tlow1 = 0.0057d0, thigh1 = 0.007d0,
     1 tlow2 = 90.0d0, thigh2 = 200.0d0 )

c..... The purpose of this subroutine is to interpolate log kappa
c      in in X, T6, R
c        izi=0 recalulate table indices to use; =1 keep previous
c        mzin=index of za(i) in block data. za(i) are metallicities
c        t6=T6=temperature in millions of degrees kelvin
c        r=R=density(g/cm**3)/T6**3
c..... to use opac insert common/e/ in the calling routine.
c      This common contains interpolated values for kappa and its
c      first derivities.
c

c      integer w
      integer*4 mx,mz,nrm,nrb,nre,nr
      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/aa1/ q(4),h(4),xxh

      common/a1/ mzz, xz(mx,mz,nt,nr),  
     . t6list(nt),alr(nr),n(mx),alt(nt),opk(nt,nr),opk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),dfsz(mz),a(3,mx),b(3),m,mf,xa(mx)
     . ,alrf(nrm),xzf(nt,nr),t6listf(ntm),za(mz)

      common/b1/ itab(mx,mz),nta(nr),x(mx,mz),y(mx,mz),
     . zz(mx,mz)
      common/d/dkap
cccccccccccccccccccc
      common/bb1/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
      common/ee1/ opl(mx,nt,nr),xx(mx),zza(mz)
      common/eee1/m1,zval
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

      save

      data (xa(i),i=1,mx-1)/0.0,0.1,0.2,0.35,0.5,.7,.8,.9,.95/
      data (za(i),i=1,mz)/.0,0.0001,.0003,.001,.002,.004,.01,.02,.03,
     x .04,.06,.08,.1/
      data (nta(i),i=1,nrm)/14*70,69,64,60,58,57/
      data (n(i),i=1,mx)/13,13,13,13,13,13,13,13,10,12/
c--------------------------------------------------------------------
      iop=1   ! provides smoothed interpolations; iop=0 gives no smoothing
      mzz=mzin
      z=za(mzz)

      if(nre .lt. 6) go to 65
c
      if((izi .eq. 0) .and. (z+xh-1.d-6 .gt. 1 )) go to 61
      if((izi .ne. 0) .and. (zval+xh-1.d-6 .gt. 1 )) go to 61
      xxh=xh
      xxi=xh
      t6i=t6
      ri=r
c
      xxx=log10(.005d0+xh)
      slt=log10(t6)
      slr=log10(r)
c
      if(itime .ne. 12345678) then
         itime=12345678
         do  i=1,mx
            xx(i)=log10(0.005d0+xa(i))
         enddo
c..... this is the first time through. Calculate the decadic
c      log of the perimeter points shifted by Z. m refers to
c      xa(m); the hydrogen table value.
 
c..... read the data files
         call readco1
         xamx1=xa(mx-1)
         xxmx1=xx(mx-1)
         dfsxmx1=dfsx(mx-1)
      endif
      mxend=mx
      xa(mx)=1.-z
      xa(mx-1)=xamx1
      xx(mx-1)=xxmx1
      dfsx(mx-1)=dfsxmx1
      if (xa(mx) .lt. xa(mx-1)) then
         mxend=mx-1
         xa(mxend)=xa(mx)
      endif
      if (xh .ge. 0.8 ) then
         xx(mxend)=log10 (0.005+xa(mxend))
         dfsx(mxend)=1./(xx(mxend)-xx(mxend-1))
      endif
c
c
c..... Determine log R and log T6 grid points to use in the
c      interpolation.
      if((slt .lt. alt(1)).or.(slt .gt. alt(nt))) go to 62
      if((slr .lt. alr (1)).or.(slr .gt. alr(nre))) go to 62
c
c
c
      if (izi .eq. 0) then      ! freeze table indices
         ilo=2
         ihi=mx
 8       if(ihi-ilo .gt. 1) then
            imd=(ihi+ilo)/2
            if(xh .le. xa(imd)+1.e-7) then
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
         if (xh .lt. 1.d-6) then
            mh=1
            mg=1
            mi=2
            mf2=1
         endif
         if((xh .le. xa(2)+1.d-7) .or. (xh .ge. xa(mx-2)-1.d-7)) mf2=mh
c
         ilo=2
         ihi=nre
 12      if(ihi-ilo .gt. 1) then
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
c
         ilo=2
         ihi=nt
 11      if(ihi-ilo .gt. 1) then
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
         l3s=l3+nrb-1
         k3s=k3+ntb-1
      endif

      kmin=0
      k1in=k1
      iadvance=0
      mfin=mf
      if ((mfin .eq. 1) .and. (xz(1,mzz,k1,l1) .gt. 9.)) then ! no data
         do i=1,6
            if (xz(1,mzz,i,l1) .gt. 9.0d0)  then
               if (xh .lt. .1) then
                  kmin=i+1
               else
                  
                  if (iadvance .eq. 0) then
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
            if ((xz(1,mzz,kmin,l1+1) .lt. 9.0d0) .and.
     x           ((slr+.01) .gt. alr(l1+1))) then
               l1=l1+1
               kmin=0
               k1=k1in
               do i=1,6
                  if (xz(1,mzz,i,l1) .gt. 9.) kmin=i+1
               enddo
               if ((kmin .ne. 0) .and. (k1in .lt. kmin)) k1=kmin
            endif
         endif
         if ((slt+.001d0) .lt. alt(k1)) then
            opact=30.
            dopact=99.
            dopacr=99.
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
      do i=14,18                ! allows jagged edge at high T,rho
c     if((l3s .gt. i) .and. (k3s .gt. nta(i+1))) go to 62
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if((l3s .gt. i) .and. (k3s .gt. nta(i+1))) go to 67
      enddo
      do 123 m=mf,mf2
         ip=3
         iq=3
         ntlimit=nta(l3s)
         if((k3. eq. ntlimit) .or. (iop .eq. 0)) then 
            ip=2
            iq=2
         endif
         if(t6 .le. t6list(2)+1.d-7) ip=2

         if((l3 .eq. nre) .or. (iop .eq. 0)) then 
            iq=2
            ip=2
         endif
         if ((l4 .le.nr) .and. (xz(m,mzz,k3,l4) .eq. .0d0)) iq=2
         if(slr .le. alr(2)+1.d-7) iq=2
c
         is=0
c
         do ir=l1,l1+iq
            do it=k1,k1+ip
               opl(m,it,ir)=xz(m,mzz,it,ir)
               is=1
            enddo
         enddo
 123  continue
      if((zz(mg,mzin) .ne. zz(mf,mzin)) .or.
     x     (zz(mh,mzin) .ne. zz(mf,mzin))) then
         write(*,'("Z does not match Z in GN93hz files you are"
     x        ," using")')
         stop
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      if(z .ne. zz(mf,mzin)) go to 66
c                  with return
c..revised wda 11-18-05
      if( abs(z-zz(mf,mzin)) .gt. 1.0d-6 )go to 66
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      is=0
      iw=1
      do 45 ir=l1,l1+iq
        do it=k1,k1+ip
           if (mf2 .eq. 1) then
              opk(it,ir)=opl(mf,it,ir)
              go to 46
           endif
           opk(it,ir)=quad(is,iw,xxx,opl(mf,it,ir),opl(mg,it,ir)
     x          ,opl(mh,it,ir),xx(mf),xx(mg),xx(mh))
           is=1
 46        continue

        enddo
 45   continue

      if (mi .eq. mf2) then     ! interpolate between quadratics
         is=0
         iw=1
         dixr=(xx(mh)-xxx)*dfsx(mh)
         do 47 ir=l1,l1+iq
            do it=k1,k1+ip
               opk2(it,ir)=quad(is,iw,xxx,opl(mg,it,ir),opl(mh,it,ir)
     x              ,opl(mi,it,ir),xx(mg),xx(mh),xx(mi))
               opk(it,ir)=opk(it,ir)*dixr+opk2(it,ir)*(1.-dixr)
               is=1

            enddo
 47      continue
c     interpolate X between two overlapping quadratics
      endif

      is=0
c
c..... completed H,Z interpolation. Now interpolate T6 and log R on a
c      4x4 grid. (log(T6(i)),i=i1,i1+3),log(R(j)),j=j1,j1+3)).Procedure
c      mixes overlapping quadratics to obtain smoothed derivatives.
c
c
      call t6rinterp1(slr,slt)

c<<<<<<< .mine
cccccccccccccccccccccccccccccccccccc
c      if( t6 .gt. 50.90e0 )
c     1     write(*,'(1p8e12.3)')t6,r,slr,slt,opact,dopact,dopacr,
c     2     10.0d0**opact     
ccccccccccccccccccccccccccccc
c      if( kzone .eq. 169 )then
c         write(*,'(a12,2i5,1p8e12.3)')'opac1-norm',lstep,kzone,
c     1        dlog10(r),dlog10(t6)+6.0d0,opact
c      endif

c=======
c>>>>>>> .r98
      return

   61 write(*,'(" Mass fractions exceed unity")')
c                  with a return
      stop'61'
   62 write(*,'("opac1: T6/LogR outside of table range",1p3e11.3)')
     1     t6,r,dlog10(t6)+6.0d0
      write(*,'(5(a10,1pe12.3))')'T6',t6,'ropal',r,'log T',
     1     dlog10(t6)+6.0d0,'log rho',dlog10(r)+3.0d0*dlog10(t6),
     2     'logR',dlog10(r)
      write(*,'(5(a10,1pe12.3))')'z',z,'xh',xh,'xxc',xxc,'xxo',xxo
      write(*,'("slt,alt(1),alt(nt),slr,alr(1),alr(nre),l3s,i,k3s,
     1 nta(i+1)",6e12.5,4i5)') slt,alt(1),alt(nt),slr,alr(1),
     2 alr(nre),l3s,i,k3s,nta(i+1)
c                  with a return 
      stop'62'
   64 write(*,'(" X not equal to zero: To run this case it
     .is necessary"/ "to recompile with parameter (mx=1)")')
      stop'64'
   65 write(*,'("Too few R values; NRE+1-NRB < 6")')
      return
   66 write(*,'(" Z does not match Z in codata* files you are",
     . " using, 66 opac1.f")')
      write(*,'(1p10e12.3)')z,zz(mf,mzin),z-zz(mf,mzin)
      write(*,'(1p10e12.3)')zz
      write(*,*)mf,mzin,izi
ccccccccccccccccccccccccccccc
 67   continue
c..   conduction region (wda 4-2-01)
c..   the entries (log T6, log ropal) are in the table range 
c..   but the tabular entries for log opacity are blank, so we
c..   approximate the opacity near thomson (klein-nishina) limit
c        write(*,*)'fill 67, out of range'
c       write(*,'(a12,1p9e12.3)')'67 fill',z,xh,xxc,xxo,t6,r,slt,slr,opact
c       write(*,'(10i10)')i,l3s,k3s,mf2
c       write(*,'(5(a10,1pe12.3))')'logT',dlog10(t6)+6.0d0,'ropal',r,
c     1     'slt',slt,'slr',slr,'xxx',xxx,'xh',xh,'z',z,'opact',opact
c       write(*,'(5(a10,i10))')'l3s',l3s,'i',i,'k3s',k3s,
c     1     'nta(i+1)',nta(i+1)

c..   keep opact estimate sane
c..   keep opacity above klein-nishina limit
c      theta = 2.2d-3*t6
c      ak2 = ( 0.39d0*(xh + 1.0d0 )*0.5d0
c     1     /(1.0d0 + theta) )
c      ak = ak2

c..   there is a discontinuity between opact from tables
c..   and log ak2, which gives recurring nonconvergence.
c...................................................................
c..   opact is not defined yet!
c.....................................................wda 4-1-08
c.....................................................wda 12-1-08
c..   acceptable values
c<<<<<<< .mine
c        opact = dlog10( ak )
c.. approximate high density radiative opacity
c.. conduction should dominate
c        write(*,'(a12,1pe12.3,5i5,1p8e12.3)')'opact',opact,
c     1    l3s,k3s,nta(i+1),nta(i),i,xz(1,6,nta(i+1),i)
c      write(*,'(10i10)')mf,mf2,l1,l1+iq,k1,k1+ip,iq,ip,mzz
c..   make a choice of representative opacity near boundary
c..   interpolation would be more accurate
c..   concuction is supposed (?) to dominate here
      opact = xz(mf,mzz,nta(i+1),i)
c      m = mf
c      do it = k1,k1+1
c      do ir = l1,l1+1
c      write(*,'(a12,4i5,1p8e12.3)')'OPACT ',m,mzz,it,ir,
c     1  alt(it),alr(ir),
c     1  xz(m,mzz,it-1,ir),xz(m,mzz,it,ir),xz(m,mzz,it+1,ir),
c     2  xz(m,mzz,it+2,ir)
c      enddo
c      enddo
c      write(*,*)opact
c      write(*,*)t6,r,slt,slr
c     write(*,'(1p10e12.3)')t6list
c     write(*,'(1p10e12.4)')alr
c     write(*,'(10i12)')nta
c      stop'opac1aaa'

ccccccccccc
c=======
c        opact = dlog10( ak )
c.  . approximate high density radiative opacity
c..   conduction should dominate
c..   make a choice of representative opacity near boundary
c..   interpolation would be more accurate
c..   concuction is supposed (?) to dominate here
        opact = xz(mf,mzz,nta(i+1),i)
ccccccccccc
c>>>>>>> .r98
         dopact=99.0d0
         dopacr=99.0d0
c.....................................................wda 1/12/08
c<<<<<<< .mine
c      if( lstep .ge. 39 .and. kzone .eq. 155 )then
c      if( kzone .ge. 165 .and. kzone .le. 200 )then
c         write(*,'(a12,2i5,1p8e12.3)')'opac1-67',lstep,kzone,
c     1        dlog10(r),6.0d0+dlog10(t6),ak1,ak2,ak,opact
c      endif
c=======
c>>>>>>> .r98
      return


      stop
      end

