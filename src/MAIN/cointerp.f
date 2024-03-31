      subroutine cointerp(xxc,xxo)
c     The purpose of this subroutine is to interpolate in C and O abund-
c     ances.
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

      save

      integer w
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb
     . ,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)

c      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),
c     . opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
ccccccccccccccccccccccc
      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),
     . opl(mx,nt,nr),cxdf(mx,mc),oxdf(mx,mc),itime(mx)

      common/aa/ q(4), h(4), xcd(mc),xod(mc), xc(mc),xo(mo)
     . ,xcs(mc),xos(mo), cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx)
     . ,nc,no

c      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr), index(101),
c     . t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr)
c     . ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) 
c     .,t6listf(ntm),opk2(nt,nr),dfsx(mx)

      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr),
     . t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr)
     . ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) 
     .,t6listf(ntm),opk2(nt,nr),dfsx(mx), index(101)


      common/bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq,xodp,xcdp,xxco,cxx,oxx
c-------------------------------------------------------------------
       is=0
      if(xxco .lt. 1.e-6) then
        do ir=l1,l1+iq
          do it=k1,k1+ip
            opl(m,it,ir)=co(m,1,1,it,ir)
          enddo
        enddo
            is=1
            go to 123
      endif
c     include boundaries that could later cause division by 0!
      if(xxc .gt. xcd(3)-1.d-6) then
c__________
      oxdp=log10(zzz+xodp)
c     handle possibility that xodp=0
      fac=max(min((oxx-ox(1))/max(oxdp-ox(1),1.d-6),1.0d0),0.0d0)
      do 40 ir=l1,l1+iq
      do 41 it=k1,k1+ip
c
c                    interpolation in region c1
c
c     include boundaries that could later cause division by 0!
      if(xxc .gt. xcd(2)-1.d-6) then
      iw=1
      a(1,m)=quad(is,iw,cxx,co(m,nc-2,1,it,ir),co(m,nc-1,1,it,ir),
     . diag(m,1,it,ir),cx(nc-2),cx(nc-1),cx(nc))
      iw=iw+1
      a(2,m)=quad(is,iw,cxx,diag(m,1,it,ir),diag(m,2,it,ir),
     . diag(m,3,it,ir),cxd(1),cxd(2),cxd(3))
         do w=1,2
           b(w)=a(w,m)
         enddo
c     handle possibility that xodp=0
           opl(m,it,ir)=b(1)+(b(2)-b(1))*fac
      is=1
      go to 41
      endif
c                    interpolation in region c2
c
      iw=1
      a(1,m)=quad(is,iw,cxx,co(m,nc-2,1,it,ir),co(m,nc-1,1,it,ir),
     . diag(m,1,it,ir),cx(nc-2),cx(nc-1),cx(nc))
      iw=iw+1
      a(2,m)=quad(is,iw,cxx,co(m,n(m,2)-2,2,it,ir),co(m,n(m,2)-1,2,it,
     . ir),diag(m,2,it,ir),cx(n(m,2)-2),cx(n(m,2)-1),cxd(2))
      iw=iw+1
      a(3,m)=quad(is,iw,cxx,diag(m,1,it,ir),diag(m,2,it,ir)
     .,diag(m,3,it,ir),cxd(1),cxd(2),cxd(3))
        do w=1,3
          b(w)=a(w,m)
        enddo
      iw=iw+1
      opl(m,it,ir)=quad(is,iw,oxx,b(1),b(2),b(3),ox(1),ox(2),oxdp)
      is=1
   41 continue
   40 continue
      if(is .eq. 1) go to 123
c__________
      endif
c
c                    interpolation in region c3 to c6
      is=0
c
      if(nc .ge. 5) then
c__________
      do 44 i=4,nc-1
c     do not go beyond middle (where c3-c6 overlaps o3-o6), and
        if((xxc .gt. xcd(i)-1.d-6) .and. (xxo .gt. xo(i-1)-1.d-6) .and. 
     $        (xcd(i-1) .gt. xc(i-1))) then
      do 42 ir=l1,l1+iq
      do 43 it=k1,k1+ip
        oxdp=log10(zzz+xodp)
        iw=1
        m1=i-1
        m2=i-2
        a(1,m)=quad(is,iw,cxx,co(m,n(m,m2)-2,m2,it,ir),co(m,n(m,m2)-1,m2
     x  ,it,ir),diag(m,m2,it,ir),cx(n(m,m2)-2),cx(n(m,m2)-1),cxd(m2))
        iw=iw+1
        a(2,m)=quad(is,iw,cxx,co(m,n(m,m1)-2,m1,it,ir),co(m,n(m,m1)-1,m1
     x  ,it,ir),diag(m,m1,it,ir),cx(n(m,m1)-2),cx(n(m,m1)-1),cxd(m1))
        iw=iw+1
        a(3,m)=quad(is,iw,cxx,diag(m,m2,it,ir),diag(m,m1,it,ir),
     x  diag(m,i,it,ir),cxd(m2),cxd(m1),cxd(i))
         do w=1,3
           b(w)=a(w,m)
         enddo
        iw=iw+1
        opl(m,it,ir)=quad(is,iw,oxx,b(1),b(2),b(3),ox(i-2),ox(i-1),oxdp)
        is=1
   43 continue
   42 continue
      if (is .eq. 1) go to 123
        endif
   44 continue
c__________
      endif
c
      if (is .eq. 1) go to 123
c
c     include boundaries that could later cause division by 0!
      if(xxo .gt. xod(3)-1.d-6) then
c__________
      cxdp=log10(zzz+xcdp)
c     handle possibility that xcdp=0
      fac=max(min((cxx-cx(1))/max(cxdp-cx(1),1.d-6),1.0d0),0.0d0)
      do 140 ir=l1,l1+iq
      do 141 it=k1,k1+ip
c
c                    interpolation in region  o1
c
c     include boundaries that could later cause division by 0!
      if(xxo .gt. xod(2)-1.d-6) then
      iw=1
      a(1,m)=quad(is,iw,oxx,co(m,1,no-2,it,ir),co(m,1,no-1,it,ir),
     . diago(m,no-1,it,ir),ox(no-2),ox(no-1),ox(no))
      iw=iw+1
      a(2,m)=quad(is,iw,oxx,diago(m,no-1,it,ir),diago(m,no-2,it,ir),
     .diago(m,no-3,it,ir),oxd(1),oxd(2),oxd(3))
        do w=1,2
          b(w)=a(w,m)
        enddo
c     handle possibility that xcdp=0
      opl(m,it,ir)=b(1)+(b(2)-b(1))*fac
      is=1
      go to 141
      endif
c                    interpolation in region  o2
c
      iw=1
      a(1,m)=quad(is,iw,oxx,co(m,1,no-2,it,ir),co(m,1,no-1,it,ir),
     . diago(m,no-1,it,ir),ox(no-2),ox(no-1),ox(no))
      iw=iw+1
      a(2,m)=quad(is,iw,oxx,co(m,2,n(m,2)-2,it,ir),co(m,2,n(m,2)-1,it,
     . ir),diago(m,no-2,it,ir),ox(n(m,2)-2),ox(n(m,2)-1),oxd(2))
      iw=iw+1
      a(3,m)=quad(is,iw,oxx,diago(m,no-1,it,ir),diago(m,no-2,it,ir),
     .diago(m,nc-3,it,ir),oxd(1),oxd(2),oxd(3))
        do w=1,3
          b(w)=a(w,m)
        enddo
      iw=iw+1
      opl(m,it,ir)=quad(is,iw,cxx,b(1),b(2),b(3),cx(1),cx(2),cxdp)
      is=1
  141 continue
  140 continue
      if(is .eq. 1) go to 123
c__________
      endif
c
c                    interpolation in region  o3 to o6
      is=0
      if(no .ge. 5) then
c__________
      do 144 i=4,no-1
c     do not go beyond middle (where o3-o6 overlaps c3-c6), and
      if((xxo .gt. xod(i)-1.d-6) .and. (xxc .gt. xc(i-1)-1.d-6) .and.
     $     (xod(i-1) .gt. xo(i-1)-1.d-6)) then
      do 142 ir=l1,l1+iq
      do 143 it=k1,k1+ip
      cxdp=log10(zzz+xcdp)
      iw=1
      m2=i-2
      m1=i-1
      a(1,m)=quad(is,iw,oxx,co(m,m2,n(m,m2)-2,it,ir),co(m,m2,n(m,m2)-
     .1,it,ir),diago(m,no-m2,it,ir),ox(n(m,m2)-2),ox(n(m,m2)-1),oxd(m2))
      iw=iw+1
      a(2,m)=quad(is,iw,oxx,co(m,m1,n(m,m1)-2,it,ir),co(m,m1,n(m,m1)-1, 
     . it,ir),diago(m,no-m1,it,ir),ox(n(m,m1)-2),ox(n(m,m1)-1),oxd(m1)) 
      iw=iw+1
      a(3,m)=quad(is,iw,oxx,diago(m,no-m2,it,ir),diago(m,no-m1,it,ir),
     .diago(m,no-i,it,ir),oxd(m2),oxd(m1),oxd(i))
        do w=1,3
          b(w)=a(w,m)
        enddo
      iw=iw+1
      opl(m,it,ir)=quad(is,iw,cxx,b(1),b(2),b(3),cx(m2),cx(m1),cxdp)
      is=1
  143 continue
  142 continue
      if (is .eq. 1) go to 123
      endif
  144 continue
c__________
      endif
c
      if (is .eq. 1) go to 123
c
c.....find index of C grid.
   52 ie=100*int(xxc)+1
      iei=index(ie)+1
c     must also allow index = nc, to avoid extrapolation
      if (iei .gt. nc) iei=nc
c
        if(iei .gt. 3) then
          i1=iei-2
          i2=iei-1
          i3=iei
        else
          i1=1
          i2=2
          i3=3
        endif
c
c.....find index of O grid
      ie=100*int(xxo)+1
      iej=index(ie)+1
c     must also allow index = no, to avoid extrapolation
      if (iej .gt. no) iej=no
c
        if(iej .gt. 3) then
          j1=iej-2
          j2=iej-1
          j3=iej
        else
          j1=1
          j2=2
          j3=3
        endif
c
c     lower-O part of grid: interpolate C before O
      if(j3.lt.no .and. i3.le.n(m,j3) .and.
     $       (xxc.lt.xcd(j3)+1.d-6 .or. xxc.ge.xxo))then
      do 20 ir=l1,l1+iq
      do 21 it=k1,k1+ip
      iw=0
        do jx=j1,j1+2
          iw=iw+1
c     if i3=n(m,jx), then must replace cx(i3) with cxd(jx)
          a(iw,m)=quad(is,iw,cxx,co(m,i1,jx,it,ir),co(m,i2,jx,it,ir),
     x    co(m,i3,jx,it,ir),cx(i1),cx(i2),min(cx(i3),cxd(jx)))
        enddo
        do w=1,3
          b(w)=a(w,m)
        enddo
      iw=iw+1
      opl(m,it,ir)=quad(is,iw,oxx,b(1),b(2),b(3),ox(j1),ox(j2),ox(j3))
      is=1
   21 continue
   20 continue
c     else: high-O part of grid: must interpolate O before C
      else
       do ir=l1,l1+iq
       do it=k1,k1+ip
        iw=0
        do ix=i1,i1+2
          iw=iw+1
          if(j3.lt.n(m,ix))then
            a(iw,m)=quad(is,iw,oxx,co(m,ix,j1,it,ir),co(m,ix,j2,it,ir), 
     $      co(m,ix,j3,it,ir),ox(j1),ox(j2),ox(j3))
          else
            a(iw,m)=quad(is,iw,oxx,co(m,ix,j1,it,ir),co(m,ix,j2,it,ir), 
     $      diago(m,no-ix,it,ir),ox(j1),ox(j2),oxd(ix))
          endif
        enddo
        do w=1,3
          b(w)=a(w,m)
        enddo
      iw=iw+1
      opl(m,it,ir)=quad(is,iw,cxx,b(1),b(2),b(3),cx(i1),cx(i2),cx(i3))
      is=1
       enddo
       enddo
      endif
  123 continue
      return
      end

c***********************************************************************
