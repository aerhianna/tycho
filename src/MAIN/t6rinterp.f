      subroutine t6rinterp(slr,slt)
c     The purpose of this subroutine is to interpolate in logT6 and logR
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

      save
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


      common/d/dkap
      common/bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq,xodp,xcdp,xxco,cxx,oxx
      common/e/ opact,dopact,dopacr,dopactd
c------------------------------------------------------------------------
      is=0
      iu=0
      do kx=k1,k1+ip
          iw=1
        iu=iu+1
        h(iu)=quad(is,iw,slr,opk(kx,l1),opk(kx,l2),opk(kx,l3),
     x  alr(l1),alr(l2),alr(l3))
          if(iq. eq. 3) then
            iw=2
            q(iu)=quad(is,iw,slr,opk(kx,l2),opk(kx,l3),opk(kx,l4),
     x      alr(l2),alr(l3),alr(l4))
          endif
        is=1
      enddo

      is=0
      iw=1
c.....k and Dlog(k)/dlog(T6) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      opact=quad(is,iw,slt,h(1),h(2),h(3),alt(k1),alt(k2),alt(k3))
      dopact=dkap
      dkap1=dkap
        if (iq. eq. 3) then
c.....k and Dlog(k)/Dlog(T6) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2) 
          opactq=quad(is,iw,slt,q(1),q(2),q(3),alt(k1),alt(k2),alt(k3)) 
          dkapq1=dkap
        endif
        if(ip .eq. 3) then
c.....k and Dlog(k)/Dlog(T6) in lower-left 3x3.
          opact2=quad(is,iw,slt,h(2),h(3),h(4),alt(k2),alt(k3),alt(k4)) 
          dkap2=dkap
c.....k and Dlog(k)/Dlog(T6) smoothed in left 3x4
          dix=(alt(k3)-slt)*dfs(k3)
          dopact=dkap1*dix+dkap2*(1.d0-dix)
          opact=opact*dix+opact2*(1.d0-dix)
        if(iq .eq. 3) then

c.....k and Dlog(k)/Dlog(T6) in upper-right 3x3.
          opactq2=quad(is,iw,slt,q(2),q(3),q(4),alt(k2),alt(k3),alt(k4))
          dkapq2=dkap
          dopactq=dkapq1*dix+dkapq2*(1.d0-dix)
          opactq=opactq*dix+opactq2*(1.d0-dix)
         endif
        endif

      iu=0
      do lx=l1,l1+iq
        iw=1
        iu=iu+1
        h(iu)=quad(is,iw,slt,opk(k1,lx),opk(k2,lx),opk(k3,lx),
     x  alt(k1),alt(k2),alt(k3))
          if(ip .eq. 3) then
            iw=2
            q(iu)=quad(is,iw,slt,opk(k2,lx),opk(k3,lx),opk(k4,lx),
     x      alt(k2),alt(k3),alt(k4))
          endif
        is=1
      enddo

      is=0
      iw=1
c.....k and Dlog(k)/Dlog(R) in lower-left 3x3
      opacr=quad(is,iw,slr,h(1),h(2),h(3),alr(l1),alr(l2),alr(l3))
      dopacr=dkap
        if(ip .eq. 3) then
          opacrq=quad(is,iw,slr,q(1),q(2),q(3),alr(l1),alr(l2),alr(l3)) 
c.....k and Dlog(k)/Dlog(R) in upper-left 3x3.
          dopacrq=dkap
        endif
        if(iq .eq. 3) then
c.....k and Dlog(k)/Dlog(R) in lower-right 3x3.
          opact2=quad(is,iw,slr,h(2),h(3),h(4),alr(l2),alr(l3),alr(l4)) 
          dix2  =(alr(l3)-slr)*dfsr(l3)
          dopacr=dopacr*dix2+dkap*(1.d0-dix2)
c.....k and Dlog(k)/Dlog(T6) smoothed in both log(T6) and log(R)
              dopact=dopact*dix2+dopactq*(1.d0-dix2)
              opact =opact*dix2 +opactq*(1.d0-dix2)
         endif
        if(ip .eq. 3) then
         if(iq .eq. 3) then
c.....k and Dlog(k)/Dlog(R) in upper-right 3x3.
          opacrq=quad(is,iw,slr,q(2),q(3),q(4),alr(l2),alr(l3),alr(l4)) 
c.....Dlog(k)/Dlog(R) smoothed in both log(T6) and Log(R).
              dopacrq=dopacrq*dix2+dkap*(1.d0-dix2)
            endif
              dopacr=dopacr*dix+dopacrq*(1.d0-dix)
        endif
      dopactd = dopact - 3.0d0*dopacr
        if (opact .gt. 1.d+15) then
          write(*,'("Interpolation indices out of range",
     x              ";please report conditions.")')
          stop
        endif
      if (opact .gt. 9) then
      opact  =30.d0
      dopact =99.d0
      dopactr=99.d0
      dopactd=99.d0
      endif

c      if( slt .ge. 7.1d0-6.0d0 )then
c 6-6-02 testing
c         write(*,*)' '
c         write(*,'(1p8e12.3)')opact,opactq,dix,dix2
c         write(*,'(a5,12i3,1p8e12.3)')'t6r',ip,iq,is,iw,k1,k2,k3,k4,
c     1     l1,l2,l3,l4,
c     1     slt,opact,dopact,dopacr,dopactd,slr,slt
c
c         write(*,'(4i5,1p8e12.3)')ip,iq,is,iw,co(5,1,1,k4,l4),xa(5)
c
c         do jj=k1,k4
c            
c         write(*,'(f4.1,19f7.3)')alt(jj),(co(5,1,1,jj,j),j=1,19)
c         enddo
c         if( slt .ge. 7.35d0-6.0d0 )  stop't6rinterp'
c      endif

      return
      end

