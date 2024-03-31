
      subroutine t6rinterp1(slr,slt)

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

c     The purpose of this subroutine is to interpolate in logT6 and logR
      save
      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/ee1/ opl(mx,nt,nr),xx(mx),zza(mz)
      common/aa1/ q(4),h(4),xxh

      common/a1/ mzz, xz(mx,mz,nt,nr),  
     . t6list(nt),alr(nr),n(mx),alt(nt),opk(nt,nr),opk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),dfsz(mz),a(3,mx),b(3),m,mf,xa(mx)
     . ,alrf(nrm),xzf(nt,nr),t6listf(ntm),za(mz)

      common/d/dkap
cccccccccccccccccccc
      common/bb1/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
      common/e/ opact,dopact,dopacr,dopactd
c------------------------------------------------------------
      iu=0
      is=0
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
c
      is=0
      iw=1
c..... k and Dlog(k)/dlog(T6) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      opact=quad(is,iw,slt,h(1),h(2),h(3),alt(k1),alt(k2),alt(k3))
      dopact=dkap
      dkap1=dkap

        if (iq. eq. 3) then
c.....    k and Dlog(k)/Dlog(T6) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
          opactq=quad(is,iw,slt,q(1),q(2),q(3),alt(k1),alt(k2),alt(k3))
          dkapq1=dkap
        endif
        if(ip .eq. 3) then
c.....    k and Dlog(k)/Dlog(T6) in lower-left 3x3.
          opact2=quad(is,iw,slt,h(2),h(3),h(4),alt(k2),alt(k3),alt(k4))
          dkap2=dkap
c.....    k and Dlog(k)/Dlog(T6) smoothed in left 3x4
          dix=(alt(k3)-slt)*dfs(k3)
          dopact=dkap1*dix+dkap2*(1.-dix)
          opact=opact*dix+opact2*(1.-dix)
        endif
        if(iq .eq. 3) then
c.....    k and Dlog(k)/Dlog(T6) in upper-right 3x3.
          opactq2=quad(is,iw,slt,q(2),q(3),q(4),alt(k2),alt(k3),alt(k4))
          dkapq2=dkap
          dopactq=dkapq1*dix+dkapq2*(1.-dix)
          opactq=opactq*dix+opactq2*(1.-dix)
        endif
c
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
c
      is=0
      iw=1
c..... k and Dlog(k)/Dlog(R) in lower-left 3x3
      opacr=quad(is,iw,slr,h(1),h(2),h(3),alr(l1),alr(l2),alr(l3))
      dopacr=dkap
        if(ip .eq. 3) then
          opacrq=quad(is,iw,slr,q(1),q(2),q(3),alr(l1),alr(l2),alr(l3))
c.....    k and Dlog(k)/Dlog(R) in upper-left 3x3.
          dkapq3=dkap
        endif
        if(iq .eq. 3) then
c.....    k and Dlog(k)/Dlog(R) in lower-right 3x3.
          opact2=quad(is,iw,slr,h(2),h(3),h(4),alr(l2),alr(l3),alr(l4))
          dix2=(alr(l3)-slr)*dfsr(l3)
          dopacr=dopacr*dix2+dkap*(1.-dix2)
            if(ip .eq. 3) then
c.....        k and Dlog(k)/Dlog(T6) smoothed in both log(T6) and log(R)
              dopact=dopact*dix2+dopactq*(1.-dix2)
              opact=opact*dix2+opactq*(1-dix2)
            endif
         endif
        if(ip .eq. 3) then
c.....    k and Dlog(k)/Dlog(R) in upper-right 3x3.
          opacrq=quad(is,iw,slr,q(2),q(3),q(4),alr(l2),alr(l3),alr(l4))
            if(iq .eq. 3) then
c.....        Dlog(k)/Dlog(R) smoothed in both log(T6) and Log(R).
              dopacrq=dkapq3*dix2+dkap*(1.-dix2)
              dopacr=dopacr*dix+dopacrq*(1.-dix)
            endif
        endif
      dopactd=dopact-3.*dopacr
        if (opact .gt. 1.e+15) then
          write(*,'("Interpolation indices out of range",
     x              ";please report conditions.")') 
          stop
        endif
          if (opact .gt. 9.) then
            dopact=99.
            dopacr=99.
            dopactd=99.
          endif
      return
      end

