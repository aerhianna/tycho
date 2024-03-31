      function quad(ic,i,x,y1,y2,y3,x1,x2,x3)
c..... this function performs a quadratic interpolation.
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

      save
      common/d/dkap
      common/coquad/ xx(3),yy(3),xx12(30),xx13(30),xx23(30),xx1sq(30)
     . ,xx1pxx2(30)
      xx(1)=x1
      xx(2)=x2
      xx(3)=x3
      yy(1)=y1
      yy(2)=y2
      yy(3)=y3
        if(ic .eq. 0) then
          xx12(i)=1.d0/(xx(1)-xx(2))
          xx13(i)=1.d0/(xx(1)-xx(3))
          xx23(i)=1.d0/(xx(2)-xx(3))
          xx1sq(i)=xx(1)*xx(1)
          xx1pxx2(i)=xx(1)+xx(2)
        endif
      c3=(yy(1)-yy(2))*xx12(i)
      c3=c3-(yy(2)-yy(3))*xx23(i)
      c3=c3*xx13(i)
      c2=( yy(1)-yy(2) )*xx12(i) - xx1pxx2(i)*c3
      c1=yy(1)-xx(1)*c2-xx1sq(i)*c3
      dkap=c2+(x+x)*c3
      quad=c1+x*(c2+x*c3)
      return
      end
c
c********************************************************************** 

      block data
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb
     . ,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)
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


c      common/b/ itab(mx,ntabs),nta(nrm),x(mx,ntabs),y(mx,ntabs),
c     . zz(mx,ntabs),xca(mx,ntabs),xoa(mx,ntabs)
ccccccccccccccccccc
      common/b/ x(mx,ntabs),y(mx,ntabs),zz(mx,ntabs),
     1     xca(mx,ntabs),xoa(mx,ntabs),itab(mx,ntabs),nta(nrm)

c      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),
c     . opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
ccccccccccccccccccccccccccccccc
      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),
     . opl(mx,nt,nr),cxdf(mx,mc),oxdf(mx,mc),itime(mx)


      common/recoin/ itimeco,mxzero
      data itime/mx*0/,itimeco/0/
      data ( index(i),i=1,101)/1,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,
     . 4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,
     . 6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
     . 7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7/
      data (xcs(i),i= 1,mc)/ 0.0d0,0.01d0,0.03d0,0.1d0,0.2d0,0.4d0,
     1 0.6d0,1.0d0/
      data (xos(i),i= 1,mc)/0.0d0,0.01d0,0.03d0,0.1d0,0.2d0,0.4d0,
     1 0.6d0,1.0d0/
      data (xa(i),i=1,5)/0.0d0,0.03d0,0.1d0,0.35d0,0.7d0/
      data (nta(i),i=1,nrm)/14*70,69,64,60,58,57/
      end

C*********************************************************************
