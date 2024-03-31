
      subroutine leqsnr(a,b,n,idim)
c..leqs look-alike using numerical recepies routines
      implicit none
      include 'dimenfile'
      integer idim, n
      real*8 a(ndim,ndim),b(ndim), foo
      integer indx(ndim)

      call ludcmpnr(a,n,ndim,indx,foo)
      call lubksbnr(a,n,ndim,indx,b)

      return
      end

      SUBROUTINE ludcmpnr(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX)
      d=1.D0
      do 12 i=1,n
        aamax=0.D0
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.0.D0) stop 'singular matrix in ludcmp'
        vv(i)=1.D0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue

        aamax=0.D0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue

        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif

        indx(j)=imax
        if(a(j,j).eq.0.D0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.D0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END

      SUBROUTINE lubksbnr(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.D0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END


      SUBROUTINE mprovenr(a,alud,n,np,indx,b,x)
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION a(np,np),alud(np,np),b(n),x(n)
      PARAMETER (NMAX=500)
CU    USES lubksb
      INTEGER i,j
      DOUBLE PRECISION r(NMAX)
      DOUBLE PRECISION sdp
      do 12 i=1,n
        sdp=-b(i)
        do 11 j=1,n
          sdp=sdp+dble(a(i,j))*dble(x(j))
11      continue
        r(i)=sdp
12    continue
      call lubksbnr(alud,n,np,indx,r)
      do 13 i=1,n
        x(i)=x(i)-r(i)
13    continue
      return
      END




