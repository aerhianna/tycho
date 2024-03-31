      subroutine saver(jj,jsq,k,jdel,ja1,ja2,kmax,aa,qkh,qlh)
      implicit real*8(a-h,o-z)

      dimension qkh(jj),qlh(jsq),aa(ja1,ja2,kmax)

c..save recursion coeficients in array aa

      do 214 i = 1,jj
        ip = i + jj*jdel
        aa(ip,1,k) = qkh(i)
        do 214 j = 1,jj
          jp1 = j + 1
          ijk = jj*(j-1) + i
          aa(ip,jp1,k) = qlh(ijk)
  214 continue

      return
      end


