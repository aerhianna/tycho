      subroutine mult(a,b,c,na,mb,n)
      implicit real*8(a-h,o-z)

      dimension a(na,n),b(n,mb),c(na,mb)

c..vector and matrix multiply a*b = c

      zero = 0.0d0

      do 10 i = 1,na
        do 10 j = 1,mb
          sum = zero
          do 11 k = 1, n
            fact = a(i,k)*b(k,j)
            sum = sum + fact
   11     continue
        c(i,j) = sum
   10 continue

      return
      end


