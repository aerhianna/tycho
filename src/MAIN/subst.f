      subroutine subst(w,b,x,ipivot,n)
      implicit real*8(a-h,o-z)

c...............................................
c forward and backward substitution

c w is working matrix
c x gets solution of ax = b
c ipivot contains pivoting strategy

c conte and de boor, elem. num. analysis:
c      an algorithmic approach, p. 128-130
c...............................................

      dimension w(n,n),b(n),x(n),ipivot(n)

      zero = 0.0d0

      if( n .gt. 1 )  go to 10
      x(1) = b(1)/w(1,1)
      return

   10 ip = ipivot(1)
      x(1) = b(ip)
      do 15 k = 2,n
        ip = ipivot(k)
        km1 = k - 1
        sum = zero
        do 14 j = 1,km1
   14     sum = w(ip,j)*x(j) + sum
   15   x(k) = b(ip) - sum

      x(n) = x(n)/w(ip,n)
      k = n
      do 20 np1mk = 2,n
        kp1 = k
        k = k - 1
        ip = ipivot(k)
        sum = zero
        do 19 j = kp1,n
   19     sum = w(ip,j)*x(j) + sum
   20   x(k) = (x(k) - sum)/w(ip,k)

      return
      end


