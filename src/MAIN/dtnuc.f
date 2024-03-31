      subroutine dtnuc(dth0,dth,nflag,ncycle,k)

c..   time step choice for subcycling
c..   dth is the updated value
c..   dth0 is the previously used value
c..   nflag = 0, before leqs, b(n) = dY(n)/dt
c..   else        after leqs, b(n) = dY(n)

      implicit none

      include 'dimenfile'
      include 'crate'
      include 'comcsolve'
      include 'cdtnuc'

      real*8    dth0,fak,dth,ydsum,aaa,taug

      integer*4 n, nflag, ncycle,k

c      data delchi/0.5d0/, chimin/1.0d-8/, fdtn/3.0d0/,
c     1 fdysum/5.0d-2/
c-------------------------------------------------------------
c..DELCHI is maximum fractional change allowed in number density
c..CHIMIN is smallest y for delchi restriction
c..FDTN is the largest increase factor allowed for timestep dth
c..NUCLEU is the fastest changing nucleus in dth evaluation
c..TAUS is fdysum/sum(abs(Ai*dYi/dt))
c..FDYSUM is of order 0.01 or less, in practice
c..this keeps the matrix solution happy (dth < taus)
c-------------------------------------------------------------

c..limit fractional increase in time step: tau(1)
      tau(1) = fdtn * dth0
      nucleu = 0

c..limit fractional change in abundance:   tau(2)
      tau(2) = tau(1)
      taug   = tau(1)

      if( nflag .eq. 0 )then
c..   before leqs, so b(n) is dY(n)/dt
c..nuclei
         do  n = 1, itot-3
            if( y(n) .gt. chimin  )then
               fak = b(n)
               if( fak .ne. 0.0d0 )then
                  taug  =  dabs( y(n)/fak )*delchi
               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
            endif
         enddo
         do  n = itot-2, itot
c            if( y(n) .gt. chimin*chimin  )then
            if( y(n) .gt. chimin )then
c..ignore neutrons on explicit guess (it has large errors)
               fak = b(n)
               if( fak .ne. 0.0d0 )then
                  taug  =  dabs( y(n)/fak )*delchi
               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
            endif
         enddo

      else
c..   after leqs, so b(n) is dY(n)
c..   nuclei
         do  n = 1, itot-3
c..   guess at best time step
            if( y(n) .gt. chimin )then
               fak = b(n)
               if( fak .ne. 0.0d0 )then
c..delchi
                  taug  =  dabs( ( 1.0d-3 + y(n) )/fak )*delchi*dth0

               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
            elseif( y(n) .gt. chimin**2 )then
               fak = b(n)
               if( fak .ne. 0.0d0 )then
c..delchi => 1
                  taug  =  dabs( ( 1.0d-3 + y(n) )/fak )*dth0
               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
            endif
         enddo

c         if( k .eq. 465 )then
c            write(*,'(a5,3i5,1p8e12.3)')'a ',ncycle,k,nucleu,tau(2),taug
c         endif



c..   n,p,alphas
         do  n =  itot-2,itot
c..   guess at best time step
            if( y(n) .gt. chimin )then
               fak = b(n)
               if( fak .ne. 0.0d0 )then
                  taug  =  dabs( ( 1.0d-3 + y(n) )/fak )*delchi*dth0
               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
c            elseif( y(n) .gt. 1.0d-15 )then
            elseif( y(n) .gt. 1.0d-6 )then
               fak = b(n)
               if( fak .ne. 0.0d0 )then
c..delchi = 1 for very low abundances (try to avoid negative y)
                  taug  =  dabs( ( y(n) )/fak )*dth0
               else
                  taug = tau(1)
               endif
               if( nucleu .eq. 0 .or. taug .lt. tau(2) )then
                  tau(2) = taug
                  nucleu = n
               endif
            endif
         enddo

c         if( k .eq. 465 )then
c            write(*,'(a5,3i5,1p8e12.3)')'b ',ncycle,k,nucleu,tau(2),taug
c         endif
c         if( ncycle .ge. 10 )then
c            write(*,*)ncycle,k
c            stop'dtnuc'
c         endif
ccccccccccccccccccc


      endif

c..limit fractional change in identity of nucleons:   tau(3)
c..limit total fractional number change to fdysum per step
      ydsum = 0
      do n = 1,itot
        aaa   = dble( nz(n) + nn(n) )
        ydsum = ydsum + abs( b(n)   )*aaa
      enddo

      if( abs(ydsum) .ne. 0.0d0 )then
        taus = fdysum/ydsum*dth0
      else
        taus = fdtn*dth0
      endif
      tau(3) = taus

c..decide on new value
      dth = tau(1)
      ndth = 1
      if( tau(2) .lt. dth )then
         dth =  tau(2)
         ndth = 2
      endif

      return
      end

