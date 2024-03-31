      subroutine haltwda(l,ll)

      implicit none
      integer*4 l,ll

      logical tobe
      character*8 lockfile

      data lockfile/'starlock'/

c unix style test to signal an orderly termination of run
c if lockfile exists, then halt

      inquire(file=lockfile,exist=tobe)
      if( tobe )then
        ll=l
        write(*,*)' starlock exists, terminate run at l = ', l
      endif
      return
      end

