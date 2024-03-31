      subroutine seconds(runt)
c  timing routine for unix systems
      implicit none
      real*8 runt, cputime
c uses c real*8 function cputime which
c returns elapsed time in seconds
      runt = cputime()
c   runt is elapsed cpu time, in seconds
c   tested  by wda 6/9/98
      return
      end
