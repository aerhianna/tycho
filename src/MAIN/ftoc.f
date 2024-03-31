      subroutine ftoc(xxx,axistag)

c..converts floating point (real*4) to character string
      implicit none

      real*4 xxx
      character*10 axistag

      real*4 axxx
c----------------------------------------------

      axxx = abs( xxx )

      rewind 10

      if(    axxx .ge. 99.9 )then
c..to avoid rounding up to 3 digits to the left of decimal
        write(10,12) xxx

      elseif( axxx .ge. 9.99 )then
c..adjusts number of significant figures after the decimal point
        if( xxx .ge. 0.0 )then
          write(10,14) xxx
        else
          write(10,15) xxx
        endif
      elseif( axxx .ge. 0.1  )then
        if( xxx .ge. 0.0 )then
          write(10,13) xxx
        else
          write(10,14) xxx
        endif
      elseif( axxx .eq. 0.0  )then
        write(10,15) xxx
      else
        write(10,12) xxx
      endif


      rewind 10
c..read back as character string

      read(10,11) axistag

      return

11    format(a10)
12    format(1pe9.2,1x)
13    format(f5.3,5x)
14    format(f5.2,5x)
15    format(f5.1,5x)

      end



