      subroutine dout3(nc,period)

c..   detailed output for a given model (time step)

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'cbug'
      include 'cburn'
      include 'cqdv.h'
      include 'cruntm.h'
      include 'cnabla'
      include 'cgtr.h'
      include 'ceoset'
      include 'cenv'
      include 'cconst'
      include 'crate'
      include 'csurface'

      real*8 xx(ndim,kdm)
      real*8 capn(kdm),d(kdm),d2(kdm)
      real*8 odeep(kdm),gdeep(kdm)
    
      real*8 scr0(kdm),scr1(kdm),scr2(kdm)

      integer*4 ionuc

      integer*4 k1,i,j,k,n,kmax,nc,jnuc
      real*8    xminwr,period,d1sum,d2sum,xmfk,fact,critical,thomson
      real*8    fakk,ofak,gfak,fak,delgm
      real*8    eperiod, eenvel
      real*8    fak1,fak2,fak3,fak4,fak5
c-------------------------------------------------------------
c..   stub

      return
      end





