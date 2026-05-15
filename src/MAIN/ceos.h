
c..ceos.inc include file

      integer*4 ndd, ntt, nerfd

      parameter( ndd=117, ntt=101, nerfd=29 )

      real*8 asi,aj2,aj3,aphi,bsi,bj2,bj3,bphi,
     1 dd, tt, array, ddd, ttd

      common/eosdata/asi(nerfd),aj2(nerfd),aj3(nerfd),aphi(nerfd),
     1   bsi(nerfd),bj2(nerfd),bj3(nerfd),bphi(nerfd),
     2   dd(ndd),tt(ntt),array(ndd,ntt,3), 
     3   ddd(4,ndd),ttd(4,ntt)

