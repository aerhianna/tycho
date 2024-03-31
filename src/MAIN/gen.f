c
c

      subroutine gen
      
c     ------------------------------------
c     GENERATE INITIAL MODEL FOR RUN:
c     READ AND WRITES FORMATTED ONLY
c     READS: new.model from:
c     unit 8, directory = local
c     ------------------------------------
c     ORIGINALY WRITTEN: 
c     7-29-06
c     REVISION HISTORY:
c     1-16-00
c     9-29-08 (C.Meakin, spring cleaning)
c     ------------------------------------

      
      implicit none
      
      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'ctmstp'
      include 'cburn'
      include 'cruntm.h'
      include 'cbug'
      include 'czone'
      include 'cgen'
      include 'cbound'
      include 'cpgplot'
      include 'cenv'
      include 'cconst'
      include 'cdtnuc'
      include 'cnabla'
      include 'crate'

      integer*4 j, k, nc, n, kedge

      real*8 dth1, sum, yesum, xmsol, zfak

      character*5  cmodel
      character*7  filehr, filecv
c..   allows prefix 0 for binary
      character*9  filein, filout
      character*40 label
      character*72 text
      character*44 txt
      character*8  labl, labl1
      character*2 txtxt

      character*7  cdummy
      character*25 fdate,cdate

      character*1 clause

      logical tobe

      data txtxt/'  '/

      
c     ----------------------------
c     FORMAT STATEMENT DEFINITIONS
c     ----------------------------
 6    format(1x,' redefined composition, leaving gen')
 7    format(a1)
 10   format(12i6)
 11   format(  4e15.7)
 12   format(i6,  4e15.7)
 13   format(i6,f10.3)
 15   format(8a8)
 21   format(2x,' ktot, l4, dlnv, xmmax, xmmin, dmmax, dmmaxc')
 22   format(2i6,5e13.5)
 23   format(3e15.7,i5)
 24   format(1h0,' later input models scratched')
 30   format(a4,1x,a5,a1)
 40   format(' ifm,ieos,nbug,nup')
 41   format(2x,' mode,  l1,  l2,  l3, ncond, iter, nfall, nshel,',
     .     ' modec, modes, modex, modez')
 50   format(1x,70('-'),1x)




c     --------------------
c     READ "params.d" FOR
c     CONTROL PARAMETERS
c     --------------------
      
      read (2,*)text
      write(*,*)text
      write(3,*)text
c..   dummy read
      read (2,*)text
      write(*,*)ctstamp
      write(3,*)ctstamp
c     
      read (2,*)text
      write(*,*)text
      write(3,*)text
      read (2,'(1x,a44,3x,a6,4x,2a1)')txt,labl,prefix
      write(*,*)txt,txtxt,labl,prefix
      write(3,*)txt,txtxt,labl,prefix
      labl1 = 'prefix'
      if( labl .ne. labl1 )goto 2000
c..   replace possible blanks with underscores in filenames
      if( prefix(1) .eq. ' ' )prefix(1)='_'
      if( prefix(2) .eq. ' ' )prefix(2)='_'
c     
      read (2,*)txt,labl,newnet
      write(*,*)txt,txtxt,labl,newnet
      write(3,*)txt,txtxt,labl,newnet
      labl1 = 'newnet'
      if( labl .ne. labl1 )goto 2000
c     
c..   warning for new network
      if( newnet .ne. 0 )then
         write(*,'(//a40//)')'******* new network ********'
      endif

      read (2,*)txt,labl,nbin
      write(*,*)txt,txtxt,labl,nbin
      write(3,*)txt,txtxt,labl,nbin
      labl1 = 'nbin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)text
      write(*,*)text
      write(3,*)text
      
      
      
c     -------------------------
c     CREATE FILENAMES AND OPEN 
c     -------------------------
c     file for HR diagram, mass loss data
      filehr = 'hr.'//prefix(1)//prefix(2)
      inquire(file=filehr,exist=tobe)
      if( tobe )then
         open(77,file=filehr,status='old')
      else
         open(77,file=filehr,status='new')
      endif
      rewind 77
      write(*,*)'HR file on unit 77 is ',filehr
      
c..   file for convection diagram
      filecv = 'cv.'//prefix(1)//prefix(2)
      open(11,file=filecv,status='unknown')
      rewind 11
      write(*,*)'CV file on unit 11 is ',filecv
      


c     ------------------------------------
c     CHECK FOR MODEL INPUT FILE: "imodel"
c     (formatted model data)
c     ------------------------------------
      filein = 'imodel'
      write(*,*)'Input model on unit 8 (filein) is ',filein
      
      inquire(file=filein,exist=tobe)
      if( tobe )then
         
c     call system('file imodel >dummy')
c     call system('file imodel')
         
c..   query to see what format convention is used for imodel file
c..   get first character (blank for formatted)
c..   unit 12 is opened and immediately closed here
         open(12,file=filein,status='OLD')
         rewind(12)
         read(12,'(a72)')note
         write(*,'(a72)')note
         close(12)
         clause = ' '
         clause = note(1:1)
         write(*,*)'header letter is ',note(1:1)         
         if( clause .eq. 'T' )then            
c..   formatted
            write(*,*)'reading formatted file ',filein,' on unit 8'
            open(8,file=filein,form='formatted',status='old')
            rewind 8
            nc = 1
            read(8,'(a72)')note
            write(*,*)note
            rewind 8
            if( note(7:7) .eq. '7' .or. note(7:7) .eq. '8' )then
c..   newer format (TYCHO 7.xx)
               write(*,*)'CALLING READF to read formatted file imodel'
               call readf(note,nc)
            else
c..   older format
               write(*,*)'CALLING NREAD to read formatted file imodelF'
               call nreadf(note,nc)
            endif
            
         elseif( clause .eq. 'H' )then
c     -----------------------------------
c     workaround: inquire not identifying 
c     unformatted files binary
c     -----------------------------------
            write(*,*)'reading binary file ',filein
            open(8,file=filein,form='unformatted',status='old')
            
            rewind 8
            nc = 1
            
c     ---------------------------------
c     READ IN MODEL FILE: (from unit=8)
c     ---------------------------------
            call readb(note,nc)
         else
            stop'gen format read error'
         endif
      else
         write(*,*)'gen: no file ',filein,' in directory'
         stop'gen'
      endif
      close(8)
c     ------------------------
c     CLOSE MODEL FILE: imodel
c     ------------------------



      
c     CHECK TIME STEP SETTING
c     
c     
      if( dth(2) .le. 0.0d0 )then
         write(*,'(a20,1p8e12.3)')'dth(2)',dth(2)
         stop'GEN: bad imodel time step'
      endif
      dth(1) = dth(2)
      dti(1) = dth(2)



c...........................................................
c..   reads revised parameters from unit params.d (default)
c     
c..   nbug.ne.0 gives diagnostic run
c..   ieos zero gives silent read of eos data file
c..   ieos negative gives write on 6 of eos data tables
c..   nup = 0 gives nbug = 1 on first time step
c..   nup = 1 stops following first failure
c     idt.ne.0 gives constant time step
c     newt = 1 gives newtonian gravity
c     nouter = 1 gives constant r(kk)
c--------------------------------------------------------
c     ktot .gt. 1 in zone gives normal logic for rezoning
c     otherwise, zone reads
c     iff = 0 add, 1 delete, -1 end
c     zq  = value in qqn array
c--------------------------------------------------------
c     ktot .eq. 0 in main gives skip of rezone call
c     ktot .lt. 0 in main following rezone call gives ktot = 0
c--------------------------------------------------------
c     ktot .lt. -1 here gives equal mass zoning (-ktot zones)
c--------------------------------------------------------
c     modez .eq. 0 gives old "zone division and merging"
c     modez .ne. 0 gives relaxation rezoning
c--------------------------------------------------------
c     ibnd = 0 gives orgin boundary condition for inner zone
c--------------------------------------------------------

      istep    = 0
      ncond    = 12
      nshel    = 0
      ibnd     = 0
      tenvelop = 0.0d0
      idt      = 0
      qcons    = 1.0d0
      etar     = 1.0d0
      temin    = 0.0d0

      
c     ------------------------
c     START READING "params.d"
c     ------------------------
      
      read (2,*)txt,labl,stime
      write(*,*)txt,txtxt,labl,stime
      write(3,*)txt,txtxt,labl,stime
      labl1   = 'stime'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ixstop
      write(*,*)txt,txtxt,labl,ixstop
      write(3,*)txt,txtxt,labl,ixstop
      labl1 = 'ixstop'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xstop
      write(*,*)txt,txtxt,labl,xstop
      write(3,*)txt,txtxt,labl,xstop
      labl1 = 'xstop'
      if( labl .ne. labl1 )goto 2000

      if( xstop .le. 0.0d0 .and. ixstop .ne. 0)then
         write(*,*)'error  gen: xstop =',xstop
         write(3,*)'error  gen: xstop =',xstop
         stop'gen: xstop'
      endif

      read (2,*)txt,labl,ll
      write(*,*)txt,txtxt,labl,ll
      write(3,*)txt,txtxt,labl,ll
      labl1 = 'll'

      if( labl .ne. labl1 )goto 2000
      
      if( stime .le. time .and. stime .gt. 0.0d0 .and. ll .ne. 1 )then
         write(*,*)'gen: stime error ',stime,time
         stop'gen: stime'
      endif

      if( newnet .ne. 0 )then
c..   reset number of time steps for a new network
         ll = 2
         write(*,'(/a40,i5/)')'Resetting ll for new network',ll
      endif
      
      read (2,*)txt,labl,ipause
      write(*,*)txt,txtxt,labl,ipause
      write(3,*)txt,txtxt,labl,ipause
      labl1 = 'ipause'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,runtmx
      write(*,*)txt,txtxt,labl,runtmx
      write(3,*)txt,txtxt,labl,runtmx
      labl1 = 'runtmx'
      if( labl .ne. labl1 )goto 2000
c     ' Cpu time in seconds for this run (0=ignore)' 'runtmx'    0.0

      read (2,*)text
      write(*,*)text
      write(3,*)text

      read (2,*)txt,labl,nsoltest
      write(*,*)txt,txtxt,labl,nsoltest
      write(3,*)txt,txtxt,labl,nsoltest
      labl1 = 'nsoltest'
      if( labl .ne. labl1 )goto 2000
c     ' Flag to freeze abundances(0=change)........' 'nsoltest'   0

      read (2,*)txt,labl,mset
      write(*,*)txt,txtxt,labl,mset
      write(3,*)txt,txtxt,labl,mset
      labl1 = 'mset'
      if( labl .ne. labl1 )goto 2000
c     ' Flag to freeze abundances(0=change)........' 'mset'      0

      read (2,*)txt,labl,nopac
      write(*,*)txt,txtxt,labl,nopac
      write(3,*)txt,txtxt,labl,nopac
      labl1 = 'nopac'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nkscale
      write(*,*)txt,txtxt,labl,nkscale
      write(3,*)txt,txtxt,labl,nkscale
      labl1 = 'nkscale'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,fthoul
      write(*,*)txt,txtxt,labl,fthoul
      write(3,*)txt,txtxt,labl,fthoul
      labl1 = 'fthoul'
      if( labl .ne. labl1 )goto 2000
      if( fthoul .lt. 0.0d0 .or. fthoul .gt. 1.0d1 )then
         write(*,*)'this is a large factor, fthoul = ',fthoul
         write(*,*)'if you really want this,'
         write(*,*)'revise circa line 244 in gen.f'
         stop'gen.f: fthoul'
      endif

      read (2,*)txt,labl,noion
      write(*,*)txt,txtxt,labl,noion
      write(3,*)txt,txtxt,labl,noion
      labl1 = 'noion'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,norad
      write(*,*)txt,txtxt,labl,norad
      write(3,*)txt,txtxt,labl,norad
      labl1 = 'norad'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nomole
      write(*,*)txt,txtxt,labl,nomole
      write(3,*)txt,txtxt,labl,nomole
      labl1 = 'nomole'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nocoul
      write(*,*)txt,txtxt,labl,nocoul
      write(3,*)txt,txtxt,labl,nocoul
      labl1 = 'nocoul'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nburn
      write(*,*)txt,txtxt,labl,nburn
      write(3,*)txt,txtxt,labl,nburn
      labl1 = 'nburn'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,neutro
      write(*,*)txt,txtxt,labl,neutro
      write(3,*)txt,txtxt,labl,neutro
      labl1 = 'neutro'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,jnb
      write(*,*)txt,txtxt,labl,jnb
      write(3,*)txt,txtxt,labl,jnb
      labl1 = 'jnb'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nopaleos
      write(*,*)txt,txtxt,labl,nopaleos
      write(3,*)txt,txtxt,labl,nopaleos
      labl1 = 'nopaleos'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,istep
      write(*,*)txt,txtxt,labl,istep
      write(3,*)txt,txtxt,labl,istep
      labl1 = 'istep'
      if( labl .ne. labl1 )goto 2000
      if( istep .lt. 0 .or. istep .gt. kk+1
     1     .or. istep .eq. 1 )then
         write(*,*)istep, ' istep error in gen'
      endif

      read (2,*)txt,labl,izams
      write(*,*)txt,txtxt,labl,izams
      write(3,*)txt,txtxt,labl,izams
      labl1 = 'izams'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,fmxenv
      write(*,*)txt,txtxt,labl,fmxenv
      write(3,*)txt,txtxt,labl,fmxenv
      labl1 = 'fmxenv'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,fmnenv
      write(*,*)txt,txtxt,labl,fmnenv
      write(3,*)txt,txtxt,labl,fmnenv
      labl1 = 'fmnenv'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,tmenv
      write(*,*)txt,txtxt,labl,tmenv
      write(3,*)txt,txtxt,labl,tmenv
      labl1 = 'tmenv'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,epsilon
      write(*,*)txt,txtxt,labl,epsilon
      write(3,*)txt,txtxt,labl,epsilon
      labl1 = 'epsilon'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,delchi
      write(*,*)txt,txtxt,labl,delchi
      write(3,*)txt,txtxt,labl,delchi
      labl1 = 'delchi'
      if( labl .ne. labl1 )goto 2000
c     ' network: fractional change.................' 'delchi'   0.2d0

      read (2,*)txt,labl,chimin
      write(*,*)txt,txtxt,labl,chimin
      write(3,*)txt,txtxt,labl,chimin
      labl1 = 'chimin'
      if( labl .ne. labl1 )goto 2000
c     ' network: minimum abundance used in dtime...' 'chimin'   1.0d-6

      read (2,*)txt,labl,fdtn
      write(*,*)txt,txtxt,labl,fdtn
      write(3,*)txt,txtxt,labl,fdtn
      labl1 = 'fdtn'
      if( labl .ne. labl1 )goto 2000
c     ' network: fractional increase in dtime......' 'fdtn'     1.1414d0

      read (2,*)txt,labl,fdysum
      write(*,*)txt,txtxt,labl,fdysum
      write(3,*)txt,txtxt,labl,fdysum
      labl1 = 'fdysum'
      if( labl .ne. labl1 )goto 2000
c     ' network: fractional of nucleons switching..' 'fdysum'   2.0d-2

      read (2,*)txt,labl,ncymax
      write(*,*)txt,txtxt,labl,ncymax
      write(3,*)txt,txtxt,labl,ncymax
      labl1 = 'ncymax'
      if( labl .ne. labl1 )goto 2000
c     ' network: subcycles allowed.................' 'ncymax'   200

      read (2,*)txt,labl,ncytest
      write(*,*)txt,txtxt,labl,ncytest
      write(3,*)txt,txtxt,labl,ncytest
      labl1 = 'ncytest'
      if( labl .ne. labl1 )goto 2000
c     ' network: write if subcycles .ge.this value.' 'ncytest'  195

      read (2,*)txt,labl,mrot
      write(*,*)txt,txtxt,labl,mrot
      write(3,*)txt,txtxt,labl,mrot
      labl1 = 'mrot'
      if( labl .ne. labl1 )goto 2000
c     ' rotation: no rotation=0, 1 and 2=yes.......' 'nrot'      0

      read (2,*)txt,labl,nsweep
      write(*,*)txt,txtxt,labl,nsweep
      write(3,*)txt,txtxt,labl,nsweep
      labl1 = 'nsweep'
      if( labl .ne. labl1 )goto 2000
c     ' mixing flag (cmix.f).......................' 'nsemi'     -1

      read (2,*)text
      write(*,*)text
      write(3,*)text

c..   debugging flag
      read (2,*)txt,labl,nbug
      write(*,*)txt,txtxt,labl,nbug
      write(3,*)txt,txtxt,labl,nbug
      labl1 = 'nbug'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nup
      write(*,*)txt,txtxt,labl,nup
      write(3,*)txt,txtxt,labl,nup
      labl1 = 'nup'
      if( labl .ne. labl1 )goto 2000
c     ' Maximum number of updates which fail.......' 'nup'        100

      read (2,*)txt,labl,newt
      write(*,*)txt,txtxt,labl,newt
      write(3,*)txt,txtxt,labl,newt
      labl1 = 'newt'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,nouter
      write(*,*)txt,txtxt,labl,nouter
      write(3,*)txt,txtxt,labl,nouter
      labl1 = 'nouter'
      if( labl .ne. labl1 )goto 2000
c     ' Constant outer radius (=1; ignore=0).......' 'nouter'    0

      read (2,*)txt,labl,mode
      write(*,*)txt,txtxt,labl,mode
      write(3,*)txt,txtxt,labl,mode
      labl1 = 'mode'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,l2
      write(*,*)txt,txtxt,labl,l2
      write(3,*)txt,txtxt,labl,l2
      labl1 = 'l2'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,l3
      write(*,*)txt,txtxt,labl,l3
      write(3,*)txt,txtxt,labl,l3
      labl1 = 'l3'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,iter
      write(*,*)txt,txtxt,labl,iter
      write(3,*)txt,txtxt,labl,iter
      labl1 = 'iter'
      if( labl .ne. labl1 )goto 2000
c     ' Number of iterations allowed...............' 'iter'      30

      read (2,*)txt,labl,modec
      write(*,*)txt,txtxt,labl,modec
      write(3,*)txt,txtxt,labl,modec
      labl1 = 'modec'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,alphaml
      write(*,*)txt,txtxt,labl,alphaml
      write(3,*)txt,txtxt,labl,alphaml
      labl1 = 'alphaml'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,uuml
      write(*,*)txt,txtxt,labl,uuml
      write(3,*)txt,txtxt,labl,uuml
      labl1 = 'uuml'
      if( labl .ne. labl1 )goto 2000

      if( uuml .le. 0.0d0 )then
         write(*,*)'error in uuml ',uuml
         goto 2000
      endif

c     overshoot parameter in pressure scale heights
      read (2,*)txt,labl,hmlfak
      write(*,*)txt,txtxt,labl,hmlfak
      write(3,*)txt,txtxt,labl,hmlfak
      labl1 = 'hmlfak'
      if( labl .ne. labl1 )goto 2000

c     switch for PAY mixing in radiative regions
      read (2,*)txt,labl,mixmode
      write(*,*)txt,txtxt,labl,mixmode
      write(3,*)txt,txtxt,labl,mixmode
      labl1 = 'mixmode'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,modes
      write(*,*)txt,txtxt,labl,modes
      write(3,*)txt,txtxt,labl,modes
      labl1 = 'modes'
      if( labl .ne. labl1 )goto 2000


      if( mode .eq. 0 .and. modes .eq. 2 )then
         write(*,'(/a60)')'GEN(ERROR) mode and modes are inconsistent'
         write(*,*)'GEN: mode = ',mode,' but modes = ',modes
         write(*,'(a60)')
     1 'hydrodynamic interior and hydrostatic envelope inconsistency'
         stop'GEN: mode, modes'
      endif

      if( modec .eq. -1 .and. modes .eq. 2 )then
         write(*,'(/a60)')'GEN(ERROR) modec and modes are inconsistent'
         write(*,*)'GEN: modec = ',modec,' but modes = ',modes
         write(*,'(a60)')
     1 'envelope requires convection be possible'
         stop'GEN: modec, modes'
      endif
      modex = 0
c     read (2,*)txt,labl,modex
c     write(*,*)txt,txtxt,labl,modex
c     write(3,*)txt,txtxt,labl,modex
c     labl1 = 'modex'
c     if( labl .ne. labl1 )goto 2000
c     ' Mixing mode(vc=0.1vs;vc=calculated)........' 'modex'     0

      read (2,*)txt,labl,modez
      write(*,*)txt,txtxt,labl,modez
      write(3,*)txt,txtxt,labl,modez
      labl1 = 'modez'
      if( labl .ne. labl1 )goto 2000
c     ' Zoning mode(0=binary;1=relax;2=smooth).....' 'modez'     0

      read (2,*)txt,labl,dth1
      write(*,*)txt,txtxt,labl,dth1
      write(3,*)txt,txtxt,labl,dth1
      labl1 = 'dth1'
      if( labl .ne. labl1 )goto 2000

c..   resid = 1e-6 works fine until RGB for 10 solar mass
c..   then 1e-5 is needed for convergence
      read (2,*)txt,labl,resid
      write(*,*)txt,txtxt,labl,resid
      write(3,*)txt,txtxt,labl,resid
      labl1 = 'resid'
      if( labl .ne. labl1 )goto 2000
c     ' iteration residual (fractional)............' 'resid'    1.0d-5 

c..   cdelt = 0.02 (and cdelv = 0.06) give trouble in pre-main sequence 
c..   contraction phase, cdelt = 0.01 works fine for M = one solar mass
      read (2,*)txt,labl,cdelt
      write(*,*)txt,txtxt,labl,cdelt
      write(3,*)txt,txtxt,labl,cdelt
      labl1 = 'cdelt'
      if( labl .ne. labl1 )goto 2000
c     ' fractional temperature change.....0.01.....' 'cdelt'     0.01

c..   cdelv=1.0 helps m3 in shell flashes ccccccccccccccccccc
      read (2,*)txt,labl,cdelv
      write(*,*)txt,txtxt,labl,cdelv
      write(3,*)txt,txtxt,labl,cdelv
      labl1 = 'cdelv'
      if( labl .ne. labl1 )goto 2000
c     ' fractional volume change..........0.03.....' 'cdelv'     0.03

c..   cdeln = 0.04 gives trouble in core H burning (ignition) in M6 (6 Msun)
      read (2,*)txt,labl,cdeln
      write(*,*)txt,txtxt,labl,cdeln
      write(3,*)txt,txtxt,labl,cdeln
      labl1 = 'cdeln'
      if( labl .ne. labl1 )goto 2000
c     ' fractional abundance change.......0.05.....' 'cdeln'     0.03

      read (2,*)txt,labl,ktot
      write(*,*)txt,txtxt,labl,ktot
      write(3,*)txt,txtxt,labl,ktot
      labl1 = 'ktot'
      if( labl .ne. labl1 )goto 2000

c..   force single zone to be added
      read (2,*)txt,labl,kforce
      write(*,*)txt,txtxt,labl,kforce
      write(3,*)txt,txtxt,labl,kforce
      labl1 = 'kforce'
      if( labl .ne. labl1 )goto 2000

c
      read (2,*)txt,labl,l4
      write(*,*)txt,txtxt,labl,l4
      write(3,*)txt,txtxt,labl,l4
      labl1 = 'l4'
      if( labl .ne. labl1 )goto 2000
      if( l4 .le. 0 )then
         write(*,*)'gen.f: l4 error, l4=',l4
         stop'GEN: L4'
      endif

      read (2,*)txt,labl,dlnv
      write(*,*)txt,txtxt,labl,dlnv
      write(3,*)txt,txtxt,labl,dlnv
      labl1 = 'dlnv'
      if( labl .ne. labl1 )goto 2000
c     ' d log rho/d zone.(zones per decade rho 0.2)' 'dlnv'      0.1

      read (2,*)txt,labl,xmmin
      write(*,*)txt,txtxt,labl,xmmin
      write(3,*)txt,txtxt,labl,xmmin
      labl1 = 'xmmin'
      if( labl .ne. labl1 )goto 2000
c     ' fractional mass in outer zone......1.0d-4..' 'xmmin'     2.0d-4 

      read (2,*)txt,labl,dmmax
      write(*,*)txt,txtxt,labl,dmmax
      write(3,*)txt,txtxt,labl,dmmax
      labl1 = 'dmmax'
      if( labl .ne. labl1 )goto 2000
c     ' biggest fractional zone mass...............' 'dmmax'     0.01d0 

      read (2,*)txt,labl,facte
      write(*,*)txt,txtxt,labl,facte
      write(3,*)txt,txtxt,labl,facte
      labl1 = 'facte'
      if( labl .ne. labl1 )goto 2000
c     ' rezone for shell (edge of flame s5/s5max)..' 'facte'     2.0e-2 

      read (2,*)txt,labl,vline
      write(*,*)txt,txtxt,labl,vline
      write(3,*)txt,txtxt,labl,vline
      labl1 = 'vline'
      if( labl .ne. labl1 )goto 2000
c     factor for v linearity rezone (gtintp.f)

      read (2,*)txt,labl,pline
      write(*,*)txt,txtxt,labl,pline
      write(3,*)txt,txtxt,labl,pline
      labl1 = 'pline'
      if( labl .ne. labl1 )goto 2000
c     factor for p linearity rezone (gtintp.f)

      read (2,*)txt,labl,drmax
      write(*,*)txt,txtxt,labl,drmax
      write(3,*)txt,txtxt,labl,drmax
      labl1 = 'drmax'
      if( labl .ne. labl1 )goto 2000
c     factor for radius increase in rezone (gtintp.f)


      read (2,*)txt,labl,ismoo
      write(*,*)txt,txtxt,labl,ismoo
      write(3,*)txt,txtxt,labl,ismoo
      labl1 = 'ismoo'
      if( labl .ne. labl1 )goto 2000
c     ' smooth zone masses in rezone(yes=1)........' 'ismoo'     0

      read (2,*)txt,labl,mapenv
      write(*,*)txt,txtxt,labl,mapenv
      write(3,*)txt,txtxt,labl,mapenv
      labl1 = 'mapenv'
      if( labl .ne. labl1 )goto 2000
      if( mapenv .ne. 0 )then
         ll = 2
         write(*,*)'GEN: mapenv = ',mapenv,' so force LL = ',ll
         if( modes .ne. 2 )then
           write(*,'(/a60)')
     1     'GEN(ERROR): trying to map envelope with no envelope'
           write(*,*)'modes = ',modes,' not equal to 2'
           write(*,*)'Change mapenv to zero (no mapping of envelope)'
           write(*,*)'or change modes = 2 (compute envelope)'
           stop'GEN: mapenv modes'
         endif
      endif
c     ' map envelope onto Henyey grid (0=no).......' 'mapenv'    0

      read (2,*)text
      write(*,*)text
      write(3,*)text

      read (2,*)txt,labl,peryear
      write(*,*)txt,txtxt,labl,peryear
      write(3,*)txt,txtxt,labl,peryear
      labl1 = 'peryear'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,mloss
      write(*,*)txt,txtxt,labl,mloss
      write(3,*)txt,txtxt,labl,mloss
      labl1 = 'mloss'
      if( labl .ne. labl1 )goto 2000

c..   read extra descriptor
      read (2,*)text
      write(*,*)text
      write(3,*)text

      read (2,*)text
      write(*,*)text
      write(3,*)text

c..   switch for Wolf-Rayet mass loss from Nugis & Lamers '02
c..   and quick & dirty pulsational mass loss 
      read (2,*)txt,labl,altloss
      write(*,*)txt,txtxt,labl,altloss
      write(3,*)txt,txtxt,labl,altloss
      labl1 = 'altloss'
      if( labl .ne. labl1 )goto 2000

c..   initialize for no mass loss
      if( mloss .eq. 0 )then
         peryear = 0.0d0
      endif

      read (2,*)text
      write(*,*)text
      write(3,*)text

c..   mass loss from blackowen
      read (2,*)txt,labl,bomloss
      write(*,*)txt,txtxt,labl,bomloss
      write(3,*)txt,txtxt,labl,bomloss
      labl1 = 'bomloss'
      if( labl .ne. labl1 )goto 2000

      read (2,*)text
      write(*,*)text
      write(3,*)text

      read (2,*)txt,labl,igraf
      write(*,*)txt,txtxt,labl,igraf
      write(3,*)txt,txtxt,labl,igraf
      labl1 = 'igraf'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,device1
      write(*,*)txt,txtxt,labl,device1
      write(3,*)txt,txtxt,labl,device1
      labl1 = 'device1'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ixflag
      write(*,*)txt,txtxt,labl,ixflag
      write(3,*)txt,txtxt,labl,ixflag
      labl1 = 'ixflag'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,nabflg
      write(*,*)txt,txtxt,labl,nabflg
      write(3,*)txt,txtxt,labl,nabflg
      labl1 = 'nabflg'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,nxflg
      write(*,*)txt,txtxt,labl,nxflg
      write(3,*)txt,txtxt,labl,nxflg
      labl1 = 'nxflg'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,fbot
      write(*,*)txt,txtxt,labl,fbot
      write(3,*)txt,txtxt,labl,fbot
      labl1 = 'fbot'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ftop
      write(*,*)txt,txtxt,labl,ftop
      write(3,*)txt,txtxt,labl,ftop
      labl1 = 'ftop'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xlogmin
      write(*,*)txt,txtxt,labl,xlogmin
      write(3,*)txt,txtxt,labl,xlogmin
      labl1 = 'xlogmin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xlogmax
      write(*,*)txt,txtxt,labl,xlogmax
      write(3,*)txt,txtxt,labl,xlogmax
      labl1 = 'xlogmax'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,cvsc
      write(*,*)txt,txtxt,labl,cvsc
      write(3,*)txt,txtxt,labl,cvsc
      labl1 = 'cvsc'
      if( labl .ne. labl1 )goto 2000

c..   HR plot

      read (2,*)txt,labl,device2
      write(*,*)txt,txtxt,labl,device2
      write(3,*)txt,txtxt,labl,device2
      labl1 = 'device2'
      if( labl .ne. labl1 )goto 2000

c..   plot limits for HR plot, used in hrplt.f
      read (2,*)txt,labl,gtmin
      write(*,*)txt,txtxt,labl,gtmin
      write(3,*)txt,txtxt,labl,gtmin
      labl1 = 'gtmin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,gtmax
      write(*,*)txt,txtxt,labl,gtmax
      write(3,*)txt,txtxt,labl,gtmax
      labl1 = 'gtmax'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,glmin
      write(*,*)txt,txtxt,labl,glmin
      write(3,*)txt,txtxt,labl,glmin
      labl1 = 'glmin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,glmax
      write(*,*)txt,txtxt,labl,glmax
      write(3,*)txt,txtxt,labl,glmax
      labl1 = 'glmax'
      if( labl .ne. labl1 )goto 2000

c..   cv (convection) plot

      read (2,*)txt,labl,device3
      write(*,*)txt,txtxt,labl,device3
      write(3,*)txt,txtxt,labl,device3
      labl1 = 'device3'
      if( labl .ne. labl1 )goto 2000

c..   plot limits for cv plot, used in cvplt.f
      read (2,*)txt,labl,pg3xmin
      write(*,*)txt,txtxt,labl,pg3xmin
      write(3,*)txt,txtxt,labl,pg3xmin
      labl1 = 'pg3xmin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,pg3xmax
      write(*,*)txt,txtxt,labl,pg3xmax
      write(3,*)txt,txtxt,labl,pg3xmax
      labl1 = 'pg3xmax'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,pg3ymin
      write(*,*)txt,txtxt,labl,pg3ymin
      write(3,*)txt,txtxt,labl,pg3ymin
      labl1 = 'pg3ymin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,pg3ymax
      write(*,*)txt,txtxt,labl,pg3ymax
      write(3,*)txt,txtxt,labl,pg3ymax
      labl1 = 'pg3ymax'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,ixpg3
      write(*,*)txt,txtxt,labl,ixpg3
      write(3,*)txt,txtxt,labl,ixpg3
      labl1 = 'ixpg3'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,iypg3
      write(*,*)txt,txtxt,labl,iypg3
      write(3,*)txt,txtxt,labl,iypg3
      labl1 = 'iypg3'
      if( labl .ne. labl1 )goto 2000


c..   db (iteration) plot

      read (2,*)txt,labl,device4
      write(*,*)txt,txtxt,labl,device4
      write(3,*)txt,txtxt,labl,device4
      labl1 = 'device4'
      if( labl .ne. labl1 )goto 2000

c..   plot limits for db plot, used in dbplt.f
      read (2,*)txt,labl,pg4xmin
      write(*,*)txt,txtxt,labl,pg4xmin
      write(3,*)txt,txtxt,labl,pg4xmin
      labl1 = 'pg4xmin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,pg4xmax
      write(*,*)txt,txtxt,labl,pg4xmax
      write(3,*)txt,txtxt,labl,pg4xmax
      labl1 = 'pg4xmax'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,pg4ymin
      write(*,*)txt,txtxt,labl,pg4ymin
      write(3,*)txt,txtxt,labl,pg4ymin
      labl1 = 'pg4ymin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,pg4ymax
      write(*,*)txt,txtxt,labl,pg4ymax
      write(3,*)txt,txtxt,labl,pg4ymax
      labl1 = 'pg4ymax'
      if( labl .ne. labl1 )goto 2000


      read (2,*)txt,labl,ixpg4
      write(*,*)txt,txtxt,labl,ixpg4
      write(3,*)txt,txtxt,labl,ixpg4
      labl1 = 'ixpg4'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,iypg4
      write(*,*)txt,txtxt,labl,iypg4
      write(3,*)txt,txtxt,labl,iypg4
      labl1 = 'iypg4'
      if( labl .ne. labl1 )goto 2000
      
c     ----------------------
c     END OF "params.d" FILE
c     ----------------------

 1002 continue



      
c     ----------------------
c     MORE TIME STEP STUFF
c     ----------------------
      if(dth1 .gt. 0.0) then
         dth(2) = dth1
         dth(1) = dth1
         dti(1) = dth1
      endif



c     ----------------------
c     LOG MODEL SETTINGS 
c     ---------------------
      write(3,50)
      write(6,50)
      
      if( mode .eq. 0 )then
         label = 'explicit 2nd-order hydrodynamics, mode ='
      elseif( mode .eq. 1 )then
         label = 'implicit damped hydrodynamics,    mode ='
      elseif( mode .eq. 2 )then
         label = 'implicit 1st-order hydrodynamics, mode ='
      elseif( mode .eq. 3 )then
         label = 'explicit thermal evolution,       mode ='
      elseif( mode .eq. 4 )then
         label = 'hydrostatic testing,              mode ='
      else
         write(*,*)mode
         stop'mode error'
      endif

      write(3,*)label,mode
      write(6,*)label,mode

      if( newt .eq. 1 )then
         label = 'Newtonian dynamics,              newt ='
      else
         label = 'GTR dynamics,                    newt ='
      endif
      write(3,*)label,newt
      write(6,*)label,newt

      if( modec .eq. 0 )then
         label = 'convection: dE + PdV = 0,       modec ='
         stop'modec gen: work in progress'
         
      elseif( modec .eq. 1 )then
         label = 'convection: Schwarzschild,      modec ='
      elseif( modec .eq. 2 )then
         label = 'convection: Ledoux+mixing,      modec ='
      elseif( modec .eq. -1 )then
         label = 'No convection,                  modec ='
      elseif( modec .eq. 3 )then
         label = 'convection: Richardson mix,     modec ='
      else
         label = 'error:                          modec = '
         write(6,*)label,modec
         stop
      endif

      write(3,*)label,modec
      write(6,*)label,modec

      if( nopac .eq. 0 )then
         label = 'OPAL type 2 opacity tables,   nopac = '
      elseif( nopac .eq. -1 )then
         label = 'for thomson opacity,          nopac = '
      elseif( nopac .eq. 1 )then
         label = 'OPAL type 1 opacity tables,   nopac = '
      else
         label = 'error:                        nopac = '
         write(6,*)label,nopac
         stop
      endif
      write(3,*)label,nopac
      write(6,*)label,nopac

      if( modes .eq. 0 )then
         label = 'Photosphere: outer zone kk,     modes ='
      elseif( modes .eq. 1 )then
         label = 'Photosphere: zone ks .le. kk,   modes ='
      elseif( modes .eq. 2 )then
         label = 'Photosphere: zone ks at kk+1,   modes ='
      else
         label = 'error:                          modes = '
         write(*,*)label
         write(*,*)'gen: bad parameter'
         stop
      endif
      write(3,*)label,modes
      write(*,*)label,modes

      if( modex .eq. 0 )then
         label = 'Mixing: conv. vel. = 0.1*sound, modex ='
      elseif( modex .eq. 1 )then
         label = 'Mixing: conv. vel. calculated,  modex ='
      else
         label = 'error:                          modex ='
         write(3,*)label,modex
         stop
      endif
      write(3,*)label,modex
      write(*,*)label,modex

      if( modez .eq. 0 )then
         label = 'Zoning: division and merging,   modez ='
      elseif( modez .eq. 1 )then
         label = 'Zoning: relaxation rezoning,    modez ='
      elseif( modez .eq. 2 )then
         label = 'Zoning: smoothing rezone,       modez ='
      elseif( modez .eq. 3 )then
         label = 'Zoning: smoothing2 rezone,      modez ='
      elseif( modez .eq. 4 )then
         label = 'Zoning: doubling rezone,        modez ='
      else
         label = 'error:                          modez ='
         write(*,*)label,modez
         stop
      endif
      write(3,*)label,modez
      write(*,*)label,modez


      if( nouter .eq. 1 )then
         write(3,*)'constant outer radius, nouter =',nouter
         write(6,*)'constant outer radius, nouter =',nouter
      elseif( nouter .ne. 0 )then
         write(3,*)'error: nouter = ',nouter
         write(6,*)'error: nouter = ',nouter
         stop
      endif

      write(3,50)
      write(*,50)



c     ---------------------------------
c     READ SOLAR SYSTEM ABUNDANCE DATA
c     reset abundances for new network:
c     if(newnet .ne. 0)
c     ---------------------------------
      call abinit
      
      

c     ---------------------------------
c     CHECK MASS FRACTIONS SUM TO UNITY
c     FOR EACH ZONE
c     ---------------------------------
      call xcheck(kk,x)



      if( newnet .eq. 0 )then
c..   normal case, uses existing network
c..   test abundance consistency
         do k = 2, kk
            do j = 1, ndim-1
               if( x(j,k) .lt. 0.0d0 )then
                  write(*,*)'GEN: negative abundance in initial model'
                  write(*,'(a8,2i5,a5,1pe12.3)')
     .                 "zone",k,j,cnuc(j),x(j,k)
                  stop'gen'
               endif
            enddo
         enddo
         write(*,*)'GEN: no negative abundances'
c..   always renormalize abundances
         if( modes .eq. 2 )then
            kedge = kk+1
         else
            kedge = kk
         endif
         do k = 2, kedge
            sum = -1.0d0
            do j = 1, ndim-1
               sum = sum + x(j,k)
            enddo
            if( abs(sum) .gt. 1.0d-8 )then
c..   more error than truncation of format would imply
               write(*,'(i5,1pe12.3,a30)')k,sum,
     .              "renormalization in GEN"
            endif
            do j = 1, ndim-1
               x(j,k) = x(j,k)/( 1.0d0 + sum )
            enddo
c..   reset Ye
            yesum = 0.0d0
            sum   = -1.0d0
            do j = 1, ndim-1
               yesum = yesum + x(j,k) * dble( lz(j) )/xa(j)
               sum = sum + x(j,k)
            enddo
            x(ndim,k) = yesum
         enddo
      else
         write(*,*)'gen: skipping renormalization'
      endif



c     --------------------------
c     CHECK MASS FRACTIONS AGAIN
c     --------------------------
      call xcheck(kk,x)



c     ---------------------
c     SETUP MASS COORDINATE
c     ---------------------
      dmi(1)    = dmh(2)
      dmi(kk+1) = dmh(kk+1)
      do k = 2, kk
         xm(k)  = xm(k-1) + dmh(k)
         dmi(k) = 0.5d0*( dmh(k+1) + dmh(k) )
      enddo

      if( modes .eq. 2 )then
c..   for consistency with hstat and fitenv
         dmi(kk) = dmh(kk)
      endif
      xm(kk+1) = xm(kk) + dmh(kk+1)
c..   gravity and area
      if( r(1,1) .gt. 0.0d0 )then
         k = 1
         g(k) = grav * xm(k) / r(1,k)**2
         a(k) = pi4 * r(1,k)**2
      else
         g(1) = 0.0d0
         a(1) = 0.0d0
      endif
      do k = 2, kk
 1       g(k) = grav * xm(k) / r(1,k)**2
         a(k) = pi4 * r(1,k)**2
      enddo
      
      xmsol = xm(kk)/sol
      
      
      write(*,'(a10,2a1,0pf12.5,a21)')
     .     'SEQUENCE ',prefix,xmsol,' solar masses on grid'
      write(3,'(a10,2a1,0pf12.5,a21)')
     .     'SEQUENCE ',prefix,xmsol,' solar masses on grid'
      if( modes .eq. 2 )then
         write(*,'(a12,0pf12.5,9x,a12,0pf12.5)')
     .        'envelope:',dmh(kk+1)/sol,'total mass:',
     .        xmsol+dmh(kk+1)/sol
         write(3,'(a12,0pf12.5,9x,a12,0pf12.5)')
     .        'envelope:',dmh(kk+1)/sol,'total mass:',
     .        xmsol+dmh(kk+1)/sol
      endif
      


c     ------------------------
c..   REMOVE "starlock" FILE
c     ------------------------
      inquire(file='starlock',exist=tobe)
      if( tobe )then
         call rename('starlock','tmp')
         write(6,*)'starlock removed by gen'
      endif



c..   initialize "note", a header for model dump files
c..   for first model or undefined ctstamp
      if( model .le. 1 .or.  note(1:10) .eq. '            '
     1     .or. note(1:10) .ne. ctstamp
     1     .or.  prefix(1) .ne. note(12:12)
     2     .or.  prefix(2) .ne. note(13:13) )then
         write(*,*)'GEN: redefining (note): old header is'
         write(*,*)note
         
c..   convert to character data for note
c..   number of nuclei in network = nnuc
c..   total mass in solar units = xm(kk+1)/sol
c..   mixing length parameter in pressure scale heights = alphaml
c..   metallicity (relative to solar system = zpop/0.015d0
c..   construct note (character*72, do not add beyond this!)
         
         if( model .le. 1 )then
c..   construct new zpop, using new abundances
            zpop = 0.0d0
            zsol = 0.0d0
            do n = 1, nnuc
               if( lz(n) .ne. 1 .and. lz(n) .ne. 2 )then
c..   uses outer zone on grid to define metallicity (k=kk)
                  zpop = zpop + x(n,kk)
                  zsol = zsol + solarx(n)
               endif
            enddo
c..   construct new initial hydrogen
            zhyd = 0.0d0
            do n = 1, nnuc
               if( lz(n) .eq. 1 )then
                  zhyd = zhyd + x(n,kk)
               endif
            enddo
c..   save for reference
            zhyd0 = zhyd
            zpop0 = zpop
            do n = 1, ndim
               xxin(n) = x(n,kk)
            enddo
         else
c..   use existing values
c..   this is a comparison sequence for a different parameter
            zhyd = zhyd0
            zpop = zpop0
         endif

         write(*,'(4(a10,1pe15.5))')'zpop',zpop,
     .        'zhyd',zhyd,'zhe',1.0d0-zpop-zhyd,'zsol',zsol
         
         
         if( zsol .le. 0.0d0 )then
            write(*,*)'gen: solar metallicity error'
            stop'gen: zsol'
         endif
         
         
         write(*,*)'GEN: redefining note header to'
         open(10,file='scratch')
c..   revised to include mixmode
         write(10,
     .        '(a10,1x,2a1,0pf7.2,i4,1pe9.1,0pf5.2,1x,0pf4.3,i2,a28)')
     .        ctstamp,prefix,xm(kk+1)/sol,netsize,zpop,
     .        alphaml,zhyd,mixmode
         
         rewind 10
         read(10,'(a72)')note
         close(10)
         call system('rm scratch')
         
c..   add present date/time at end of header
         cdate = fdate()
         note(47:72)=cdate
         write(*,*)note
         cdummy = prefix(1)//prefix(2)//'?????'
         write(*,*)'You may wish to delete any old models ',cdummy
      else
         cdummy = note(27:33)
         open(10,file='scratch')
         write(10,'(a7)')cdummy
         rewind 10
         read(10,'(f7.1)')zfak
c..   construct new zsol
         zsol = 0.0d0
         do n = 1, nnuc
            if( lz(n) .ne. 1 .and. lz(n) .ne. 2 )then
               zsol = zsol + solarx(n)
            endif
         enddo
         zpop = zfak
         close(10)

         call system('rm scratch')
c     
         write(*,'(3a10,1pe11.3,2(a10,1pe11.3))')
     .        'cdummy',cdummy,'zpop',zpop,'zsol',zsol,
     .        'zpop/sol',zpop/zsol
      endif
      


c     --------------------
c     WRITE THE MODEL FILE
c     --------------------
      call modflgo(model,cmodel)
      filout = prefix(1)//prefix(2)//cmodel
      
      inquire(file=filout,exist=tobe)
      if( .not. tobe )then
         open(8,file=filout,form='formatted',status='new')
         rewind 8
         nc = 1
         call ritef(note,nc)
         close(8)
         
         write(3,92)filout,model
         write(6,92)filout,model
 92      format(1x,'file ',a11,' is written, containing model',i6)
      endif


c     SUCCESS
      return
      


      
c     ERROR MESSAGE AND STOP
 2000 continue
      write(*,*)'ERR(gen): error in file named ',filein
      stop'ERR(gen):'
      
      end

