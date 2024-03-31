      program hrplot

c..   reads hr.?? files from local directory
c..   pgplot 

      implicit none

      include 'dimenfile'
      include 'cconst'

      integer*4 it, l,model1, iradii
      integer*4 journal
      integer*4 pgopen, linesty
      integer*4 istat,lastmodel

      real*8 xej(ndim)
      real*8 time,dth,rphot,rhoph,uphot,xlol,telog,xmkk,
     1     peryear,speryear,omegkk
      real*8 radius,tlum,teff

      integer*4 ip,i,j, n, nfiles, numfiles, iobs
      integer*4 lw
      real*4    fsize
      parameter( ip = 2000, numfiles = 40)

      real*4       xx(ip),yy(ip)
      real*4       xmin,xmax,ymin,ymax
      character*5  hrfile(numfiles)
      character*80 toplab
      character*5  device
      character*72 text
      character*44 txt
      character*8  labl, labl1
      character*2 txtxt
      character*2 cprefix
      character*5 chr
ccccccccccccccc


      logical tobe,debug

      data l/0/, i/0/
      data lw/3/,fsize/1.2/
      data radius/1.0d0/,linesty/0/,iobs/0/

c..   special for jim liebert (8/12/05)
c..   if .true. this echos each header line to standard output
c..   to allow sm use if piped to a file
c      data debug /.true./
      data debug /.false./
      data txtxt /'  '/
c--------------------------------------------------------------------
      lastmodel = 0
c..   read graphics options

      inquire(file='hrplot.in',exist=tobe)
      if( tobe )then
         open(2,file='hrplot.in')
      else
         write(*,*)'hrplot: no file hrplot.in in directory'
         stop'hrplot: no hrplot.in input file error'
      endif

c..   dummy reads for labelling header lines

      read (2,*)text
      write(*,*)text
      read (2,*)text
      write(*,*)text
      read (2,*)text
      write(*,*)text

      read (2,*)txt,labl,device
      write(*,*)txt,txtxt,labl,device
      labl1 = 'device'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,journal
      write(*,*)txt,txtxt,labl,journal
      labl1 = 'journal'
      if( labl .ne. labl1 )goto 2000
      if( journal .ne. 0 )then
c..heavier line for journal
         lw = 5
      endif

      read (2,*)txt,labl,nfiles
      write(*,*)txt,txtxt,labl,nfiles
      labl1 = 'nfiles'
      if( labl .ne. labl1 )goto 2000

      if( nfiles .gt. numfiles )then
         write(*,*)'too many file:',nfiles,numfiles
         stop'hrplot'
      endif

      do n = 1 , nfiles
         read (2,*)txt,labl,hrfile(n)
         write(*,*)txt,txtxt,labl,hrfile(n)
         labl1 = 'hrfile'
         if( labl .ne. labl1 )goto 2000
      enddo

      read (2,*,err=2001)txt,labl,iradii
      write(*,*)txt,txtxt,labl,iradii
      labl1 = 'iradii'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xmin
      write(*,*)txt,txtxt,labl,xmin
      labl1 = 'xmin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,xmax
      write(*,*)txt,txtxt,labl,xmax
      labl1 = 'xmax'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ymin
      write(*,*)txt,txtxt,labl,ymin
      labl1 = 'ymin'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,ymax
      write(*,*)txt,txtxt,labl,ymax
      labl1 = 'ymax'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,linesty
      write(*,*)txt,txtxt,labl,linesty
      labl1 = 'linesty'
      if( labl .ne. labl1 )goto 2000

      read (2,*)txt,labl,iobs
      write(*,*)txt,txtxt,labl,iobs
      labl1 = 'iobs'
      if( labl .ne. labl1 )goto 2000

c..   dummy reads for labelling trailer lines
      read (2,*)text
      write(*,*)text

      close(2)


c      toplab = hrfile(1)
c      do n = 2, nfiles
c         toplab = hrfile(n)//' '//toplab
c      enddo
c      write(*,*)toplab
cccccccccccccccccccccccccccccccccccccccccccccc

c..   THE FOLLOWING VARIABLES WERE WRITTEN IN EDIT.F IN TYCHO
c..   model1   corresponding model number
c..   it       iterations used
c..   time     seconds elapsed in whole sequence
c..   dth(2)   time step
c..   rphot    photospheric radius
c..   rhoph    photospheric density (g/cc)
c..   uphot    photospheric velocity (cm/s)
c..   xlol     log( L / Lsun )              log => base 10
c..   telog    log( Teff )                  ln  => base e
c..   xm(kk)   mass (gm) interior to outer fitting radius r(1,kk)
c..   peryear  mass loss(-) or gain(+) in solar masses per year
c..   speryear solar masses ejected in this time step
c..   omeg(kk) angular rotational velocity (1/s)
c..   xej      ejected nucleon fraction array of size ndim

c..   initialize graphics

c..   scratch file for character conversion
      open(10,file='scratch.pg')

c      if ( pgbeg(0,device,1,1) .NE. 1 ) STOP' pgbeg error'
c..pgopen allows multiple pgplot windows
      if ( pgopen(device) .NE. 1 ) STOP' pgbeg error'

c..   no query for device
      call PGASK (.FALSE.)
c..   roman font = 2
c..   sans serif = 1
      call PGSCF(1)
c..   font scaling (size)
      call PGSCH(fsize)
c..   new page
      call PGPAGE

c..   window
      call PGSWIN(xmin,xmax,ymin,ymax)
c..   reset color
      call PGSCI(1)
c..   line width
      call PGSLW(lw)
c..   tickmarks on x(bottom=B,top=C) and y(left)
      call PGBOX('BCNST',0.0,0,'BCNST',0.0,0)

c     call PGBOX('A',0.0,0,'CMST',0.0,0)
c..   Left label, dispacement outside in character heights,
c..   location in fractions of edge, justification(0.5=centered)
c..   and character string

      call PGMTXT('L',2.0,0.5,0.5,'log L/sol')
      if( iradii .eq. 0 )then
         call PGMTXT('B',2.0,0.5,0.5,'log Te(K)')
      else
         call PGMTXT('B',2.0,0.5,0.5,'log R/sol')
      endif
c      call PGMTXT('T',2.0,0.5,0.5,toplab)

      do n = 1, nfiles

c..   l is running count of models read in this sequence
         l = 0
c..   reset color for this sequence
         if( n .le. 15 )then
c..   n=0 background
c..   n=1 white on black and vise versa
c..   up to n=15 light gray
            call PGSCI(n)
         elseif( n .le. 30 )then
c..   n=16 white on black and v.v.
c..   n=30 -->15 light gray
            call PGSCI(n-15)
         else
            write(*,*)'hrplot: too many colors ',n
         endif
c..   line width
         call PGSLW(lw)

c..   hr.? data files
         inquire(file=hrfile(n),exist=tobe)
         if( tobe )then
            open(8,file=hrfile(n))


 100        continue
            read(8,70,end=101,iostat=istat) model1, it, time, dth,
     1           rphot, rhoph, uphot, xlol, telog, xmkk, peryear,
     2           speryear, omegkk
     
            if( istat .ne. 0 )then
               write(*,*)'read error in hrplot.f '
               write(*,*)'file ',hrfile(n),' number ',n,' in line'
               write(*,*)'istat ',istat,' lastmodel ',lastmodel
               stop'hrplot.f'
            else
               lastmodel = model1
            endif

            if( debug )then
            write(*,70) model1, it, time, dth,
     1           rphot, rhoph, uphot, xlol, telog, xmkk, peryear,
     2           speryear, omegkk
            endif

c..   require self-consistency
            teff = 10.0d0**telog
            tlum = 10.0d0**xlol * sollum
            radius = sqrt( tlum /
     1           ( pi4*sigma*teff**4) )

            read(8,71,iostat=istat) xej
            if( istat .ne. 0 )then
               write(*,*)'2 read error in hrplot.f '
               write(*,*)'file ',hrfile(n),' number ',n,' in line'
               write(*,*)'istat ',istat,' lastmodel ',lastmodel
               stop'hrplot.f'
            endif

            l = l+1
c..   group ploted points in ip clumps for speed
            i = mod(l,ip)
            if( i .eq. 0 )then
               i     = ip
               if( iradii .eq. 0 )then
                  xx(i) = telog
               else
                  xx(i) = dlog10(radius/solrad)
               endif
               yy(i) = xlol
c..   draw curve
               if( linesty .eq. 0 .or. linesty .eq. 2)then
                  call PGPT(ip,xx,yy,20)
               elseif( linesty .eq. 1 .or. linesty .eq. 2 )then
                  call PGLINE(ip,xx,yy)
               endif
            else
c..   fill arrays
               if( iradii .eq. 0 )then
                  xx(i) = telog
               else
                  xx(i) = dlog10(radius/solrad)
               endif
               yy(i) = xlol
            endif
            goto 100
 101        continue
c..   fill in last values
            j     = i-1
            if( iradii .eq. 0 )then
               xx(j) = telog
            else
               xx(j) = dlog10(radius/solrad)
            endif
            yy(j) = xlol
c..   draw curve
            if( linesty .eq. 0 .or. linesty .eq. 2)then
               call PGPT(j,xx,yy,20)
            elseif( linesty .eq. 1  .or. linesty .eq. 2)then
               call PGLINE(j,xx,yy)
            endif
c..sequence finished; label
            if( journal .eq. 0)then
               xx(1) = xmin + float(n)*(xmax-xmin)*0.05
               yy(1) = ymax + 0.05*(ymax-ymin)
               chr = hrfile(n)
               cprefix = chr(4:5)
c               call pgptxt(xx(1),yy(1),0.0,0.5,cprefix)
            else
c..add labels for journal version
               if( n .eq. 4 )then
                  xx(1) = xmin + (xmax-xmin)*0.75
                  yy(1) = ymax - 0.15*(ymax-ymin)
                  call pgptxt(xx(1),yy(1),0.0,0.5,'1.643')
c..last label
                  call PGSCI(1)
                  xx(1) = xmin + (xmax-xmin)*0.65
                  yy(1) = ymax - 0.1*(ymax-ymin)
                  call pgptxt(xx(1),yy(1),0.0,0.5,'\\ga\\dML\\u=')

               elseif( n .eq. 3 )then
                  xx(1) = xmin + (xmax-xmin)*0.83
                  yy(1) = ymax - 0.1*(ymax-ymin)
                  call pgptxt(xx(1),yy(1),0.0,0.5,'2.323')
                  
               elseif( n .eq. 2 )then
                  xx(1) = xmin + (xmax-xmin)*0.9
                  yy(1) = ymax - 0.15*(ymax-ymin)
                  call pgptxt(xx(1),yy(1),0.0,0.5,'3.286')
                  
               elseif( n .eq. 1 )then
                  xx(1) = xmin + (xmax-xmin)*0.95
                  yy(1) = ymax - 0.1*(ymax-ymin)
                  call pgptxt(xx(1),yy(1),0.0,0.5,'5.190')
        
               else
                  write(*,*)'n is ',n
                  stop'hrplot: bad journal labels'
               endif
               
            endif

            close(8)

         else
            write(*,'(/a30)')'HRPLOT:  WARNING!'
            write(*,'(a10,i3,a15,a10,a15)')'number',n,' in list, 
     1           no file ',hrfile(n),' in directory'
            write(*,*)'Please generate this file or adjust hrplot.in'
         endif
      enddo

c..   reset color
      call PGSCI(1)

      if( iobs .ne. 0 )then
         if( iradii .eq. 0 )then
c..   luminosity and effective temperature  ......................
c..   font scaling (size)
            call PGSCH(fsize*0.8)
c..   present day sun
            call pgpt1(3.762,0.0,12)
            call pgptxt(3.762,0.0,0.0,0.0,' Sol')

            call PGSCI(2)
c..   andersen points for EM Car (L vs Te) 
c..   error in log Te of component A 22.3
            call pgerrx(1,4.505,4.557,5.02,1.0)
c..   error in log L of component A
            call pgerry(1,4.531,4.92,5.12,1.0)
c..   error in log Te of component B 20.3
            call pgerrx(1,4.505,4.561,4.92,1.0)
c..   error in log L of component B
            call pgerry(1,4.531,4.82,5.02,1.0)
c            call PGARRO(4.531,5.02,4.520,4.919)
c            call PGARRO(4.531,4.92,4.520,4.816)
            call pgptxt(4.55,5.02,0.0,1.0,'EM Car')

            call PGSCI(3)
c..   andersen points for V478 Cyg (L vs Te)
c..   error in log Te of component A 16.67
            call pgerr1(5,4.484,4.63,0.015,1.0)
c..   error in log L of component A
            call pgerr1(6,4.484,4.63,0.06,1.0)
c..   error in log Te of component B 16.31
            call pgerr1(5,4.485,4.63,0.015,1.0)
c..   error in log L of component A
            call pgerr1(6,4.485,4.63,0.06,1.0)
            call pgptxt(4.504,4.63,0.0,1.0,'V478Cyg')

            call PGSCI(4)
c..   andersen points for EK Cep (L vs Te)
c..   error in log Te of component A 2.03
            call pgerrx(1,3.944,3.964,1.17,1.0)
c..   error in log L of component A
            call pgerry(1,3.954,1.21,1.13,1.0)
c..   error in log Te of component B 1.124
            call pgerrx(1,3.741,3.771,0.21,1.0)
c..   error in log L of component B
            call pgerry(1,3.756,0.27,0.15,1.0)
            call pgptxt(3.97,1.17,0.0,1.0,'EK Cep')

            call PGSCI(5)
c..   andersen points for CW Cep (L vs Te)
c..   error in log Te of component A
            call pgerrx(1,4.467,4.437,4.27,1.0)
c..   error in log L of component A
            call pgerry(1,4.452,4.21,4.33,1.0)
c..   error in log Te of component B
            call pgerrx(1,4.458,4.426,4.15,1.0)
c..   error in log L of component B
            call pgerry(1,4.442,4.08,4.22,1.0)
            call pgptxt(4.47,4.27,0.0,1.0,'CW Cep')

            call PGSCI(6)
c..   andersen points for U Oph (L vs Te)
c..   error in log Te of component A
            call pgerr1(5,4.22,2.91,0.020,1.0)
c..   error in log L of component A
            call pgerr1(6,4.22,2.91,0.08,1.0)
c..   error in log Te of component B
            call pgerr1(5,4.119,2.70,0.02,1.0)
c..   error in log L of component B 
            call pgerr1(6,4.119,2.70,0.08,1.0)
            call pgptxt(4.24,2.91,0.0,1.0,'U Oph')


            call PGSCI(7)
c..   andersen points for TZ For (L vs Te)
c..   error in log Te of component A
            call pgerrx(1,3.690,3.708,1.59,1.0)
c..   error in log L of component A
            call pgerry(1,3.699,1.55,1.63,1.0)
c..   error in log Te of component B
            call pgerrx(1,3.796,3.810,1.36,1.0)
c..   error in log L of component B
            call pgerry(1,3.803,1.33,1.39,1.0)
            call pgptxt(3.68,1.59,0.0,0.0,'TZ For')

c..   reset color
            call PGSCI(8)
c..   andersen points for AI Hya (L vs Te)
c..   error in log Te of component A
            call pgerrx(1,3.822,3.830,1.44,1.0)
c..   error in log L of component A
            call pgerry(1,3.826,1.42,1.46,1.0)
c..   error in log Te of component B
            call pgerrx(1,3.847,3.855,1.24,1.0)
c..   error in log L of component B
            call pgerry(1,3.851,1.22,1.26,1.0)
            call pgptxt(3.83,1.44,0.0,1.0,'AI Hya')

            call PGSCI(9)
c..   ribas points for AI Hya (L vs Te)
c..   error in log Te of component A
            call pgerrx(1,3.842,3.860,1.54,1.0)
c..   error in log L of component A
            call pgerry(1,3.851,1.52,1.56,1.0)
c..   error in log Te of component B
            call pgerrx(1,3.860,3.878,1.31,1.0)
c..   error in log L of component B
            call pgerry(1,3.869,1.29,1.33,1.0)

c..   reset color
            call PGSCI(10)

c..   andersen points for zeta Phe (L vs Te)
c..   error in log Te of component A
            call pgerr1(5,4.163,2.51,0.010,1.0)
c..   error in log L of component A
            call pgerr1(6,4.163,2.51,0.04,1.0)
c..   error in log Te of component B
            call pgerr1(5,4.076,1.79,0.007,1.0)
c..   error in log L of component A
            call pgerr1(6,4.076,1.79,0.04,1.0)
            call pgptxt(4.183,2.51,0.0,1.0,'zeta Phe')

            call PGSCI(11)
c..   andersen points for IQ Per (L vs Te)
c..   error in log Te of component A
            call pgerr1(5,4.090,2.09,0.008,1.0)
c..   error in log L of component A
            call pgerr1(6,4.090,2.09,0.03,1.0)
c..   error in log Te of component B
            call pgerr1(5,3.885,0.85,0.008,1.0)
c..   error in log L of component A
            call pgerr1(6,3.885,0.85,0.04,1.0)
            call pgptxt(4.110,2.09,0.0,1.0,'IQ Per')

            call PGSCI(12)
c..   griffin points for HR6902 (V2291 Oph) 
c..   Schroeder, Pols, Eggleton (1997) MNRAS 285, 696 "best test"
c..   iwamoto and saio apj 521:297 1999
c..v2291 Oph  A=3.86+-0.15 Msun
c..   error in log Te of component A
            call pgerr1(5,3.686,2.73,0.009,1.0)
c..   error in log L of component A
            call pgerr1(6,3.686,2.73,0.08,1.0)
c..v2291 Oph  B= 2.95+-0.09 Msun
c..   error in log Te of component B
            call pgerr1(5,4.041,2.07,0.020,1.0)
c..   error in log L of component A
            call pgerr1(6,4.041,2.07,0.08,1.0)
            call pgptxt(3.66,2.73,0.0,0.0,'V2291 Oph')

            call PGSCI(13)
c..alpha Aur  A=2.69+-0.06 Msun
c..   error in log Te of component A
            call pgerr1(5,3.694,1.895,0.004,1.0)
c..   error in log L of component A
            call pgerr1(6,3.694,1.895,0.007,1.0)
c..alpha Aur  B= 2.56+-0.04 Msun
c..   error in log Te of component B
            call pgerr1(5,3.756,1.890,0.008,1.0)
c..   error in log L of component A
            call pgerr1(6,3.756,1.890,0.014,1.0)
            call pgptxt(3.674,1.895,0.0,0.0,'alpha Aur')

            call PGSCI(14)
c..zeta And  A=2.39+-0.14 Msun
c..   error in log Te of component A
            call pgerr1(5,3.703,1.81,0.007,1.0)
c..   error in log L of component A
            call pgerr1(6,3.703,1.81,0.02,1.0)
c..zeta And  B= 2.26+-0.14 Msun
c..   error in log Te of component B
            call pgerr1(5,3.700,1.60,0.007,1.0)
c..   error in log L of component A
            call pgerr1(6,3.700,1.60,0.02,1.0)
            call pgptxt(3.703,1.7,0.0,1.0,'zeta And')


            call PGSCI(2)
c..Procyon A: Te = 6543+- 84 K, MA=1.497 Msol, L=2.614e34 erg/s to 5 %
c..log Te = 3.816+5-6
c..Log L/Lsol = log(2.614e34/3.8515e33) = 0.8317 +0.021-0.022
c..   error in log Te of component A
            call pgerr1(5,3.816,0.8317,0.005,1.0)
c..   error in log L of component A
            call pgerr1(6,3.816,0.8317,0.022,1.0)
            call pgptxt(3.825,0.83,0.0,1.0,'Proc A')

c..   Sirius A: L/Lsol = 25.4+-1.3, log()=1.4048(1.4265,1.3820)+-0.02
c..   R/Rsol = 1.711+-0.013 
c..   M/Msol = 2.02+-0.03
c..   log(l/R**2)/4=0.2346
c..   log Te(sol) = log(5777 K) = 3.7617
c..   log Te = 3.9963
c..   error in log Te of component A
            call pgerr1(5,3.996,1.405,0.005,1.0)
c..   error in log L of component A
            call pgerr1(6,3.996,1.405,0.02,1.0)
            call pgptxt(4.01,1.405,0.0,1.0,'Sirius A')

c..   Tau Ceti: L/Lsol = 0.52+-?, log()=-0.284+-?
c..   R/Rsol = 0.793+-? 
c..   M/Msol = 0.84+-?
c..   log(l/R**2)/4=0.2346
c..   log Te = log(5375+-25 K) = 3.7304
c..   error in log Te of component A
            call pgerr1(5,3.7304,-0.284,0.005,1.0)
c..   error in log L of component A
            call pgerr1(6,3.7304,-0.284,0.02,1.0)
c            call pgptxt(4.01,1.405,0.0,1.0,'Sirius A')

c.. cepheids from  E. Schmidt (1984) ApJ 287, 261, table 2.....
            call PGSCI(1)
c..   EV sct NGC 6664
c..   error in log Te 
            call pgerr1(5,3.802,3.07,0.02,1.0)
c..   error in log L 
            call pgerr1(6,3.802,3.07,0.1,1.0)
c            call pgptxt(3.802,3.07,0.0,1.0,'6664')
c..   CE Cas b NBC 7790
c..   error in log Te 
            call pgerr1(5,3.771,3.06,0.02,1.0)
c..   error in log L 
            call pgerr1(6,3.771,3.06,0.1,1.0)
c            call pgptxt(3.771,3.06,0.0,1.0,'7790')
c..   CF Cas NBC 7790
c..   error in log Te 
            call pgerr1(5,3.747,2.99,0.02,1.0)
c..   error in log L 
            call pgerr1(6,3.747,2.99,0.1,1.0)
c            call pgptxt(3.747,2.99,0.0,1.0,'7790')
c..   CE Cas a NBC 7790
c..   error in log Te 
            call pgerr1(5,3.757,3.06,0.02,1.0)
c..   error in log L 
            call pgerr1(6,3.757,3.06,0.1,1.0)
c            call pgptxt(3.757,3.06,0.0,1.0,'7790')

c..   CV Mon anon
c..   error in log Te 
            call pgerr1(5,3.765,3.14,0.02,1.0)
c..   error in log L 
            call pgerr1(6,3.765,3.14,0.1,1.0)
c            call pgptxt(3.765,3.14,0.0,1.0,'cv mon')
c..   U Sgr M25
c..   error in log Te 
            call pgerr1(5,3.782,3.45,0.02,1.0)
c..   error in log L 
            call pgerr1(6,3.782,3.45,0.1,1.0)
c
c..   DL Cas NGC 129
c..   error in log Te 
            call pgerr1(5,3.774,3.50,0.02,1.0)
c..   error in log L 
            call pgerr1(6,3.774,3.50,0.1,1.0)
c
c..   S Nor NGC 6087
c..   error in log Te 
            call pgerr1(5,3.740,3.43,0.02,1.0)
c..   error in log L 
            call pgerr1(6,3.740,3.43,0.1,1.0)
c
c..   TW Nor Ly 6
c..   error in log Te 
            call pgerr1(5,3.732,3.20,0.02,1.0)
c..   error in log L 
            call pgerr1(6,3.732,3.20,0.1,1.0)
c


c..   reset color index
            call PGSCI(1)
c..   reset font scaling (size)
            call PGSCH(fsize)
         else
c..   luminosity and radius  .....................................
c..   andersen points for EM Car
            call pgerrx(1,0.9782,0.9626,5.02,1.0)
            call pgerry(1,0.9704,4.92,5.12,1.0)
            call pgerrx(1,0.927,0.913,4.92,1.0)
            call pgerry(1,0.9206,4.82,5.02,1.0)

            call PGARRO(0.9703,5.02,0.9438,4.919)
            call PGARRO(0.9207,4.92,0.8924,4.816)

c..   andersen points for EK Cep
c            call pgerrx(1,0.1965,0.2003,1.17,1.0)
c            call pgerry(1,0.1984,1.13,1.21,1.0)

c            call pgerrx(1,0.1169,0.1209,0.21,1.0)
c            call pgerry(1,0.1189,0.15,0.27,1.0)

c..   andersen points for AI Hya
c            call pgerrx(1,0.5892,0.5960,1.44,1.0)
c            call pgerry(1,0.5926,1.42,1.46,1.0)

c            call pgerrx(1,0.4392,0.4445,1.24,1.0)
c            call pgerry(1,0.4419,1.22,1.26,1.0)


c..   Tau Ceti: L/Lsol = 0.52+-?, log()=-0.284+-?
c..   R/Rsol = 0.793+-? logR = -0.100 
c..   M/Msol = 0.84+-?
c..   log(l/R**2)/4=0.2346
c..   log Te = log(5375+-25 K) = 3.7304
c..   error in R of component A
            call pgerr1(5,-0.1,-0.284,0.05,1.0)
c..   error in log L of component A
            call pgerr1(6,-0.1,-0.284, 0.02,1.0)
c            call pgptxt(4.01,1.405,0.0,1.0,'Sirius A')

c..   griffin points for HR6902 (V2291 Oph)
c..   error in log R/Rsun of component A
            call pgerr1(5,1.517,2.73,0.046,1.0)
c..   error in log L of component A
            call pgerr1(6,1.517,2.73,0.08,1.0)
c..   error in log R/Rsun of component B
            call pgerr1(5,0.477,2.07,0.067,1.0)
c..   error in log L of component A
            call pgerr1(6,0.477,2.07,0.08,1.0)

c..Procyon A: Te = 6543+- 84 K, MA=1.497 Msol, L=2.614e34 erg/s to 5 %
c..log Te = 3.816+5-6
c..log R = 0.3077 +0.002 -0.003
c..Log L/Lsol = log(2.614e34/3.8515e33) = 0.843 +0.021-0.022
c..   error in log R of component A
            call pgerr1(5,0.3077,0.843,0.001,1.0)
c..   error in log L of component A
            call pgerr1(6,0.3077,0.8317,0.028,1.0)
            call pgptxt(0.367,0.843,0.0,1.0,'Proc A')
         endif
      endif


c..temporary: dots for envelope debugging
c      open(34,file='fort.34')
c      n = 0
c      call PGSCI(1)
c 3001 continue
c      n = n+1
c      read(34,FMT='(2i5,1p8e12.3)',END=3000)i,j,
c     1                    xx(n),yy(n)
c      write(*,'(2i5,1p8e12.3)')i,j,
c     1                    xx(n),yy(n)
c
c      call PGPT(1,xx(n),yy(n),22)
c      goto 3001
c
c 3000 continue

c..   pause if on screen, else go on to clean-up
      if( device .eq. '/xwin' )then
         pause
      endif

      write(*,*)'normal termination'

      close(10)
      call PGEND
      stop

 70   format(2i6,1pe14.6,1p8e12.4,1p2e11.3)
 71   format(1p10e11.3)

 2000 continue

      write(*,'(//a40)')'hrplot: error in  input file: hrplot.in'
      stop'hrplot error'
 2001 continue
      write(*,'(//a60//)')
     1     'Does the nfiles value match the number of files given?'
      stop'hrplot files error'

      end


