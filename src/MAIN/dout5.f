      subroutine dout5(nc)

c..   writes model to file; allows debugging dumps

      implicit none
      include 'dimenfile'

      character*5 cmodel
      character*9 filout

      logical tobe

      include 'comod'
      include 'ctmstp'
      include 'cruntm.h'
      include 'czone'
      include 'cbug'
      include 'cgen'
      include 'cconst'

      integer*4 nd5, nc
      integer*4 ndumps
      parameter( ndumps=20 )
      character*2 ci(ndumps)

      data nd5/0/
      data ci/'00','01','02','03','04','05','06','07','08','09',
     1     '10','11','12','13','14','15','16','17','18','19'/

      save nd5, ci
c---------------------------------------------------------

      if( l .lt. 0 )then
         write(*,*)'dout5: l error: ',l
         stop'dout5'
      elseif( l .gt. 0 )then
c..   default: saves model on disk, in file named prefix//model

         if( l .lt. ll )then
            
            call modflgo(model,cmodel)
         else
            cmodel = '_last'
         endif

         if( nbin .eq. 0 )then
c..   this needs _last logic (see below)
cccccccccccccccccccccccccc
c..   binary output
            filout = '0'//prefix(1)//prefix(2)//cmodel

            call riteb(note,nc)

            return

         else
c..   formatted output
            filout = prefix(1)//prefix(2)//cmodel
            inquire(file=filout,exist=tobe)
            if( tobe .and. cmodel .ne. "_last" )then
c..   do not overwrite
               write(6,90)filout
               write(3,90)filout
 90            format(' file ',a9,' already exists, not overwritten')
               return
            elseif( cmodel .eq. "_last" )then
c..   last model in this run, ??_last file may exist
               write(*,'(a20,2a1,i5)')'last model is ',prefix(1),
     1              prefix(2),model
               open(8,file=filout,form='formatted',status='unknown')
            else
c..   new model
               open(8,file=filout,form='formatted',status='new')
            endif
         endif
         rewind 8
      else
c..   formatted output for debug
c..   debug dumps; place call dout5(0) where needed; only 10 allowed
         nd5 = nd5 + 1
         if( nd5 .gt. ndumps )then
            write(*,*)'dout5: too many debug dumps ',nd5
            stop'dout5'
         endif
         filout =  prefix(1)//prefix(2)//'dmp'//ci(nd5)
         write(*,*)'dump file ',filout,model,nd5
         open(8,file=filout,form='formatted')
         rewind 8
      endif

c
c      call nritef(note,nc)

      call ritef(note,nc)

      close(8)

      write(6,93)' file ',filout,' created, ',
     1     'to restart from this file: cp ',filout,' imodel'
      write(3,93)' file ',filout,' created, ',
     1     'to restart from this file: cp ',filout,' imodel'
 93   format(a6,a9,a10,a30,a7,a7)


      return

      end

