      program program
      implicit none
      integer nmax, nzmax
      parameter (nmax = 30000, nzmax = 800000)
      integer i, ia(nmax+1), ja(nzmax)
      integer iout
      real*8 A(nzmax), rhs(1)
      integer job, ipos

      character title*72, key*8,type*3, guesol*2
      logical valued


      real*8 A1(nzmax)
      integer n, nrhs,nnz, ncol, nrow,ierr
      integer ia1(nmax+1), ja1(nzmax)

c Leer Matriz (CSR format)

      job = 2
      nrhs = 0

      open (unit=5,file='Matriz1.rsa')

      call readmt (nzmax,nzmax,job,5,A,ja,ia, rhs, nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)

c---- if not readable return 
      if (ierr .ne. 0) then
         write (iout,100) ierr
 100     format(' **ERROR: Unable to read matrix',/,
     *        ' Message returned fom readmt was ierr =',i3)
         stop
      endif

      CLOSE(5)

c Transformar de CSR => CSC
      write(*,*)'nrow',nrow
      n=nrow

      job = 1
      ipos = 1

      call csrcsc(n,job,ipos,A,ja,ia,A1,ja1,ia1)

      open (unit=7,file='matrixInfoCheck.out')
      title = 'Matrix4info'
      type  = 'RUA'
      key   = 'PR1'
      iout  = 7
      guesol='NN'
      valued = .TRUE.
      call dinfo1(n,iout,A1,ja1,ia1,valued,title,key,type,A,ja,ia)

      end program
