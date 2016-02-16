      program program

      implicit none
      Integer nmax, nzmax
      parameter (nmax = 30000, nzmax = 800000)

      Integer ierr
      Integer nrow, ncol

      Character title*72, key*8,type*3, guesol*2
      Logical valued   

c Time
      real etime, t(2), t1, t2
c Pruebas
      Integer nTimes
      
c     Vector X and Y
      real*8 x(nmax),y(nmax)

c CSR
      Integer iout
      Integer i, j, ia(nzmax+1), ja(nzmax)
      Real*8 A(nzmax), rhs(1)
      Integer job, ipos
c Almaceno D para calcular la matriz A completa
      Real*8 D(nzmax)
      Integer id(nzmax+1), jd(nzmax)
      Integer ioff(nzmax)

c A=A+A^T
      Real iw(nmax), w(nzmax)

c Temporal Matrix to process
      Real*8 AT(nzmax)
      Integer n, nrhs,nnz
      Integer iat(nzmax), jat(nzmax)

c DIA
      Real diag(nzmax)
      Integer ndiag, idiag

C JAGGED
      Integer iperm(nzmax)


c Leer Matriz (CSR format)

      job = 2
      nrhs = 0
      open (unit=5,file='Matriz3.rua')

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

      write(*,*)'nrow',nrow
      write(*,*)'ncol',ncol
      n=nrow

c Initialize x  
      do 1 j=1, n
          x(j) = real(j)
 1    continue

c Transformar de CSR => CSC
      job = 1
      ipos = 1

      call csrcsc(n,job,ipos,A,ja,ia,AT,jat,iat)

      open (unit=7,file='matrix3Info.out')
      title = 'Matrix3info'
      type  = 'RUA'
      key   = 'PR1'
      iout  = 7
      guesol='NN'
      valued = .TRUE.
      call dinfo1(n,iout,AT,jat,iat,valued,title,key,type,A,ja,ia)
      

      nTimes = 3000

c Multiplicar Matriz *  Vecotr
      t1 = etime(t)         !  Startup etime - do not use result
      do 2 j=1, nTimes
            call amux(n,x,y, A, ja, ia)
 2    continue

      t2 = etime( t )
      print *, 'CSR: ', (t2-t1)/nTimes

c Tomar tiempo M*C CSC
      t1 = etime(t)         !  Startup etime - do not use result
      do 3 j=1, nTimes
            call atmux (n, x, y, AT, jat, iat)
 3    continue
      t2 = etime( t )
      print *, 'CSC: ', (t2-t1)/nTimes

c Transformar de CSR => MSR

      call csrmsr(n, A, ja, ia, AT, jat, At, jat)
      t1 = etime(t)         !  Startup etime - do not use result
      do 4 j=1, nTimes
            call amuxms (n, x, y, AT,jat)
 4    continue
      t2 = etime( t )
      print *, 'MSR: ',(t2-t1)/nTimes

c Transformar de CSR => DIA

      idiag = 30
      call csrdia(n,idiag,10,A,ja,ia,nmax,AT,ioff,AT,jat,iat,w) 
      t1 = etime(t)         !  Startup etime - do not use result
      do 5 j=1, nTimes
            call amuxd (n,x,y,AT,nmax,idiag,ioff) 
 5    continue     
      t2 = etime( t )
      print *, 'DIA: ', (t2-t1)/nTimes

c Transformar de CSR => JAGGED

      call csrjad (n, A, ja, ia, idiag, iperm, AT, jat, iat)
      t1 = etime(t)         !  Startup etime - do not use result
      do 6 j=1, nTimes
            call amuxj(n, x, y, idiag, AT, ja, iat)
 6    continue     
      t2 = etime( t )
      print *, 'JAD: ', (t2-t1)/nTimes


      end program