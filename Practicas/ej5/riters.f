      program riters
c-----------------------------------------------------------------------
c test program for iters -- the basic iterative solvers
c
c     this program generates a sparse matrix using
c     GEN57PT and then solves a linear system with an 
c     artificial rhs (the solution is a vector of (1,1,...,1)^T).
c-----------------------------------------------------------------------
c      implicit none
c      implicit real*8 (a-h,o-z)
      integer nmax, nzmax, maxits,lwk
      parameter (nmax=5000,nzmax=100000,maxits=60,lwk=nmax*40)

c Leer Matriz

      Integer ierr
      Integer nrow, ncol
      Integer job
      Integer n, nrhs,nnz

      Character title*72, key*8,type*3, guesol*2
      Logical valued

      real*8  A(nzmax)
      integer ia(nmax),ja(nzmax)

c Vectores x y b
      real*8 x(nmax),b(nmax)

c Temporal Matrix to process
      Real*8 AT(nzmax)
      Integer iat(nzmax), jat(nzmax)
c Almaceno D para calcular la matriz A completa
      Real diag(nzmax)
      Real*8 D(nzmax)
      Integer id(nzmax+1), jd(nzmax)
      Integer ioff(nzmax)
      Integer ndiag, idiag



c =========================
      integer jau(nzmax),ju(nzmax),iw(nmax*3)
      integer jr(nzmax),jwl(nzmax),jwu(nzmax)
      integer ipar(16),nx,ny,nz,i

      integer lfil,nwk
      
      real*8  sol(nmax),rhs(nmax),au(nzmax),wk(nmax*40)
      real*8  xran(nmax), fpar(16), al(nmax)
      real*8  gammax,gammay,alpha,tol
      external cg,bcg,dbcg
      external cgnr, fom, runrc, ilut
c     
      common /func/ gammax, gammay, alpha
c-----------------------------------------------------------------------  
c pde to be discretized is :
c---------------------------
c
c -Lap u + gammax exp (xy)delx u + gammay exp (-xy) dely u +alpha u
c
c where Lap = 2-D laplacean, delx = part. der. wrt x,
c dely = part. der. wrt y.
c gammax, gammay, and alpha are passed via the commun func.
c 
c-----------------------------------------------------------------------  
c
c data for PDE:
c
      nx = 6
      ny = 6
      nz = 1
      alpha = 0.0
      gammax = 0.0
      gammay = 0.0
c
c     set the parameters for the iterative solvers
c
      ipar(2) = 2
      ipar(3) = 1
      ipar(4) = lwk
      ipar(5) = 10 
      ipar(6) = maxits
      fpar(1) = 1.0D-5
      fpar(2) = 1.0D-10
c--------------------------------------------------------------
c call GEN57PT to generate matrix in compressed sparse row format
c
c     al(1:6) are used to store part of the boundary conditions
c     (see documentation on GEN57PT.)
c--------------------------------------------------------------
c      al(1) = 0.0
c      al(2) = 0.0
c      al(3) = 0.0
c      al(4) = 0.0
c      al(5) = 0.0
c      al(6) = 0.0
c      nrow = nx * ny * nz
c      call gen57pt(nx,ny,nz,al,0,nrow,a,ja,ia,ju,rhs)
c      print *, 'RITERS: generated a finite difference matrix'
c      print *, '        grid size = ', nx, ' X ', ny, ' X ', nz
c      print *, '        matrix size = ', nrow


c Leer Matriz (CSR format)

      job = 2
      nrhs = 0
      open (unit=5,file='./../ej4/Matriz1.rsa')

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

c Calculo la matríz completa
      job = 1
c Almaceno D
      idiag = 1
      ndiag=nrow
      call csrdia(nrow,idiag,job,A,ja,ia,ndiag,diag,ioff,D,jd,id,w)  
c A=A+A^T
      call apmbt(nrow,ncol,job,A,ja,ia,A,ja,ia,AT,jat,iat,nzmax,iw,ierr)
c A=A-D
      job = -1
      call apmbt(nrow,ncol,job,AT,jat,iat,D,jd,id,A,ja,ia,nzmax,w,ierr)

c Initialize x  
      do 1 j=1, n
          x(j) = real(j)
 1    continue

c Calculate b
      call amux(n,x,b, A, ja, ia)


c     set-up the preconditioner ILUT(15, 1E-4)  ! new definition of lfil
      lfil = 5
      tol = 1.0D-4
      nwk = nzmax
c     Bueno
c      call ilut(nrow,A,ja,ia,lfil,tol,au,jau,ju,nwk,wk,iw,ierr)

      call ilut(n,A,ja,ia,lfil,tol,au,jau,ju,nwk,wu,iw,jr,jwl,jwu,ierr)

      write(*,*)'ilut error: ',ierr
      write(*,*)'ipar values: ',ipar
      ipar(2) = 2
c     generate a linear system with known solution
      do i = 1, nrow
          sol(i) = 1.0D0
          xran(i) = 0.d0
      end do
      call amux(nrow, sol, rhs, A, ja, ia)

      print *, ' '
      print *, '	*** CG ***'
c      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,A,ja,ia,au,jau,ju,
c     +     cg)

10    call cg(n, rhs, sol, ipar, fpar, wk);
c      write(*,*)'ipar values: ',ipar
      if (ipar(1).eq.1) then
        write(*,*)'Requiere un producto matriz vector con A'
        call amux(n,wk(ipar(8)),wk(ipar(9)),A,ja,ia)
        goto 10
      else if (ipar(1).eq.2) then
        write(*,*)'Requiere un producto matriz vector con A^T'
        call atmux(n,wk(ipar(8)),wk(ipar(9)),A,ja,ia)
        goto 10
      else if (ipar(1).eq.3 .or. ipar(1).eq.5 ) then
        write(*,*)'Requiere un precondicionado por la izquierda Ml^-1'
        call lusol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
        goto 10
      else if (ipar(1).eq.4 .or. ipar(1).eq.6 ) then
        write(*,*)'Requiere un precondicionado por la izquierda Ml^T'
        call lutsol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
        goto 10
      else if (ipar(1).eq.10) then
c        call ‘mi propia rutina para hacer el test de parada’
          write(*,*)'mi propia rutina para hacer el test de parada'
          goto 10
      else if (ipar(1).gt.0) then
        write(*,*)'ipar(1) es un código no especificado'
        goto 10
      else
        write(*,*)'El proceso termina con ipar(1)=',ipar(1)
        write(*,*)'la solución es......'
      endif

!       print *, ' '
!       print *, '	*** BCG ***'
!       call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
!      +     bcg)
!       print *, ' '
!       print *, '	*** DBCG ***'
!       call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
!      +     dbcg)
!       print *, ' '
!       print *, '	*** CGNR ***'
!       call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
!      +     cgnr)
!       stop
      end
c-----end-of-main
c-----------------------------------------------------------------------
