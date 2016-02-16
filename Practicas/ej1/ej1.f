      program program
      parameter (nmax = 1000, nzmax=20*nmax)

      integer ia(nzmax), ja(nzmax), iwk(nmax)
      integer iao(nzmax), jao(nzmax)
      
      real*8 A(nzmax), ao(nzmax)
      character title*72, key*3,type*8, guesol*2

      open (unit=7,file='ej1.mat')
      m = 100
      n = m
      ic = n/2
      index = 10
      alpha = 5.0
      nn = nzmax

      call matrf2(m,n,ic,index,alpha,nn,nz,A,ia,ja,ierr)
      print *, ierr
      title = 'Matrix from zlatev'
      type  = 'RUA'
      key   = ' ZLATEV1'
      iout  = 7
      guesol='NN'

      ifmt = 3
      job = 2

      call coocsr(m,nz,A,ia,ja,ao,jao,iao);

      call prtmt (n,n,ao,jao,iao,rhs,guesol,title,type,key,
     &  ifmt,job,iout)
      end program

