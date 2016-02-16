      program program
      parameter (nmax = 1000, nzmax=20*nmax)

      integer ia(nzmax), ja(nzmax), iwk(nmax)
      integer iao(nzmax), jao(nzmax)
      integer iao2(nzmax), jao2(nzmax)
      
      real*8 A(nzmax), AO(nzmax), AO2(nzmax) 
      character title*72, key*3,type*8, guesol*2
      logical valued
      integer ipos

      m = 100
      n = m
      ic = n/2
      index = 10
      alpha = 5.0
      nn = nzmax

      call matrf2(m,n,ic,index,alpha,nn,nz,A,ia,ja,ierr)
      print *, ierr

      call coocsr(m,nz,A,ia,ja,AO,jao,iao);

c Para funcionar dinfom hay que tener la 
c matriz en formato csc.

      job = 2
      ipos = 1

      call csrcsc(n,job,ipos,AO,jao,iao,AO2,jao2,iao2)

      open (unit=7,file='matrixInfo.out')
      title = 'Matrix_info'
      type  = 'RUA'
      key   = 'PR1'
      iout  = 7
      guesol='NN'
      valued = .TRUE.

      call dinfo1(n,iout,AO2,jao2,iao2,valued,title,key,type,A,ja,ia)

      CLOSE(7)
      end program
