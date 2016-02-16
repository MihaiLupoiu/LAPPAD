      program program

      parameter ( nzmax=49776)

      integer i, ia(nzmax), ja(nzmax)
      real*8 A(nzmax)

      real*8 A1(nzmax)
      integer ia1(nzmax), ja1(nzmax)
      integer n, nz

      real*8 A2(nzmax)
      integer ia2(nzmax), ja2(nzmax)
      integer job, ipos

      character title*72, key*3,type*8, guesol*2
      logical valued
      integer iout

c Leer Matriz

      open (unit=5,file='L11_4_ringhals.txt')
      
      do 100 i = 1, 49776
      	read (5,*) ia(i), ja(i), A(i)
100   continue

C Imprimiendo Matriz
c      do 101 i = 0, 10
c      	PRINT *,	 ia(i), ja(i), A(i)
c 101   continue	

      n = 4680
      nz = 49776 

c Transformar de COO => CSR 

      call coocsr(n,nz,A,ia,ja,A1,ja1,ia1)

c      do 102 i = 0, 9
c      	PRINT *,     ia1(i), ja1(i), A1(i)
c 102   continue

c Transformar de CSR => CSC

      job = 1
      ipos = 1
	  
	call csrcsc(n,job,ipos,A1,ja1,ia1,A2,ja2,ia2)

	open (unit=7,file='matrixInfo.out')
      title = 'Matrix_info'
      type  = 'RUA'
      key   = 'PR1'
      iout  = 7
      guesol='NN'
      valued = .TRUE.

      call dinfo1(n,iout,A2,ja2,ia2,valued,title,key,type,A,ja,ia)

      CLOSE(7)

      end program