        subroutine inv_lfe(n,a)
        use modmain

        ! inversion of epsgg'
        ! and redefinition of epsgg with the
        ! inclusion of local field effects
  
        implicit none
        ! arguments
        integer, intent(in) :: n
        complex(8), intent(inout) :: a(n,n)
        ! local variables
        integer info
        ! automatic arrays
        integer ipiv(n)
        complex(8) work(n)
       
       !subroutine zgetrf(m,n,a,lda,ipiv,info)	
       ! m,n integer rows and columns, a complex(8) input matrix
       ! lda leading dimension of a, ipiv pivot indices arrays
       ! info integer, =0 successful exit
       ! ZGETRF computes an LU factorization of a general M-by-N matrix A
       ! using partial pivoting with row interchanges.

        call zgetrf(n,n,a,n,ipiv,info)

        if (info /= 0) then
          write(*,*)
          write(*,'("Error(inv_lfe): unable to invert matrix")')
          write(*,'(" ZGETRF returned INFO = ",I8)') info
          write(*,*)
          stop
        end if

        !subroutine zgetri(n,a,lda,ipiv,work,lwork,info)
        ! n integer matrix order, a complex(8) input matrix
        ! from ZGETRF, lda integer leading dimension of the array,
        ! ipiv integer array dimension(n), work complex(8) array
        ! lwork integer dimension of the array work
        ! info integer, =0 successful exit
        ! ZGETRI computes the inverse of a matrix using the
        ! LU factorization computed by ZGETRF.

        call zgetri(n,a,n,ipiv,work,n,info)
        
        if (info /= 0) then
          write(*,*)
          write(*,'("Error(inv_lfe): unable to invert matrix")')
          write(*,'(" ZGETRI returned INFO = ",I8)') info
          write(*,*)
          stop
        end if

        end subroutine
