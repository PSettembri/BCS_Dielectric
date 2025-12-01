        subroutine read_gvec
        use modmain
        implicit none
        ! local variables
        integer ig,skip
       
        open(50,file='GVECRF.OUT',status='unknown')
        
        read(50,*)
        read(50,*)

        do ig=1,ngrf
           read(50,*) skip,ivg(1:3,ig)
        end do

        close(50)
        
        end subroutine
