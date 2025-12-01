        subroutine read_lat
        use modmain
        implicit none
        ! local variables        
        integer ios,i
        character skip

          open(50,file='LATTICE.OUT',status='OLD',form='FORMATTED',iostat=ios)

        if (ios /= 0) then
          write(*,*)
          write(*,'("Error(readinput): error opening LATTICE.OUT")')
          write(*,*)
        stop
        end if

        do i=1,5
          read(50,*)
        end do 
       
        read(50,*) skip,skip,skip,avec(:,1)
        read(50,*) skip,skip,skip,avec(:,2)
        read(50,*) skip,skip,skip,avec(:,3)
        
        do i=1,11
          read(50,*)
        end do

        read(50,*) skip,skip,skip,skip,omega
        
        do i=1,6
          read(50,*)
        end do

        read(50,*) skip,skip,skip,bvec(:,1)
        read(50,*) skip,skip,skip,bvec(:,2)
        read(50,*) skip,skip,skip,bvec(:,3)
                
        do i=1,11
          read(50,*)
        end do

        read(50,*) skip,skip,skip,skip,omegabz

        close(50)

        end subroutine 
