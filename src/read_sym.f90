        subroutine read_sym
        use modmain
        implicit none
        ! local variables
        integer isym,i

        ! reads the crystal symmetries
        open(50,file='SYMCRYS.OUT',status='old')
        
        do i=1,3
          read(50,*)
        end do 

        read(50,*) nsymcrys

        allocate(vtlsymc(3,nsymcrys))
        allocate(symlatspa(3,3,nsymcrys))
        allocate(symlatspin(3,3,nsymcrys))

        do isym=1,nsymcrys
          read(50,*)
          read(50,*)
          read(50,*)

          read(50,*) vtlsymc(:,isym)

          read(50,*)
          do i=1,3
            read(50,*) symlatspa(i,:,isym)
          end do
          
          read(50,*)
          do i=1,3
            read(50,*) symlatspin(i,:,isym)
          end do
        end do
        
        close(50)
        end subroutine
