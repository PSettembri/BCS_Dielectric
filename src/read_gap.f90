        subroutine read_gap
        use modmain
        implicit none
        integer :: ios

        open(unit=50,file='GAP.OUT',status='OLD')

        ios=0
        do while (ios.eq.0)
          read(50,fmt=*,iostat=ios) elkgap
        end do 

        close(50)

        if (scissor.ne.0) then
          elkgap=elkgap+scissor
        end if

        end subroutine
