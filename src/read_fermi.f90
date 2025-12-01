        subroutine read_fermi
        use modmain
        implicit none

        open(unit=50,file='EFERMI.OUT',status='OLD')

        read(50,*)elkfermi

        close(50)

        end subroutine
