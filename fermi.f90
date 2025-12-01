        subroutine fermi(e,occ)
        use modmain 
        !local variables
        real (8), intent(in) :: e
        real (8), intent(out) :: occ
        real (8) x 


        x=(e-elkfermi)/smear

        if (x < -50.d0) then
          occ=2.d0
        else if (x > 50.d0) then
          occ=0.d0
        else
          occ=2.d0*((exp(x)+1.d0)**(-1.d0))
        end if 

        ! else if (abs(x).lt.1.d-14) then
        ! x=1.d-14*sign(1.d0,x)
        ! occ=2.d0*((exp(x)+1.d0)**(-1.d0))

        end subroutine
