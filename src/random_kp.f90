        subroutine random_kp(ene_grid,vel_grid,nkgoal,k_rec,k_cart, & 
                                & energy,velocity,weight)
        ! for a given band this routine selects the random points according to a Monte Carlo
        ! distribution and interpolates the energy and velocity from a FFT grid to the random points.
        use modmain
        use modrand
        implicit none
        integer :: iran,ntotran,nkgoal
        real(8) :: energy(nkgoal),velocity(nkgoal,3),weight(nkgoal) 
        real(8) :: k(3),ene,xxx,Prbl_Model,Vk(3)
        real    :: ene_grid(inp%nr(1),inp%nr(2),inp%nr(3))
        real    :: vel_grid(inp%nr(1),inp%nr(2),inp%nr(3),3)
        real(8) :: k_rec(nkgoal,3),k_cart(nkgoal,3)
        real(8) :: weightsum,kc(3)
        integer :: nk
        
        nk=0
        weightsum=0.d0
        ntotran=0

        do iran=1,10**8 
           ntotran=ntotran+1
           call random_number(k(:))
           call get_eneofk(k,ene,ene_grid,vel_grid,Vk(:),kc)
           
           call random_number(xxx)
           ! Pruning by a Metropolis algorithm with a smooth distribution
           if (xxx.gt.Prbl_model(ene)) cycle 
           nk=nk+1
           k_rec(nk,:)=k(:)
           k_cart(nk,:)=kc(:)
           energy(nk)=ene
           velocity(nk,:)=Vk(:)
           weight(nk)=1d0/Prbl_model(ene)
           if (nk.eq.nkgoal) exit
        enddo

        weight=weight/sum(weight)

        end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

        subroutine get_eneofk(k,ene,ene_grid,vel_grid,Vk,kc)
        ! calculates the energy and velocity for the k point k(:),
        ! k (input) internal-kp
        ! ene (out) energy-kp
        ! ene_grid (input) FFT energy
        ! vel_grid (input) FFT velocity
        ! Vk (out) velocity-kp
        use modmain
        use modrand
        implicit none
        real(8) :: k(3),ene,dist,norm
        integer :: iis(3),is1,is2,is3,ik(3),j
        real(8) :: kbox(3),kboxc(3),Vk(3)
        real    :: ene_grid(inp%nr(1),inp%nr(2),inp%nr(3))
        real    :: vel_grid(inp%nr(1),inp%nr(2),inp%nr(3),3)
        real(8) :: kc(3) ! cartesian k

        ene=0.0d0
        Vk(:)=0.0d0
        norm=0.0d0
        
        kc(:)=0.0
        do j=1,3
           kc(:)=kc(:)+bvec(:,j)*k(j)
        end do
        
        ! runs over the corners of the box containing k to interpolate its energy
        do is1=0,1 ; do is2=0,1 ; do is3=0,1
           iis(1)=is1 ; iis(2)=is2 ; iis(3)=is3
           ik(:)=1+int(inp%nr(:)* k(:)) + iis(:)
           if (ik(1).gt.inp%nr(1)) ik(1)=1
           if (ik(2).gt.inp%nr(2)) ik(2)=1
           if (ik(3).gt.inp%nr(3)) ik(3)=1
           kbox(1)=float(ik(1)-1)/float(inp%nr(1))
           kbox(2)=float(ik(2)-1)/float(inp%nr(2))
           kbox(3)=float(ik(3)-1)/float(inp%nr(3))

           kboxc(:)=0.0
           do j=1,3
           kboxc(:)=kboxc(:)+bvec(:,j)*kbox(j)
           end do 

           dist=sqrt(dot_product(kboxc-kc,kboxc-kc))+1d-14
           norm=norm+1d0/dist
           ene=ene+ene_grid(ik(1),ik(2),ik(3))/dist
           Vk(:)=Vk(:)+vel_grid(ik(1),ik(2),ik(3),:)/dist
        enddo ; enddo ; enddo
        ene=ene/norm

        end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

! this routine sets the Monte Carlo acceptance distribution.
!samp%width sets the width of the Fermi level region
!sam%skin sets how abruptily this ends 
!samp%Pmin sets the asymptotic value. 
! For semiconductors the width should be of the order of the gap or 
! slightly larger, similarly the skin.

        real(8) function Prbl_model(xi)
        use modrand, only: samp  
        implicit none
        real(8) :: xi
        real(8) :: exp1,exp2

        exp1=exp((abs(xi)-samp%width)/samp%skin)
        exp2=exp(-samp%width/samp%skin)
        Prbl_model=(1d0-samp%Pmin)*(1d0+exp2)/(1d0+exp1) + samp%Pmin

        end function

