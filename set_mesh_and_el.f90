        subroutine set_mesh_and_el
        use modmain
        use modrand
        implicit none
        integer :: ib,i,j
        real,dimension(:,:,:,:), allocatable :: vr_vkf
        real,dimension(:,:,:), allocatable :: vr
        real(8), allocatable :: energy(:),velocity(:,:),weight(:)
        real(8), allocatable :: k_rec(:,:),k_cart(:,:)
        integer :: ind,ik,ne,ie,i1,i2,i3
        real(8) :: e0,e1,xi,samp_dist,exp1,exp2

        ! fft of the starting grid to a denser one
        allocate(vr(inp%nr(1),inp%nr(2),inp%nr(3)))
        allocate(vr_vkf(inp%nr(1),inp%nr(2),inp%nr(3),3))

        ! set up the number of random kp for each band
        if (inp%autokp) call set_nkib

        write(*,*) '! Total number of points, and per band'
        write(*,*) el%nktot, el%nkib(:) ! kp : total and for each band 
        allocate (el%k(el%nktot,3))
        allocate (el%energy(el%nktot))
        allocate (el%weight(el%nktot))
        allocate (el%istrand(el%nktot))
        allocate (el%ikmap(el%nktot,el%nktot,2))
        allocate (el%kocc(el%nktot)) 
        allocate (el%igmap(el%nktot,el%nktot,ngrf))

        call random_seed(size=i)
        allocate(inp%iseed(i))
        inp%iseed = (/ (j, j=1+inp%irand, i+inp%irand) /)  
        call random_seed (PUT=inp%iseed)

        open(801,file='kp_random.dat')
        open(803,file='kp_random_cart.dat')
        open(804,file='fv_random.dat')
        write(801,*) '#  ', el%nktot, 1 
        write(803,*) '#  ', el%nktot, 1 
        
        ind=0
        ! loop on the bands
        do ib=1,el%nob
           allocate(energy(el%nkib(ib)),velocity(el%nkib(ib),3))
           allocate(weight(el%nkib(ib)))
           allocate(k_rec(el%nkib(ib),3),k_cart(el%nkib(ib),3))

           ! read eigenvalues on a grid and perform a FFT to a denser one
           call fft_interpolation(ib,vr,vr_vkf) 

           ! Generate Random kp and interpolate energies and velocities.
           call random_kp(vr,vr_vkf,el%nkib(ib),k_rec,k_cart, &
                             & energy,velocity,weight) 

           do ik=1,el%nkib(ib)
              ind=ind+1
              el%k(ind,:)=k_rec(ik,:)
              el%energy(ind)=energy(ik)
              el%weight(ind)=weight(ik)
              el%istrand(ind)=ib

              write(801,*) k_rec(ik,:),energy(ik),weight(ik),ib      
              write(803,*) k_cart(ik,:),energy(ik),weight(ik),ib      
              write(804,*) velocity(ik,:),energy(ik),weight(ik),ib      
           enddo
           deallocate(energy,velocity,weight)
           deallocate(k_rec,k_cart)
        enddo

        close(801)
        close(803)
        close(804)
      
        do ik=1,el%nktot
                call fermi(el%energy(ik),el%kocc(ik))
        end do
        
        open(100,file='el_bands.dat')
            call write_bands_r4(real(el%emat,4),el%nk(:),el%nob,100)
        close(100)

        open(100,file='el_full_eig.dat')
             do i1=0,el%nk(1)-1
             do i2=0,el%nk(2)-1
             do i3=0,el%nk(3)-1
             do ib=1,el%nob
                write(100,*) real(i1)/el%nk(1),el%emat(i1,i2,i3,ib)
             enddo
             enddo
             enddo
             enddo
        close(100)

        open(100,file='random_bands.dat')
            call random_bands(real(el%k,4),real(el%energy,4),el%nktot,100)
        close(100)

        open(100,file='sampling_dist.dat')
        ne=1000
        do ie=1,ne
        e0=minval(el%emat)-0.05
        e1=maxval(el%emat)+0.05
        xi=e0+ie*(e1-e0)/real(ne)
        exp1=exp((abs(xi)-samp%width)/samp%skin)
        exp2=exp(-samp%width/samp%skin) 
        samp_dist=(1d0-samp%Pmin)*(1d0+exp2)/(1d0+exp1) + samp%Pmin
        write(100,*)xi,samp_dist
        end do 
        close(100)

        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine set_nkib 
        use modrand
        implicit none
        integer :: ib
        real(8) :: totw,p 
        real(8) :: band_min(el%nob),band_max(el%nob),band_avg(el%nob)

        p=0.5d0  ! this power setst the distribution function 

        do ib=1,el%nob
           band_min(ib)=minval(abs(el%emat(:,:,:,ib)))
           band_max(ib)=maxval(abs(el%emat(:,:,:,ib)))
           band_avg(ib)=(band_max(ib)+band_min(ib))/2.0
        enddo

        totw=sum( 1d0 / (  band_avg(:)+ inp%ef_window  )**p ) ! this sets the distribution function 

        el%nkib=0.0d0
        do ib=1,el%nob
           el%nkib(ib)=el%nkib(ib) + nint ( el%nktot*1d0/ & 
                      & (band_avg(ib) + inp%ef_window)**p/totw )
        enddo
        el%nkib(el%nob)=el%nkib(el%nob) - (sum(el%nkib(:)) - el%nktot)  
        ! adjust total number (for rare cases of rounding up problems)

        if (any(el%nkib.eq.0)) stop 'too few points for autokp'

        end subroutine set_nkib

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !generic band plotting routine
        subroutine write_bands_r4(mat,np,ne,io)
        implicit none
        ! io = file to write to
        ! mat= data ; np(3)=its dimension; ne=number of bands
        integer :: io,np(3),ne
        real(4) :: mat(np(1),np(2),np(3),ne)
        integer :: i
        do i=1,np(1)
           write(io,'(1000f12.6)') float(i-1)/float(np(1)),mat(i,1,1,:)
        enddo
        write(io,*)
        do i=1,np(1)
           write(io,'(1000f12.6)') float(i-1)/float(np(1)),mat(i,i,1,:)
        enddo
        write(io,*)
        do i=1,np(1)
           write(io,'(1000f12.6)') float(i-1)/float(np(1)),mat(i,i,i,:)
        enddo
        write(io,*)
        end subroutine write_bands_r4

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !band plotting routine for random kp
        subroutine random_bands(k,energy,nk,io)
        implicit none
        ! io = file to write to
        ! energy= random points eigenvalues 
        ! k = direct lattice coordinates of k points
        ! nk = number of k points
        integer :: io,nk
        real(4) :: k(nk,3),energy(nk)
        integer :: ik
        real(4) :: width

        width=0.03

        do ik=1,nk
            if (abs(k(ik,2)).le.width.or.abs(1-k(ik,2)).le.width) then 
            if (abs(k(ik,3)).le.width.or.abs(1-k(ik,3)).le.width) then
                    write(io,*) k(ik,1), energy(ik)
            end if 
            end if 
        end do 
        write(io,*)
        do ik=1,nk
            if (abs(k(ik,1)-k(ik,2)).le.width) then
            if (abs(k(ik,3)).le.width.or.abs(1-k(ik,3)).le.width) then
                    write(io,*) k(ik,1), energy(ik)
            end if 
            end if 
        end do
        write(io,*)
        do ik=1,nk
            if (abs(k(ik,1)-k(ik,2)).le.width) then
            if (abs(k(ik,1)-k(ik,3)).le.width) then    
                    write(io,*) k(ik,1), energy(ik)
            end if 
            end if 
        end do

        end subroutine 

