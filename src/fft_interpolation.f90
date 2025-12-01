        subroutine fft_interpolation(ib,vr,vr_vkf)  
        ! double Fourier transform to interpolate the eigen_grid 
        !eigenvalues to a large set nr(1),nr(2),nr(3)
        !inp%nk(:) is the input grid of eigenvalues
        !inp%nr(:) is the chosen FFT grid (120 or something like that) 
        use modmain
        use modrand
        implicit none
        real :: vr(inp%nr(1),inp%nr(2),inp%nr(3))
        real :: vr_vkf(inp%nr(1),inp%nr(2),inp%nr(3),3)
        real, dimension(:,:,:), allocatable :: MATOUT2
        real, dimension(:,:,:), allocatable :: MAT1,MAT2
        real, dimension(:,:,:,:), allocatable :: MATV
        real :: fact,Rvec
        integer :: i1,i2,i3,ii1,ii2,ii3,ib,j

        allocate(MAT1(el%nk(1),el%nk(2),el%nk(3)))
        allocate(MAT2(el%nk(1),el%nk(2),el%nk(3)))
        if (inp%fermi_velocity) then 
                allocate(MATV(el%nk(1),el%nk(2),el%nk(3),3))
        end if 

        MAT2=0.0d0
        !emat is the matrix of eigenvalues from elk
        MAT1(:,:,:)=real(el%emat(:,:,:,ib),4)

        call fft3d(el%nk(1),el%nk(2),el%nk(3),MAT1,MAT2,-1)

        allocate(MATOUT2(inp%nr(1),inp%nr(2),inp%nr(3)))
        vr_vkf=0.0d0
        vr=0.0d0

        ! the matrix is reorganized for the back transform
        do i1=1,el%nk(1)
           if(i1.le.1+el%nk(1)/2) then
              ii1=i1
           else
              ii1=i1-el%nk(1)+inp%nr(1)
           endif
           do i2=1,el%nk(2)
              if(i2.le.1+el%nk(2)/2) then
                 ii2=i2
              else
                 ii2=i2-el%nk(2)+inp%nr(2)
              endif
              do i3=1,el%nk(3)
                 if(i3.le.1+el%nk(3)/2) then
                    ii3=i3
                 else
                    ii3=i3-el%nk(3)+inp%nr(3)
                 endif
         
                   vr(ii1,ii2,ii3)=MAT1(i1,i2,i3)
                   do j=1,3
            Rvec=float(ii1-1)/inp%nr(1)*real(avec(j,1),4)+   &
                 float(ii2-1)/inp%nr(2)*real(avec(j,2),4)+   &
                 float(ii3-1)/inp%nr(3)*real(avec(j,3),4)
                    if (inp%fermi_velocity) then               
                        vr_vkf(ii1,ii2,ii3,j)=MAT1(i1,i2,i3)*Rvec
                    end if 
                   enddo
            enddo
           enddo
        enddo

        fact=real(el%nk(1)*el%nk(2)*el%nk(3))
        vr=vr/fact

        MATOUT2=0.0d0
        
        call fft3d(inp%nr(1),inp%nr(2),inp%nr(3),vr,MATOUT2,1)
       
        !     CHECK QUALITY PRINT FFT FITTED EIGENVALUES
             open(1770,file='FFT_bands.dat',position='append')
             do i1=1,inp%nr(1)
                write(1770,*) i1,vr(i1,1,1),vr(i1,i1,1),vr(i1,i1,i1)
             enddo
             write(1770,*) 
             write(1770,*) 
             close(1770)

             !open(1771,file='FFT_full_eig.dat',position='append')
             !do i1=1,inp%nr(1)
             !do i2=1,inp%nr(2)
             !do i3=1,inp%nr(3)
             !   write(1771,*) i1,vr(i1,i2,i3)
             !enddo
             !enddo
             !enddo
             !write(1771,*)
             !write(1771,*)
             !close(1771)


      
        if (inp%fermi_velocity) then
           vr_vkf=vr_vkf/fact
           MATOUT2=0.0d0
           do j=1,3
           call fft3d(inp%nr(1),inp%nr(2),inp%nr(3), MATOUT2, &
                                &vr_vkf(:,:,:,j) ,1 )
           vr_vkf(:,:,:,j)=MATOUT2
           enddo
        endif   

        deallocate(MATOUT2,MAT1,MAT2)
        if(inp%fermi_velocity) deallocate(MATV)

        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine fft3d(nx,ny,nz,x,y,sgn)
      implicit none 
      integer :: nx,ny,nz,sgn
      real    :: x(nx,ny,nz),y(nx,ny,nz)
      
!	calculates the discrete fouriertransform f(i1,i2,i3)=
!	sum_(j1,j2,j3) exp(isign*i*2*pi*(j1*i1/n+j2*i2/n+j3*i3/n)) r(j1,j2,j3)
!	input:
!		n:physical dimension of the transform. it must be a product
!		  of the prime factors 2,3,5
!		n1:first dimension of x,y. n1 must always be greater or 
!		   equal than n. it is recomended to chose n1=n if n is odd
!		   and n1=n+1 if n is even to obtain optimal execution speed.
!		x(i1,i2,i3)=real(r(i1,i2,i3))
!		y(i1,i2,i3)=imag(r(i1,i2,i3))
!	output:
!		real(f(i1,i2,i3))=x(i1,i2,i3)
!		imag(f(i1,i2,i3))=y(i1,i2,i3)

      call cft(x,y,nx*ny*nz,nx ,nx                  ,sgn)
      call cft(x,y,nx*ny*nz,ny ,nx*ny               ,sgn)
      call cft(x,y,nx*ny*nz,nz ,nx*ny*nz            ,sgn)

      return
      end subroutine
