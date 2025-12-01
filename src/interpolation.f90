        subroutine interpolation(nqsp,ngridsp,nqde,qde,epssp,epsde)
        use modmain
        use modrand
        implicit none
        ! implicit none
        integer :: iq,jq,is1,is2,is3,j
        integer :: iis(3),ik(3)
        real (8) :: kbox(3),kboxc(3),qdec(3)
        real (8) :: dist,norm
        ! input variables
        integer :: nqsp,ngridsp(3),nqde
        real (8) :: qde(3,nqde)
        complex (8) :: epssp(ngrf,ngrf,nw,nqsp)
        complex (8) :: epsde(ngrf,ngrf,nw,nqde)
 
      ! nqsp dimension of the sparse grid
      ! ngridsp(3) dimensions of the sparse grid in each direction
      ! nqde dimension of the dense grid
      ! qde(3,nqde) vector in lattice coordinates of the dense grid
      ! epssp(ngrf,ngrf,nw,nqsp) original dielectric function on the sparse grid
      ! epsde(ngrf,ngrf,nw,nqde) interpolated dielectric function on the dense grid

      ! Interpolate of the dielectric function from a sparse grid
      ! onto a denser grid


      do iq=1,nqde
       
        qdec(:)=0.d0
        do j=1,3
                qdec(:)=qdec(:)+bvec(:,j)*qde(j,iq)
        end do

        do is1=0,1 
         do is2=0,1 
          do is3=0,1
           
           iis(1) = is1 
           iis(2) = is2 
           iis(3) = is3

           ik(:) = int(ngridsp(:)*qde(:,iq)) + iis(:)
             
           if (ik(1).gt.ngridsp(1)) ik(1)=0
           if (ik(2).gt.ngridsp(2)) ik(2)=0
           if (ik(3).gt.ngridsp(3)) ik(3)=0

           kbox(1)=float(ik(1))/float(ngridsp(1)) 
           kbox(2)=float(ik(2))/float(ngridsp(2)) 
           kbox(3)=float(ik(3))/float(ngridsp(3))

           kboxc(:)=0.d0
           do j=1,3
                kboxc(:)=kboxc(:)+bvec(:,j)*kbox(j)
           end do

           dist=sqrt(dot_product(kboxc-qdec,kboxc-qdec))+1d-14
           norm=norm+1.d0/dist

           jq=ivqiqnr(ik(1),ik(2),ik(3))

           epsde(:,:,:,iq) = epsde(:,:,:,iq) + epssp(:,:,:,jq)/dist
         
          end do 
         end do 
        end do

        epsde(:,:,:,iq) = epsde(:,:,:,iq)/norm

        end do 

        end subroutine
