        subroutine istogram_rand
        use modmain
        use modrand
        integer ibin,ig,iq,iw,nbin
        real (8) veclen, kmax
        real (8) vgqc(3), vgql(3)
        integer, dimension(:), allocatable :: countbin
        real(8), dimension(:), allocatable :: k
        complex(8), dimension(:,:), allocatable :: epsisto
        real(8), dimension(:,:), allocatable :: lossisto

        ! takes the eps_gg(w,q) functions 
        ! with or without the LFE inclusion
        ! creates an istogram for the isotropic
        ! and angular dependent total loss function
        ! this will be used in the python code for
        ! DM detection rate computation 


        ! totaleps(ig,iw,iq) in input
        ! kmax and nbin from input.in

        do ig=1,ngrf
        do iq=1,qgr%nfinqpt

        vgql(:)=ivg(:,ig)+qgr%vfinql(:,iq)

        vgqc(:)=vgql(1)*bvec(:,1)+vgql(2)*bvec(:,2)+vgql(3)*bvec(:,3)

        veclen=sqrt(sum((vgqc(:))**2.d0)) 
        veclen=veclen*(hbarc_ev/br_ang)

        if (veclen.gt.kmax) then
                kmax=veclen
        end if

        end do 
        end do


        if (kisto.eq.0)then
                kisto=kmax
        end if 
        if (kisto.gt.kmax) then
                kisto = kmax  
                write(*,*)'Value of kisto changed to',kmax
        end if 


        nbin=int(kisto/dbin)+1

        allocate (k(nbin),countbin(nbin))
        allocate (epsisto(nw,nbin),lossisto(nw,nbin))

        k(:)=0.0d0
        epsisto(:,:)=0.0d0
        countbin(:)=0
        lossisto(:,:)=0.0d0

        do ibin=1,nbin
       
         k(ibin)=dbin*(ibin-1)
        
        end do 


        do ig=1,ngrf
        do iq=1,qgr%nfinqpt

        vgql(:)=ivg(:,ig)+qgr%vfinql(:,iq)
        vgqc(:)=vgql(1)*bvec(:,1)+vgql(2)*bvec(:,2)+vgql(3)*bvec(:,3)

        veclen=sqrt(sum((vgqc(:))**2.0))
        veclen=veclen*(hbarc_ev/br_ang)

        do ibin=1,nbin-1
        if ((veclen.ge.k(ibin)).and.(veclen.lt.k(ibin+1))) then 

                countbin(ibin)=countbin(ibin)+1
                
                epsisto(:,ibin)=epsisto(:,ibin)+qgr%totaleps(ig,:,iq)
                
        end if
        if (veclen.eq.k(nbin))then 
                countbin(nbin)=countbin(nbin)+1

                epsisto(:,nbin)=epsisto(:,nbin)+qgr%totaleps(ig,:,iq)
        end if
        end do
        
        end do
        end do 
               

        do ibin=1,nbin
        if (countbin(ibin).eq.0)then
           epsisto(:,ibin) =0
        else
        epsisto(:,ibin) = epsisto(:,ibin)/countbin(ibin)
        end if
        end do 

        ! This file can be read by DarkELF as input

        if (invers) then 
        open(50,file='eps_iso_lfe_rand.dat')
        else
        open(50,file='eps_iso_rand.dat')
        end if

        write(50,*)'# P.Settembri et al., ELK code & BCSDiel'
        !write(50,*)'#','Max k',kisto,'Number of k',nbin
        !write(50,*)'#','Max w',wplt(nw),'Number of w',nw
        do ibin=1,nbin
        do iw=1,nw
        write(50,*)wplt(iw),k(ibin)*kev_ev,dble(epsisto(iw,ibin)),  & 
                                    & aimag(epsisto(iw,ibin))
        end do
        end do
        close(50)
        

        ! Loss function for direct plot

        lossisto(:,:) = aimag(epsisto(:,:))/(zabs(epsisto(:,:))**2)

        if (invers) then
        open(50,file='loss_iso_lfe_rand.dat')
        else
        open(50,file='loss_iso_rand.dat')
        end if

        write(50,*)'#','Max k',kisto,'Number of k',nbin
        write(50,*)'#','Max w',wplt(nw),'Number of w',nw
        do ibin=1,nbin
        if (countbin(ibin).eq.0)then
           lossisto(:,ibin) =0
        end if
        do iw=1,nw
        write(50,*)wplt(iw),k(ibin),lossisto(iw,ibin)
        end do 
        end do
        close(50)


        end subroutine
