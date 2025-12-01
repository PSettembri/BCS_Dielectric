        subroutine combine_rand
        use modmain
        use modrand
        implicit none
        !local variables
        integer :: iw,iq,jq,ig,jg
        integer :: igindex,jgindex,isym
        !integer :: iq1,iq2,iq3,iqgrid
        real (8) xxx,yyy,skip
        complex (8) epsgrid(nw,nqpt),epsrand(nw,qgr%nfinqpt)
        complex (8) a(ngrf,ngrf,nw,nqpt)
        complex (8) abz(ngrf,ngrf,nw,nqptnr)
        complex (8) ainter(ngrf,ngrf,nw,qgr%nfinqpt)
        complex (8) b(ngrf,ngrf,nw,qgr%nfinqpt)
        complex (8) c(ngrf,ngrf,nw,qgr%nfinqpt)
        character(256) fname,fname2
       
        ! External bands contributions from Elk files

        a(:,:,:,:)=(0.d0,0.d0)

        do ig=1,ngrf
        do jg=1,ngrf

        write(fname,'("EPSELK/EPS_",3I0,"_",3I0,"_exp.OUT")') & 
          &  ivg(:,ig),ivg(:,jg)
        open(50,file=trim(fname))

        epsgrid(:,:)=(1.d0,0.d0)
        do iq=1,nqpt 
         read(50,*)
         do iw=1,nw
           !write(*,*)iw
           read(50,*)skip,xxx,yyy
           epsgrid(iw,iq)=complex(xxx,yyy)
         end do
         read(50,*)
       end do 
       close(50)

        a(ig,jg,:,:)=epsgrid(:,:)

       end do
       end do 

       ! Passing from the elk dielectric function on the IBZ to the full BZ

       do iq=1,nqptnr
        do ig=1,ngrf
          do jg=1,ngrf

              jq=ivqiq(ivq(1,iq),ivq(2,iq),ivq(3,iq))
              isym=mapqirrtoq(iq)

              igindex=rotg(isym,ig)
              jgindex=rotg(isym,ig)

              abz(igindex,jgindex,:,iq)=a(ig,jg,:,jq)

           end do 
         end do 
       end do  

        ! Internal bands contribution computed using VZZ and random
                
        do ig=1,ngrf
        do jg=1,ngrf

        write(fname2,'("EPS/eps_gg1_",3I0,"_",3I0,"_rand.dat")') &
          &  ivg(:,ig),ivg(:,jg)
        open(50,file=trim(fname2))

        epsrand(:,:)=(1.d0,0.d0)
        do iq=1,qgr%nfinqpt
         read(50,*)
         do iw=1,nw
         read(50,*)skip,xxx,yyy
         epsrand(iw,iq)=complex(xxx,yyy)
         end do
         read(50,*)
       end do
       close(50)

        b(ig,jg,:,:)=epsrand(:,:)

       end do
       end do


       ! Interpolation 
       ! abz(ngrf,ngrf,nw,nqptnr)
       ! ainter(ngrf,ngrf,nw,qgr%nfinqpt)
       ! b(ngrf,ngrf,nw,qgr%nfinqpt)
       ! c(ngrf,ngrf,nw,qgr%nfinqpt)

       ! call interpolation(nqsp,ngridsp,nqde,ngridde,qde,epssp,epsde)

         ainter(:,:,:,:)=0.d0

         call interpolation(nqptnr,ngridq,qgr%nfinqpt,qgr%nfinq, &
                                & qgr%vfinql,abz,ainter) 

         c(:,:,:,:)=ainter(:,:,:,:)+b(:,:,:,:)


      ! Nearest neighbor 
      ! abz(ngrf,ngrf,nw,nqptnr)
      ! b(ngrf,ngrf,nw,qgr%nfinqpt)
      ! c(ngrf,ngrf,nw,qgr%nfinqpt) 

      ! do iq=1,qgr%nfinqpt
      !
      !  iq1=modulo(nint(qgr%vfinql(1,iq)*ngridq(1)),ngridq(1))
      !  iq2=modulo(nint(qgr%vfinql(2,iq)*ngridq(2)),ngridq(2))
      !  iq3=modulo(nint(qgr%vfinql(3,iq)*ngridq(3)),ngridq(3))
      !
      !  iqgrid = ivqiqnr(iq1,iq2,iq3)
      !
      !  closest elk-grid q vector to the dense grid
      !  c(:,:,iq,:)=abz(:,:,iqgrid,:)+b(:,:,iq,:)
      !
      ! end do 


        do ig=1,ngrf
          c(ig,ig,:,:)=c(ig,ig,:,:)-1.d0
        end do

        !inversion of c
        
        do iq=1,qgr%nfinqpt
        do iw=1,nw
        call inv_lfe(ngrf,c(:,:,iw,iq)) 
        end do 
        end do 

        ! for the istogram of the loss function
        ! we use epsgg=1/(epsgg^-1) to include LFEs
       
        do ig=1,ngrf

        qgr%totaleps(ig,:,:)= 1.d0/(c(ig,ig,:,:))
         
        end do 

        do ig=1,ngrf
        do iq=1,qgr%nfinqpt
        do iw=1,nw
        
        write(1024,*)wplt(iw),dble(qgr%totaleps(ig,iw,iq)), & 
                              & aimag(qgr%totaleps(ig,iw,iq))

        end do
        write(1024,*)
        end do
        end do

        call symmetrization

        invers=.true.
        call istogram_rand

        end subroutine
