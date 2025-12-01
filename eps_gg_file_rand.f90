        subroutine eps_gg_file_rand
        use modmain
        use modrand
        implicit none
        !local variables
        integer iw,iq,ig
        real (8) xxx,yyy,skip
        complex (8) a(ngrf,nw,qgr%nfinqpt)
        complex (8) epsrand(nw,qgr%nfinqpt)
        character(256) fname
       
        do ig=1,ngrf

        write(fname,'("EPS/eps_gg_",3I0,"_",3I0,"_rand.dat")') & 
          &  ivg(:,ig),ivg(:,ig)
        open(50,file=trim(fname))

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

        a(ig,:,:)=epsrand(:,:)

        ! end loop over ig
       end do
                
       qgr%totaleps(:,:,:)= a(:,:,:)

      ! call symmetrization

       invers=.false.
       call istogram_rand

       end subroutine

