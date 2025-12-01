        subroutine eps_gg1_file_rand
        use modmain
        use modrand
        implicit none
        !local variables
        integer iw,iq,ig,jg
        real (8) xxx,yyy,skip
        complex (8) a(ngrf,ngrf,nw,qgr%nfinqpt)
        complex (8) epsrand(nw,qgr%nfinqpt)
        character(256) fname
       
        do ig=1,ngrf
        do jg=1,ngrf

        write(fname,'("EPS/eps_gg1_",3I0,"_",3I0,"_rand.dat")') & 
          &  ivg(:,ig),ivg(:,jg)
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

        a(ig,jg,:,:)=epsrand(:,:)

        ! end loop over jg
       end do
        ! end loop over ig
       end do 
                
        !inversion of a 
        
        do iq=1,qgr%nfinqpt
        do iw=1,nw
          call inv_lfe(ngrf,a(:,:,iw,iq)) 
        end do 
        end do 

        ! for the istogram of the loss function
        ! we use epsgg=1/(epsgg^-1) to include LFEs
       
        do ig=1,ngrf

        qgr%totaleps(ig,:,:)= 1.d0/(a(ig,ig,:,:))
         
        end do 

        call symmetrization

        invers=.true.
        call istogram_rand

        end subroutine

