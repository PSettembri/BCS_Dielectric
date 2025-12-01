        subroutine eps_gg_file
        use modmain
        implicit none
        !local variables
        integer iw,iq,ig
        real (8) a,b,skip
        character(256) fname
       

        do ig=1,ngrf

        write(fname,'("EPS/eps_gg_",3I0,"_",3I0,".dat")') &
                        & ivg(:,ig),ivg(:,ig)
        open(50,file=trim(fname))

        eps(:,:)=(0.0,0.0)

        do iq=1,nqpt 

          read(50,*)

           do iw=1,nw
             read(50,*)skip,a,b 
             eps(iw,iq)=complex(a,b)
           end do

          read(50,*)

       ! end loop over iq  
       end do 
        
        close(50)

        totaleps(ig,:,:) = eps(:,:)

        ! end loop over ig
       end do

       ! takes the epsgg(w,q) to computed
       ! the full k istogtam of the loss function

       invers=.false.
       call istogram
        
       end subroutine

