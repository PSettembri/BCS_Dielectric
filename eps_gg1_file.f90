        subroutine eps_gg1_file
        use modmain
        implicit none
        !local variables
        integer iw,iq,ig,jg
        real (8) b,c,skip
        complex (8) a(ngrf,ngrf,nw,nqpt)
        character(256) fname
       
        do ig=1,ngrf
        do jg=1,ngrf

        write(fname,'("EPS/eps_gg1_",3I0,"_",3I0,".dat")') & 
          &  ivg(:,ig),ivg(:,jg)
        open(50,file=trim(fname))

        eps(:,:)=(0.0,0.0)

        do iq=1,nqpt 

         read(50,*)

         do iw=1,nw
         
         read(50,*)skip,b,c
         eps(iw,iq)=complex(b,c)

         end do

         read(50,*)

       ! end loop over iq  
       end do 

       close(50)

        a(ig,jg,:,:)=eps(:,:)

        ! end loop over ig
       end do
        ! end loop over jg
       end do 
                
        !inversion of a 
        
        do iq=1,nqpt
        do iw=1,nw
        call inv_lfe(ngrf,a(:,:,iw,iq)) 
        end do 
        end do 

        ! for the istogram of the loss function
        ! we use epsgg=1/(epsgg^-1) to include LFEs
       
        do ig=1,ngrf

        totaleps(ig,:,:)= 1.d0/(a(ig,ig,:,:))
         
        end do 

        invers=.true.
        call istogram

        end subroutine

