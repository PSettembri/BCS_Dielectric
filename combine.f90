        subroutine combine
        use modmain
        use modrand
        implicit none
        !local variables
        integer iw,iq,ig,jg
        real (8) xxx,yyy,skip
        complex (8) a(ngrf,ngrf,nw,nqpt)
        complex (8) b(ngrf,ngrf,nw,nqpt)
        complex (8) c(ngrf,ngrf,nw,nqpt)
        character(256) fname,fname2
       
        ! External bands contributions from Elk files

        a(:,:,:,:)=(0.d0,0.d0)

        do ig=1,ngrf
        do jg=1,ngrf

        write(fname,'("EPSELK/EPS_",3I0,"_",3I0,"_exp.OUT")') & 
          &  ivg(:,ig),ivg(:,jg)
        open(50,file=trim(fname))

        eps(:,:)=(1.d0,0.d0)
        do iq=1,nqpt 
         read(50,*)
         do iw=1,nw
           !write(*,*)iw
           read(50,*)skip,xxx,yyy
           eps(iw,iq)=complex(xxx,yyy)
         end do
         read(50,*)
       end do 
       close(50)

        a(ig,jg,:,:)=eps(:,:)

       end do
       end do 

        ! Internal bands contribution computed using VZZ
                
        do ig=1,ngrf
        do jg=1,ngrf

        write(fname2,'("EPS/eps_gg1_",3I0,"_",3I0,".dat")') &
          &  ivg(:,ig),ivg(:,jg)
        open(50,file=trim(fname2))

        eps(:,:)=(1.d0,0.d0)
        do iq=1,nqpt
         read(50,*)
         do iw=1,nw
         read(50,*)skip,xxx,yyy
         eps(iw,iq)=complex(xxx,yyy)
         end do
         read(50,*)
       end do
       close(50)

        b(ig,jg,:,:)=eps(:,:)

       end do
       end do

        c(:,:,:,:)=a(:,:,:,:)+b(:,:,:,:)
        do ig=1,ngrf
          c(ig,ig,:,:)=c(ig,ig,:,:)-1.d0
        end do

        !inversion of c
        
        do iq=1,nqpt
        do iw=1,nw
        call inv_lfe(ngrf,c(:,:,iw,iq)) 
        end do 
        end do 

        ! for the istogram of the loss function
        ! we use epsgg=1/(epsgg^-1) to include LFEs
       
        do ig=1,ngrf

        totaleps(ig,:,:)= 1.d0/(c(ig,ig,:,:))
         
        end do 

        do ig=1,ngrf
        do iq=1,nqpt
        do iw=1,nw
        
        write(1024,*)wplt(iw),dble(totaleps(ig,iw,iq)), & 
                              & aimag(totaleps(ig,iw,iq))

        end do
        write(1024,*)
        end do
        end do


        invers=.true.
        call istogram

        end subroutine

