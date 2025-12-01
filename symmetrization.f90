        subroutine symmetrization
        use modmain
        use modrand
        implicit none
        complex (8) symeps(ngrf,nw,qgr%nfinqpt) 
        integer :: iq,ig,isym,i1,i2,i3,irq,irg
        real (8) :: q(3),Rq(3)

        symeps(:,:,:)=0.0d0

        write(*,*) qgr%nfinqpt
        write(*,*) qgr%nfinq(1)        
        

        do iq=1,qgr%nfinqpt

           q(:)=qgr%vfinql(:,iq)
                
           do isym=1,nsymcrys

           Rq(:) = symlatspa(1,:,isym)*q(1)    &
                  +symlatspa(2,:,isym)*q(2)    &
                  +symlatspa(3,:,isym)*q(3)

           i1 = int(Rq(1)*dble(qgr%nfinq(1)))
           i2 = int(Rq(2)*dble(qgr%nfinq(2)))
           i3 = int(Rq(3)*dble(qgr%nfinq(3)))

           write(*,*) i1,i2,i3
           irq = qgr%ivfinqiq(i1,i2,i3)

                do ig=1,ngrf

                irg = rotg(isym,ig) 

                symeps(ig,:,iq)=symeps(ig,:,iq) + & 
                                    & qgr%totaleps(irg,:,irq)

                end do
            end do 
        end do

        symeps(:,:,:)=symeps(:,:,:)/dble(nsymcrys)

        qgr%totaleps(:,:,:)=symeps(:,:,:)

       end subroutine
