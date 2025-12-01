        subroutine read_elkeig
        use modmain
        implicit none
        ! local variables        
        integer ios,i,ik,is,skip,jk,istate

           open(50,file='EIGVAL.OUT',status='OLD',form='FORMATTED',iostat=ios)

        if (ios /= 0) then
        write(*,*)
        write(*,'("Error(readinput): error opening EIGVAL.OUT")')
        write(*,*)
        stop
        end if

        elkeig(:,:)=-1000

        do i=1,2
                read(50,*)
        end do

        do ik=1,nkpt
                read(50,*)
                read(50,*)
                read(50,*)
                do istate=1,nstsv
                   if (istate.ge.nstminrf.and.istate.le.nstmaxrf) then
                         is=istate-nstminrf+1
                         read(50,*)skip,elkeig(is,ik),elkocc(is,ik)
                   else 
                         read(50,*)
                   end if 
                end do  
                read(50,*)
        end do

        close(50)

        if (scissor.ne.0) then
          do is=1,nstates
           if (any(elkeig(is,:).gt.elkfermi)) then
                elkeig(is,:)=elkeig(is,:)+scissor
           end if      
          end do
          elkfermi=elkfermi+(scissor/2.d0)
        end if


        !open(40,file='full_eig.dat',status='unknown')

        do ik=1,nkptnr

                jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))

                do is=1,nstates
                        elkeig(is,ik)=elkeig(is,jk)
                        elkocc(is,ik)=elkocc(is,jk)
                        !write(40,*)vkl(1,ik),elkeig(is,ik)
                end do 
        end do

        close(40)

        end subroutine 
