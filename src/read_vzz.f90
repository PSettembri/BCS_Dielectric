        subroutine read_vzz
        use modmain
        implicit none
        ! local variables
        integer :: iq,ik,ig,jg,ist,jst,i,skip
        character(256) fname

        do ig=1,ngrf
                do jg=1,ngrf

      write(fname,'("VZZ/VZZ_",3I0,"_",3I0,".OUT")') ivg(:,ig),ivg(:,jg)
         open(50,file=trim(fname))
        
         do i=1,11
                read(50,*)
         end do 

            do iq=1,nqpt
                do i=1,6
                  read(50,*)
                end do

              do ik=1,nkptnr

                do i=1,6
                  read(50,*)
                end do

                do ist=1,nstates

                  read(50,*)

                  do jst=1,nstates
                  
                  read(50,*) skip,vzzkq(ig,jg,ist,jst,iq,ik)
                  
                  end do

                  read(50,*)

                end do
              end do
            end do

            close(50)

          end do
        end do


        end subroutine
