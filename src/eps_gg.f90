        subroutine eps_gg
        use modmain
        implicit none
        !local variables
        integer iw,iq,ik,ist,jst,ig
        integer jk,jkq,isym,taskwrite
        real (8) ei,ej,eij,t1
        real (8) vkql(3),vqpl(3)
        complex (8) etasm,t2
        character(256) fname
        !real(8) occi,occj
       
        etasm=complex(0,eta)
        taskwrite=4

        do ig=1,ngrf

         if (task.eq.taskwrite) then
        write(fname,'("EPS/eps_gg_",3I0,"_",3I0,".dat")') &
                        & ivg(:,ig),ivg(:,ig)
        open(50,file=trim(fname))
        end if

        eps(:,:)=(0.d0,0.d0)

        do iq=1,nqpt 

         vqpl(:)=vql(:,iq)
         
         do ik=1,nkptnr
               
                ! k+q-vector in lattice coordinates
                vkql(:)=vkl(:,ik)+vqpl(:)
                ! equivalent reduced k-points for k and k+q
                jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
                call findkpt(vkql,isym,jkq)

                do ist=1,nstates
                       do jst=1,nstates
                                
                          ei=elkeig(ist,jk)
                          ej=elkeig(jst,jkq) 
                          
                          eij=ei-ej 

                          !call fermi(ei,occi)
                          !call fermi(ej,occj)
                          !t1=omega*(occi-occj)/(dble(nkptnr))

                t1=omega*(elkocc(ist,jk)-elkocc(jst,jkq))/(dble(nkptnr))

                          if (abs(t1) < 1.d-8) cycle

                          do iw=1,nw
                          
                t2=t1/(w(iw)+(eij)+etasm)
                eps(iw,iq)=eps(iw,iq)+t2*vzzkq(ig,ig,ist,jst,iq,ik)

                          ! end loop over iw
                          end do 
                      ! end loop over jst
                      end do 
               ! end loop over ist
               end do 

          ! end loop over ik
          end do 

         eps(:,iq)=1.d0-eps(:,iq)
         wloss(:,iq)=aimag(eps(:,iq))/(zabs(eps(:,iq))**2)
        
         if (task.eq.taskwrite) then 
          write(50,*)'#',iq,vql(:,iq)
           do iw=1,nw
     write(50,*)wplt(iw),dble(eps(iw,iq)),aimag(eps(iw,iq)),wloss(iw,iq)
           end do
          write(50,*)
         end if

       ! end loop over iq  
       end do 
        
        if (task.eq.taskwrite) then
        close(50)
        end if 

        totaleps(ig,:,:) = eps(:,:)

        ! end loop over ig
       end do

       ! takes the epsgg(w,q) to computed
       ! the full k istogtam of the loss function
       
       invers=.false.
       call istogram
        
       end subroutine

