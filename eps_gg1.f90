        subroutine eps_gg1
        use modmain
        implicit none
        !local variables
        integer iw,iq,ik,ist,jst,ig,jg
        integer jk,jkq,isym,taskwrite
        real (8) ei,ej,eij,t1
        real (8) vkql(3),vqpl(3)
        complex (8) etasm,t2
        real (8) totalloss
        complex (8) a(ngrf,ngrf,nw,nqpt)
        character(256) fname
        !real(8) occi,occj
       
        etasm=complex(0,eta)
        taskwrite=7


        do ig=1,ngrf
        do jg=1,ngrf

        if (task.eq.taskwrite) then
        write(fname,'("EPS/eps_gg1_",3I0,"_",3I0,".dat")') & 
          &  ivg(:,ig),ivg(:,jg)
        open(50,file=trim(fname))
        end if

        eps(:,:)=(0.0d0,0.0d0)

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
                eps(iw,iq)=eps(iw,iq)+t2*vzzkq(ig,jg,ist,jst,iq,ik)

                          ! end loop over iw
                          end do 
                      ! end loop over jst
                      end do 
               ! end loop over ist
               end do 

          ! end loop over ik
          end do 

         eps(:,iq)=-eps(:,iq)
         if (ig.eq.jg) then
         eps(:,iq)=1.d0+eps(:,iq)
         end if

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

        a(ig,jg,:,:)=eps(:,:)
        

        ! end loop over ig
       end do
        ! end loop over jg
       end do 
                

       if (.not.combination) then

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

          open(unit=50,file='epsgg_frominverse.dat',status='unknown')
          do ig=1,ngrf
          do iq=1,nqpt
          do iw=1,nw  
       totalloss=aimag(totaleps(ig,iw,iq))/(zabs(totaleps(ig,iw,iq))**2)
          write(50,*)wplt(iw),dble(totaleps(ig,iw,iq)), & 
                 & aimag(totaleps(ig,iw,iq)),totalloss
          end do 
                  write(50,*)
          end do
                  write(50,*)
          end do 
          close(50)

          invers=.true.
          call istogram
        
        end if

        end subroutine

