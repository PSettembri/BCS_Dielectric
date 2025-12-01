        subroutine eps_gg1_rand
        use modmain
        use modrand
        implicit none
        !local variables
        integer iw,iq,ik,ist,jst,ig,jg
        integer taskwrite,gloop,g1loop
        real (8) ei,ej,eij,t1
        !real (8) vkql(3),vqpl(3)
        complex (8) etasm,t2
        real (8) totalloss
        complex (8) a(ngrf,ngrf,nw,qgr%nfinqpt)
        character(256) fname
        real (8) occi,occj,qrand(3),wk,wk1
        integer iqdens,iq1,iq2,iq3,irand,jrand
       
        etasm=complex(0,eta)
        taskwrite=16

        qgr%eps(:,:)=0.0d0
        a(:,:,:,:)=0.0d0

        do irand=1,el%nktot
              
            ei  = el%energy(irand)
            ist = el%istrand(irand)
            occi= el%kocc(irand)

            do jrand=1,el%nktot
                            
                    jst=el%istrand(jrand)
                    iq=el%ikmap(irand,jrand,1)
                    ik=el%ikmap(irand,jrand,2)

                    ! qrand = k'rand-krand
                    qrand(:) = el%k(jrand,:)-el%k(irand,:)

                    iq1=modulo(nint(qrand(1)*qgr%nfinq(1)),qgr%nfinq(1))
                    iq2=modulo(nint(qrand(2)*qgr%nfinq(2)),qgr%nfinq(2))
                    iq3=modulo(nint(qrand(3)*qgr%nfinq(3)),qgr%nfinq(3))

                    iqdens=qgr%ivfinqiq(iq1,iq2,iq3)

                    ! energy difference
                    ej = el%energy(jrand)
                    eij = ei-ej
                    occj = el%kocc(jrand)
                    wk = el%weight(irand)
                    wk1 = el%weight(jrand)

                    t1=omega*(occi-occj)*(wk*wk1)

                    if (abs(t1) < 1.d-8) cycle

                    do gloop=1,ngrf
                        
                      ig = el%igmap(irand,jrand,gloop)

                      do g1loop=1,ngrf

                        jg = el%igmap(irand,jrand,g1loop)

                        do iw=1,nw
                          
                        t2=t1/(w(iw)+(eij)+etasm)
        a(gloop,g1loop,iw,iqdens) = a(gloop,g1loop,iw,iqdens)+ & 
                               & t2*vzzkq(ig,jg,ist,jst,iq,ik)

                        ! end loop over iw
                          end do 
                ! end loop over g1loop
                  end do 
                ! end loop over gloop           
                end do 
        ! end loop over jrand
          end do 
        ! end loop over irand
          end do 

        a(:,:,:,:)=-a(:,:,:,:)
        do ig=1,ngrf
          a(ig,ig,:,:)=1.d0+a(ig,ig,:,:)
        end do 

        if (task.eq.taskwrite) then

        do ig=1,ngrf
        do jg=1,ngrf

        qgr%eps(:,:) = a(ig,jg,:,:)

        write(fname,'("EPS/eps_gg1_",3I0,"_",3I0,"_rand.dat")') &
          &  ivg(:,ig),ivg(:,jg)
        open(50,file=trim(fname))
          do iqdens=1,qgr%nfinqpt
                write(50,*)'#',iqdens,qgr%vfinql(:,iqdens) 
                do iw=1,nw
                   write(50,*)wplt(iw),dble(qgr%eps(iw,iqdens)), &
                               & aimag(qgr%eps(iw,iqdens))   
                end do
           write(50,*)
           end do 
         close(50)
         end do 
         end do 

         end if 

        if (.not.combination) then
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

          open(unit=50,file='epsgg_frominverse.dat',status='unknown')
          do ig=1,ngrf
          do iq=1,qgr%nfinqpt
          do iw=1,nw  

          totalloss=aimag(qgr%totaleps(ig,iw,iq))/&
                &(zabs(qgr%totaleps(ig,iw,iq))**2)

          write(50,*)wplt(iw),dble(qgr%totaleps(ig,iw,iq)), & 
                 & aimag(qgr%totaleps(ig,iw,iq)),totalloss
          end do 
                  write(50,*)
          end do
                  write(50,*)
          end do 
          close(50)

          call symmetrization

          invers=.true.
          call istogram_rand

        end if 

        end subroutine

