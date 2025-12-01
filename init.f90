        subroutine init
        use modmain
        implicit none
        integer i,iq,ig,jg,isym,jq
        real (8) diff,Rg(3),q(3),Rq(3)
        
        ! definition of the frequencies
        nw=int((wmaxrf-wminrf)/dwrf)+1
        allocate (w(nw),wplt(nw))
        do i=1,nw
                w(i)=wminrf+dwrf*(i-1)
        end do 
        wplt(:)=w(:)*ha_ev

        ! number of states considered
        if (nstminrf.eq.0) then
           nstminrf=1
        end if
        if (nstmaxrf.eq.0) then
           nstmaxrf=nstsv
        end if    
        nstates=nstmaxrf-nstminrf+1
       
        ! allocating the VZZ matrix elements
        allocate (vzzkq(ngrf,ngrf,nstates,nstates,nqpt,nkptnr))
        ! allocating the dielectric and loss functions
        allocate (wloss(nw,nqpt),eps(nw,nqpt),totaleps(ngrf,nw,nqpt))
        ! allocating the ELK eigenvalues
        allocate (elkeig(nstates,nkptnr))
        ! allocating the ELK occupations
        allocate (elkocc(nstates,nkptnr))

        invers=.false.


        ! operations to go from qirr to q

        allocate(mapqirrtoq(nqptnr))

        do iq=1,nqptnr

          jq=ivqiq(ivq(1,iq),ivq(2,iq),ivq(3,iq))
          q(:)=vql(:,jq)

          do isym=1,nsymcrys

            Rq(:) = symlatspa(1,:,isym)*q(1)    &
                   +symlatspa(2,:,isym)*q(2)    &
                   +symlatspa(3,:,isym)*q(3)

            diff=sum(abs(Rq(:)-vql(:,iq)))

            if (diff.le.epslat) then

               mapqirrtoq(iq)=isym

            end if

         end do
       end do


       ! map of rotated Gs
        allocate (rotg(nsymcrys,ngrf))

        do ig=1,ngrf

        do isym=1,nsymcrys

            Rg(:)=symlatspa(1,:,isym)*ivg(1,ig)    &
                  +symlatspa(2,:,isym)*ivg(2,ig)   &
                  +symlatspa(3,:,isym)*ivg(3,ig)

            do jg=1,ngrf
                 if (sum(abs(ivg(:,jg)-Rg(:))).eq.0) then
                          rotg(isym,ig)=jg
                 end if
            end do
        end do
        end do


        end subroutine
