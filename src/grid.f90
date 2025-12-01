        subroutine grid
        use modmain
        use modrand
        implicit none
        ! local variables
        real (8) :: vpl(3),q(3),Rk(3)
        integer :: irand,jrand,isym,irrq
        integer :: i1,i2,i3,iq,ik,ind
        integer :: indIBZq(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1,2)
        ! el%nktot number of random k points
        ! el%k(el%nktot,3) random k points vectors
        ! vql(j,i) j=1,3 i=1,nqptnr q-point vectors
        ! el%ikmap(el%nktot,el%nktot,2) mapping in output

        do iq=1,nqptnr
                vpl(:)=vql(:,iq)
                call findkpt(vpl,isym,ik) 
                indIBZq(ivq(1,iq),ivq(2,iq),ivq(3,iq),1) = ik
                indIBZq(ivq(1,iq),ivq(2,iq),ivq(3,iq),2) = isym
        end do 

        ! run on the random points
        do irand=1,el%nktot 
        if (mod(irand,2000) == 0) then
                write(*,*)'Step',irand,'out of',el%nktot
        end if
                do jrand=1,el%nktot 

                        ! q = k'-k
                        q(:)=el%k(jrand,:)-el%k(irand,:) 

                        ! shift of the grid is given only on kgrid
                        ! but kgrid=qgrid so it is equivalent

                        q(:)=q(:)-vkloff(:)/dble(ngridk(:))
                        i1=modulo(nint(q(1)*ngridq(1)),ngridq(1))
                        i2=modulo(nint(q(2)*ngridq(2)),ngridq(2))
                        i3=modulo(nint(q(3)*ngridq(3)),ngridq(3))

                        ! get the irreducible q index
                        irrq = indIBZq(i1,i2,i3,1)   

                        ! rotate k accordingly 
                        isym = indIBZq(i1,i2,i3,2)
                        
                        Rk(:)=symlatspa(1,:,isym)*el%k(irand,1)    &
                                +symlatspa(2,:,isym)*el%k(irand,2) &
                                +symlatspa(3,:,isym)*el%k(irand,3)

                        Rk(:)=Rk(:)-vkloff(:)/dble(ngridk(:))
                        i1=modulo(nint(Rk(1)*ngridk(1)),ngridk(1))
                        i2=modulo(nint(Rk(2)*ngridk(2)),ngridk(2))
                        i3=modulo(nint(Rk(3)*ngridk(3)),ngridk(3))

                        ik=ivkiknr(i1,i2,i3)

                        el%igmap(irand,jrand,:)=rotg(isym,:)

                        !the irand and jrand pair corresponds to the 
                        !the following indices in the matrix element list
                        el%ikmap(irand,jrand,1) = irrq
                        el%ikmap(irand,jrand,2) = ik 
                end do
        end do   

        
        ! We defined the dense q-grid

        qgr%nfinqpt = qgr%nfinq(1)*qgr%nfinq(2)*qgr%nfinq(3)

        allocate (qgr%vfinql(3,qgr%nfinqpt))
        allocate (qgr%ivfinq(3,qgr%nfinqpt))
        allocate (qgr%ivfinqiq(0:qgr%nfinq(1)-1, & 
                        & 0:qgr%nfinq(2)-1,0:qgr%nfinq(3)-1))
        allocate (qgr%eps(nw,qgr%nfinqpt))
        allocate (qgr%totaleps(ngrf,nw,qgr%nfinqpt))

        ind=0

        do i1=0,qgr%nfinq(1)-1
        do i2=0,qgr%nfinq(2)-1
        do i3=0,qgr%nfinq(3)-1
        
        ind=ind+1

        qgr%vfinql(1,ind)=i1/dble(qgr%nfinq(1))
        qgr%vfinql(2,ind)=i2/dble(qgr%nfinq(2))
        qgr%vfinql(3,ind)=i3/dble(qgr%nfinq(3))

        qgr%ivfinq(1,ind)=i1
        qgr%ivfinq(2,ind)=i2
        qgr%ivfinq(3,ind)=i3

        qgr%ivfinqiq(i1,i2,i3)=ind

        end do 
        end do 
        end do 

        end subroutine

