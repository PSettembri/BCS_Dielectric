        subroutine montecarlo
        use modmain
        use modrand
        implicit none
        integer :: ik,istate,is,nmax
        integer :: jk(3)


        ! states of interest
        if (el%n1.eq.0) then
                el%n1=nstminrf
        end if 
        if (el%nob.eq.0) then
                el%nob=nstates
                allocate(el%nkib(el%nob))
        end if 

        nmax=el%n1+el%nob-1

        ! initial k point grid
        el%nk(:)=ngridk(:)

        ! eigenvalues assignment
        allocate(el%emat(0:el%nk(1)-1,0:el%nk(2)-1,0:el%nk(3)-1,el%nob))

        do ik=1,nkptnr
              jk(:)=ivk(:,ik)
              do istate=1,nstates
                 if (istate.ge.el%n1.and.istate.le.nmax) then
                         is=istate-el%n1+1
                         el%emat(jk(1),jk(2),jk(3),is)=elkeig(istate,ik)
                        
                 end if 
              enddo
        enddo

        ! Eigenvalues are shifted at the fermi level and by the
        ! additional fe_shift, which is 0 if not specified
        
        call read_gap
        write(*,*) '### GAP.OUT read ###'
        
        el%emat=el%emat-elkfermi-el%fe_shift

        if (samp%width.lt.elkgap) then
                samp%width=2*elkgap
                write(*,*) '!Saxon-Wood width changed to 2*Egap'
                write(*,*) samp%width,'Ry'
        end if


        ! Setting of the final dense q grid
        qgr%nfinq(:)=30

        call set_mesh_and_el
        
        call grid

        end subroutine

