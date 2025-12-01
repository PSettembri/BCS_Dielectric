        subroutine read_rand
        use modmain, only: T,beta,kb,ngrf
        use modrand,  only : el
        implicit none
        integer :: ik,iik
        character(len=1) :: skip
        logical :: MBgap=.false.
        real(8) ::  eps=1d-6, ee, Delt
        character(len=3) :: dummy

        open(unit=50,file='kp_random.dat',status='OLD')

        read(50,*)skip,el%nktot
        
        allocate (el%k(el%nktot,3))     ! kpoints internal coordinates
        allocate (el%energy(el%nktot))  ! Normal KS eigenvalues
        allocate (el%Delta(el%nktot))   ! Superconducting gap (KS or MB depending on 'MBgap' flag)
        allocate (el%weight(el%nktot))  ! kp weight in the random mesh 
        allocate (el%istrand(el%nktot)) ! band index (it might not match elk band indices !!!) 
        allocate (el%ikmap(el%nktot,el%nktot,2))
        allocate (el%kocc(el%nktot))
        allocate (el%igmap(el%nktot,el%nktot,ngrf))

        
        do ik=1,el%nktot
                read(50,*) el%k(ik,:),el%energy(ik),el%weight(ik),el%istrand(ik)
        end do 

        close(50)

        if (MBgap) then ! reading the physical many body gap
                open(unit=51,file='Gap_and_Z_ofk_at_wn1.dat',status='OLD')
                read(51,*) dummy,dummy,dummy,T
                el%Delta=0.0d0
                do ik=1,el%nktot
                        read(51,*,end=100) ee,Delt,iik
                        el%Delta(iik)=Delt  ! not all k points are included in this file. Those excluded are far from EF and Delta=0
                        if (abs(ee-el%energy(iik)).gt.eps) stop 'inconsistent energies in kp and gapMB files'
                end do
        !else !reading the Kohn Sham gap
        !        open(unit=51,file='gap_of_k.dat',status='OLD')
        !        read(51,*) dummy,T
        !        do ik=1,el%nktot
        !                read(51,*) ee,el%Delta(ik)
        !                if (abs(ee-el%energy(ik)).gt.eps) stop 'inconsistent energies in kp and gapKS files'
        !        end do
        endif

        100 continue
        beta=1d0/kb/T
        close(51)


        do ik=1,el%nktot
                call fermi(el%energy(ik),el%kocc(ik))
        end do 

        call grid
        
        end subroutine

