        subroutine read_var
        use modmain
        implicit none
        ! local variables        
        character(256) block
        integer ios,skip,i,j,k
        !!!!!!!
        !! Default values
        !!!!!!!
        nkptnr = 0
        nqptnr = 0
        nqpt = 0
        ngridk(:)= 0
        ngridq(:)= 0
        nstmaxrf=0
        nstminrf=0

        open(50,file='VARIABLES.OUT',status='OLD',form='FORMATTED',iostat=ios)

        if (ios /= 0) then
          write(*,*)
          write(*,'("Error(readinput): error opening VARIABLES.OUT")')
          write(*,*)
        stop
        end if

        10 continue
        read(50,*,end=30) block
        ! check for a comment
     if ((scan(trim(block),'!') == 1).or.(scan(trim(block),'#') == 1)) goto 10
        select case(trim(block))
        case('gmaxrf')
                  read(50,*,err=20)
                  read(50,*,err=20) gmaxrf
        case('ngrf')
                  read(50,*,err=20)
                  read(50,*,err=20) ngrf
                  if (.not. allocated(ivg)) allocate (ivg(3,ngrf))
        case('ivg')
                  if (ngrf.ne.0) then
                  read(50,*)
                  do i=1,ngrf
                  do j=1,3
                  read(50,*) ivg(j,i)
                  end do
                  end do
                  end if
        case('nstminrf')
                  read(50,*,err=20)
                  read(50,*,err=20) nstminrf
        case('nstmaxrf')
                  read(50,*,err=20)
                  read(50,*,err=20) nstmaxrf
        case('nstsv')
                  read(50,*,err=20)
                  read(50,*,err=20) nstsv
        case('ivk')
                  if (nkptnr.ne.0) then
                  if (.not. allocated(ivk)) then
                  allocate (ivk(3,nkptnr))
                  read(50,*)
                  do i=1,nkptnr
                  do j=1,3
                        read(50,*)ivk(j,i)
                  end do 
                  end do 
                  end if
                  end if 
         case('ivq')
                  if (nqptnr.ne.0) then
                  if (.not. allocated(ivq)) then
                  allocate (ivq(3,nqptnr))
                  read(50,*)
                  do i=1,nqptnr
                  do j=1,3
                        read(50,*)ivq(j,i)
                  end do
                  end do
                  end if
                  end if
        case('ivkik')
                  read(50,*,err=20) skip,nkptnr
                  if (.not. allocated(vkl)) allocate (vkl(3,nkptnr))
                  if (.not. allocated(ivkik)) then
                  if (ngridk(1).ne.0) then
              allocate(ivkik(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
                          do k=0,ngridk(3)-1
                          do j=0,ngridk(2)-1
                          do i=0,ngridk(1)-1
                            read(50,*)ivkik(i,j,k)
                          end do 
                          end do 
                          end do 
                  end if
                  end if  
        case('ivqiq')
                  read(50,*,err=20) skip,nqptnr
                  if (.not. allocated(vql)) allocate (vql(3,nqptnr))
                  if (.not. allocated(ivqiq)) then
                  if (ngridq(1).ne.0) then
              allocate(ivqiq(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
                          do k=0,ngridq(3)-1
                          do j=0,ngridq(2)-1
                          do i=0,ngridq(1)-1
                            read(50,*)ivqiq(i,j,k)
                          end do
                          end do
                          end do
                  end if
                  end if                
        case('ivkiknr')
                  read(50,*,err=20) skip,nkptnr
                  if (.not. allocated(vkl)) allocate (vkl(3,nkptnr))
                  if (.not. allocated(ivkiknr)) then
                  if (ngridk(1).ne.0) then
           allocate(ivkiknr(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
                          do k=0,ngridk(3)-1
                          do j=0,ngridk(2)-1
                          do i=0,ngridk(1)-1
                            read(50,*)ivkiknr(i,j,k)
                          end do
                          end do
                          end do
                  end if
                  end if
        case('ivqiqnr')
                  read(50,*,err=20) skip,nqptnr
                  if (.not. allocated(vql)) allocate (vql(3,nqptnr))
                  if (.not. allocated(ivqiqnr)) then
                  if (ngridq(1).ne.0) then
           allocate(ivqiqnr(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
                          do k=0,ngridq(3)-1
                          do j=0,ngridq(2)-1
                          do i=0,ngridq(1)-1
                            read(50,*)ivqiqnr(i,j,k)
                          end do
                          end do
                          end do
                  end if
                  end if                  
        case('vkl')
                  if (nkptnr.ne.0) then
                  read(50,*)
                  do i=1,nkptnr
                  do j=1,3
                  read(50,*) vkl(j,i)
                  end do 
                  end do
                  end if
       case('vql')
                  if (nqptnr.ne.0) then
                  read(50,*)
                  do i=1,nqptnr
                  do j=1,3
                  read(50,*) vql(j,i)
                  end do
                  end do
                  end if
        case('nkpt')
                  read(50,*,err=20)
                  read(50,*,err=20)nkpt        
        case('nqpt')
                  read(50,*,err=20)
                  read(50,*,err=20)nqpt
        case('vkloff')
                read(50,*,err=20)
                do i=1,3
                read(50,*,err=20)vkloff(i)
                end do 
        case('ngridk')
                read(50,*,err=20)
                do i=1,3
                read(50,*,err=20)ngridk(i)
                end do 
        case('ngridq')
                read(50,*,err=20)
                do i=1,3
                read(50,*,err=20)ngridq(i)
                end do
        case('scissor')
                if (scissor.eq.0)then
                read(50,*,err=20)
                read(50,*,err=20)scissor
                end if
        case('')
        goto 10
       
      !  case default
      !   write(*,*)
      !   write(*,'("Error(readinput): invalid block name : ",A)') trim(block)
      !   write(*,*)
      !  stop
        
        end select
        goto 10

        20 continue
       write(*,*)
       write(*,'("Error(readinput): error reading from VARIABLES.OUT")')
       write(*,'("Problem occurred in ''",A,"'' block")') trim(block)
       write(*,'("Check input convention in manual")')
       write(*,*)
        stop
        30 continue
        close(50)


        end subroutine 
