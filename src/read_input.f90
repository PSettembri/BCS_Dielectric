        subroutine read_input
        use modmain
        use modrand
        implicit none
        ! local variables        
        character(256) block,str
        integer i,ios

        !Default values
        ntasks=0
        epslat=1.d-6
        kisto=0.d0
        dbin=0.5d0
        scissor=0.d0
        combination=.false.
        
        !Default values for random variables 
                
        inp%nr(:)=120
        inp%irand=0
        inp%autokp=.true.
        inp%ef_window=0.0
        inp%fermi_velocity=.true.

        el%nob=0
        el%n1=0
        el%nktot=10000
        el%fe_shift=0.0

        samp%Pmin=0.01d0
        samp%width=0.002d0
        samp%skin=0.02d0

        qgr%nfinq(:)=30

       open(50,file='input.in',status='OLD',form='FORMATTED',iostat=ios)

        if (ios /= 0) then
          write(*,*)
          write(*,'("Error(readinput): error opening input.in")')
          write(*,*)
        stop
        end if

        10 continue
        read(50,*,end=30) block
        ! check for a comment
     if ((scan(trim(block),'!') == 1).or.(scan(trim(block),'#') == 1)) goto 10
        select case(trim(block))
        case('tasks')
          do i=1,maxtasks
            read(50,'(A)',err=20) str
            if (trim(str) == '') then
              if (i == 1) then
                write(*,*)
                write(*,'("Error(read_input): no tasks to perform")')
                write(*,*)
                stop
              end if
              ntasks=i-1
              goto 10
            end if
            read(str,*,iostat=ios) tasks(i)
            if (ios /= 0) then
              write(*,*)
              write(*,'("Error(read_input): error reading tasks")')
              write(*,'("(blank line required after tasks block)")')
              write(*,*)
              stop
            end if
          end do
          write(*,*)
          write(*,'("Error(read_input): too many tasks")')
          write(*,'("Adjust maxtasks in modmain and recompile code")')
          write(*,*)
          stop  
        case('wminrf')
                  read(50,*,err=20) wminrf
        case('wmaxrf')
                  read(50,*,err=20) wmaxrf
        case('dwrf')
                  read(50,*,err=20) dwrf
        case('smear')
                  read(50,*,err=20) smear
        case('eta')
                  read(50,*,err=20) eta
        case('epslat')
                  read(50,*,err=20) epslat
        case('kisto')
                  read(50,*,err=20) kisto
        case('dbin')
                  read(50,*,err=20) dbin
        case('scissor')
                  read(50,*,err=20) scissor
        case('nob')
                  read(50,*,err=20) el%nob
                  allocate(el%nkib(el%nob))
        case('n1')
                  read(50,*,err=20) el%n1                  
        case('nktot')
                  read(50,*,err=20) el%nktot
        case('nkib')
                if (allocated(el%nkib)) then
                if (.not. inp%autokp) then
                  read(50,*,err=20) el%nkib(:) 
                end if      
                end if 
        case('fe_shift')
                  read(50,*,err=20) el%fe_shift
        case('pmin')
                  read(50,*,err=20) samp%Pmin
        case('width')
                  read(50,*,err=20) samp%width
        case('skin')
                  read(50,*,err=20) samp%skin
        case('nr')
                  read(50,*,err=20) inp%nr(:)
        case('irand')
                  read(50,*,err=20) inp%irand
        case('autokp')
                  read(50,*,err=20) inp%autokp
        case('ef_window')
                  read(50,*,err=20) inp%ef_window
        case('fermi_velocity') 
                  read(50,*,err=20) inp%fermi_velocity
        case('nfinq')
                read(50,*,err=20) qgr%nfinq(:)
        case('combination')
                read(50,*,err=20) combination
        case('')
        goto 10
        
        end select
        goto 10

        20 continue
       write(*,*)
       write(*,'("Error(readinput): error reading from input.in")')
       write(*,'("Problem occurred in ''",A,"'' block")') trim(block)
       write(*,'("Check input convention in manual")')
       write(*,*)
        stop
        30 continue
        close(50)


        end subroutine 
