        Program bcsdiel
                use modmain
               

                write(*,*) '### START ###'
                call cpu_time(tstart)


                ! read input.in file
                call read_input
                write(*,*) '### Input.in read ###'

                do itask=1,ntasks
                  task=tasks(itask)

            write(*,*)
            write(*,'(" ### Current task : ",I6," ###")') task
            write(*,*)

                  select case(task)
                  case(0)
                        ! read VARIABLES.OUT file
                        call read_var
                        write(*,*) '### VARIABLES.OUT read ###'
                        ! read LATTICE.OUT
                        call read_lat
                        write(*,*) '### LATTICE.OUT read ###'
                        ! read GVECRF.OUT
                        call read_gvec
                        write(*,*) '### GVECRF.OUT read ###'
                        ! read EFERMI.OUT
                        call read_fermi
                        write(*,*) '### EFERMI.OUT read ###'
                        ! read SYMCRYS.OUT
                        call read_sym
                        write(*,*) '### SYMCRYS.OUT read ###'
                        ! initialise variables
                        call init
                        write(*,*) '### Initialized variables ###' 
                        ! read EIGVAL.OUT
                        call read_elkeig
                        write(*,*) '### EIGVAL.OUT read ###'
                  case(1)
                        ! read VZZ.OUT        
                        call read_vzz 
                        write(*,*) '### VZZ.OUT read ###'
                  case(2)
                        ! compute eps for q=0
                         call eps00_Gamma
                        write(*,*) '### eps00_Gamma.dat written ###'
                  case(3)
                        ! compute eps00
                         call eps00_q
                         write(*,*) '### eps00_q.dat written ###'
                  case(4)
                        ! compute eps_gg & writes to file
                        call eps_gg
                        write(*,*) '### eps_gg.dat written ###'
                        write(*,*) '### eps_iso.dat written ###'
                        write(*,*) '### loss_iso.dat written ###'
                  case(5)
                        ! compute eps_gg
                        call eps_gg
                        write(*,*) '### eps_iso.dat written ###'
                        write(*,*) '### loss_iso.dat written ###'
                  case(6)
                        ! compute eps_gg istogram from file
                        call eps_gg_file
                        write(*,*) '### eps_iso.dat written ###'
                        write(*,*) '### loss_iso.dat written ###'
                  case(7)
                        ! compute eps_gg' & istogram & writes to file
                        call eps_gg1
                        write(*,*) '### eps_gg1.dat written ###'
                        write(*,*) '### eps_iso_lfe.dat written ###'
                        write(*,*) '### loss_iso_lfe.dat written ###'
                  case(8)
                        ! compute eps_gg'
                        call eps_gg1
                        write(*,*) '### eps_iso_lfe.dat written ###'
                        write(*,*) '### loss_iso_lfe.dat written ###'
                  case(9)
                        ! compute eps_gg' istogram from file
                        call eps_gg1_file
                        write(*,*) '### eps_iso_lfe.dat written ###'
                        write(*,*) '### loss_iso_lfe.dat written ###'
                  case(10)
                        ! combines the eps from elk to the grid vzz one
                        call combine
                        write(*,*) '### total dielectric elk+grid & 
                                    & function computed ###'
                  case(11)
                        ! generates random k points
                        call montecarlo
                        write(*,*) '### kp_random.dat written ###'
                  case(12)  
                        ! reads the random k points file
                        call read_rand
                        write(*,*) '### kp_random.dat read ###'    
                  case(13)
                        ! compute eps_gg with random & istogram & writes to file
                        call eps_gg_rand
                        write(*,*) '### eps_gg_rand.dat written ###'
                        write(*,*) '### eps_iso_rand.dat written ###'
                        write(*,*) '### loss_iso_rand.dat written ###'
                  case(14)
                        ! compute eps_gg with random & istogram
                        call eps_gg_rand
                        write(*,*) '### eps_iso_rand.dat written ###'
                        write(*,*) '### loss_iso_rand.dat written ###'
                  case(15)
                        ! compute epso_gg istogram starting from file
                        call eps_gg_file_rand
                        write(*,*) '### eps_iso_rand.dat written ###'
                        write(*,*) '### loss_iso_rand.dat written ###'
                  case(16)
                        ! compute eps_gg' with random & istogram & writes to file
                        call eps_gg1_rand
                      write(*,*) '### eps_gg1_rand.dat written ###'
                      write(*,*) '### eps_iso_lfe_rand.dat written ###'
                      write(*,*) '### loss_iso_lfe_rand.dat written ###'
                  case(17)
                        ! compute eps_gg' with random & istogram
                        call eps_gg1_rand
                      write(*,*) '### eps_iso_lfe_rand.dat written ###'
                      write(*,*) '### loss_iso_lfe_rand.dat written ###'
                  case(18)
                        ! compute eps_gg' istogram starting from file
                        call grid
                        call eps_gg1_file_rand
                        write(*,*) '### eps_iso_lfe.dat written ###'
                        write(*,*) '### loss_iso_lfe.dat written ###'
                  case(19)
                        ! combines the eps from elk to the random vzz one
                        call combine_rand
                        write(*,*) '### total dielectric elk+random &
                                    & function computed ###'
                  case default
                    write(*,*)
                    write(*,'("Error: task not defined : ",I8)') task
                    write(*,*)
                    stop
                end select
                end do

                write(*,*)
                call cpu_time(tend)
                write(*,*) 'CPU TIME:',tend-tstart,'s'
                write(*,*) '### END ###'

        end program
