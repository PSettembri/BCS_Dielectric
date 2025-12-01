        subroutine eps00_Gamma
        use modmain
        implicit none
        !local variables
        integer iw,ik,ist,jst
        real (8) ei,ej,eij,t1
        complex (8) etasm,t2
        !real(8) occi,occj
       
        open(unit=50,file='eps00_Gamma.dat',status='unknown')

        etasm=complex(0.d0,eta)
        eps(:,:)=(0.d0,0.d0)

         do ik=1,nkptnr
                do ist=1,nstates
                       do jst=1,nstates
                                
                          ei=elkeig(ist,ik)
                          ej=elkeig(jst,ik)
                          
                          eij=ei-ej 

                          !call fermi(ei,occi)
                          !call fermi(ej,occj)
                          !t1=omega*(occi-occj)/(dble(nkptnr))

                 t1=omega*(elkocc(ist,ik)-elkocc(jst,ik))/(dble(nkptnr))

                          if (abs(t1) < 1.d-8) cycle

                          do iw=1,nw
                          
                t2=t1/(w(iw)+(eij)+etasm)
                eps(iw,1)=eps(iw,1)+t2*vzzkq(1,1,ist,jst,1,ik)
               
                          end do 
                        end do 
                end do 
         end do 

         eps(:,1)=1.d0-eps(:,1)
         wloss(:,1)=aimag(eps(:,1))/(zabs(eps(:,1))**2)
        
         open(unit=50,file='eps00_Gamma.dat',status='unknown')
        
         do iw=1,nw
        write(50,*)wplt(iw),dble(eps(iw,1)),aimag(eps(iw,1)),wloss(iw,1)
         end do
         write(50,*)

        close (50)

        end subroutine
