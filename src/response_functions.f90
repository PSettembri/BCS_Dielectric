subroutine response_functions
use modmain, only: T 
use modrand, only: el 
!el%Delta =  superconducting gap on random points
!el% 
!use mod_input, only: func
implicit none
integer :: ind
integer :: iwn
integer :: indp,iw,ind_gkkp,ind_freq,iep,ie,iomega,indq
real(8) :: wn,ei,eip,Z,D,ZN, EE, EE0, domega, EEs, EEsp
real(8), allocatable :: uk(:),vk(:),uk2(:),vk2(:),ukvk(:),fe(:),nu2(:),chi0iso(:,:,:)
real(8), allocatable :: fEN(:),chi0Niso(:,:)
complex(8), allocatable :: om(:),chi0NisoR(:,:),chi0isoR(:,:,:)
complex(8) :: ii
real(8) :: emax_isto,emin_isto,de,Zsum,NORM,tmp,Rk(3),numax,nu,eta,q(3)
real(8) :: Mkkp
real(8) :: ww,e12,e21,f12,f21,xxx,xxx0
integer :: npoints_isto, i, iRot, nk(3),ib,nb,ik1,ik2,ik3,nb1,nb2,ik1_,ik2_,ik3_
integer :: nnu,inu,nqiso,np2,ngridq(3)
integer, allocatable :: indqmat(:,:)
logical :: kp_are_there
real(8) :: fermi2 ! function
ii=cmplx(0.0d0,1d0)
xxx0=0.00
write(*,*) 'calc Response'
flush(6)

!playing with data
!el%Delta=0.000d0
!where(el%ekf(:).gt.0.0d0) el%ekf(:)=0.005 !introducing a gap
!where(el%ekf(:).lt.0.0d0) el%ekf(:)=-0.005 !introducing a gap

! Preliminary
write(*,*) 'preparing uk, vk and fk'
open(5000,file="uk2vk2.dat")
open(5001,file="ukvk.dat")
allocate(uk2(el%npoints),vk2(el%npoints),uk(el%npoints),vk(el%npoints),ukvk(el%npoints))
do ind=1,el%npoints
   ie=el%ind2imesh(ind)
   !ei=ene1D%emesh(ie)
   ei=el%ekf(ind)
   EEs=sqrt(ei**2 + el%Delta(ind)**2) 
   !EE=sqrt(eip**2 + func%gamma3_SE*el%Delta(indp)**2)  ! adjusted as in the SE code
   uk2(ind)=(1d0+ei/EEs)/2 ; uk(ind)=sqrt(uk2(ind))
   vk2(ind)=(1d0-ei/EEs)/2 ; vk(ind)=sqrt(vk2(ind))
   ukvk(ind)=abs(el%Delta(ind))/EEs/2
   write(5000,'(4f12.6)') ei, vk2(ind),uk2(ind)
   write(5001,'(4f12.6)') ei, vk(ind)*uk(ind),ukvk(ind)
enddo
close(5000)
close(5001)

open(5000,file='fermifunctions.dat')
allocate (fE(-el%npoints:el%npoints),fEN(el%npoints))
do ind=1,el%npoints
   ie=el%ind2imesh(ind)
   !ei=ene1D%emesh(ie)
   ei=el%ekf(ind)
   EEs=sqrt(ei**2 + el%Delta(ind)**2)
   fE(ind)=fermiII(EEs) 
   fE(-ind)=fermiII(-EEs)
   feN(ind)=fermiII(ei)
   write(5000,*) ei,feN(ind)
enddo
close(5000)


! preparing frequencies
!chi0=0.0d0
!nq(:)=100 
nnu=50
numax=0.04d0
eta=0.0001d0
write(*,*) 'kbT=', 1d0/T%beta 
write(*,*) 'eta=', eta
eta=20d0/T%beta  !  a bit high but is multiplied by nu/numax
write(*,*) 'fixed eta=', eta
write(*,*) 'd_omega=',numax/float(nnu)

allocate(nu2(0:nnu)) 
allocate(om(0:nnu)) 
do inu=0,nnu
   nu=float(inu)/float(nnu)*numax
   if (nu.eq.0.0d0) nu=1d-14
   nu2(inu)=nu**2
   !om(inu)=nu + ii*eta
   ! frequency dependent nu (kind of a rotation on the imaginary plane)
   om(inu)=nu + ii*eta*nu/numax
enddo


ngridq=1
ngridq(1)=10
allocate(indqmat(el%npoints,el%npoints))
indqmat=0

! Setting q-grid
open(5000,file='testq.dat')
do ind=1,el%npoints
   ei=el%ekf(ind)
   if (abs(ei).gt.2*numax) cycle
   EEs=sqrt(ei**2 + el%Delta(ind)**2)
   do indp=1,el%npoints
      eip=el%ekf(indp)
      if (abs(eip).gt.2*numax) cycle
      EEsp=sqrt(eip**2 + el%Delta(indp)**2)
      np2=np2+1 ! counter
      q(:)=el%k(ind,:)-el%k(indp,:)
      call bring2gammac(q(:))
      if (abs(q(3)).gt.0.001.or.abs(q(2)).gt.0.001) cycle
      if (abs(q(3)).lt.0.001) write(5000,'(3f10.5)')  q(1),q(2),q(3)
      call bring2ucII(q(:))
      call vect2grid(q,ngridq,indq)
      !write(*,*) indq,q(1)
      indqmat(ind,indp)=indq
   enddo
enddo
close(5000)
write(*,*) np2, el%npoints*el%npoints, el%npoints

!allocate(chi0iso(10,0:nnu,4)
chi0iso=0.0d0
allocate(chi0isoR(10,0:nnu,8)) 
allocate(chi0Niso(10,0:nnu)) 
allocate(chi0NisoR(10,0:nnu))

chi0isoR=0.0d0
chi0Niso=0.0d0
chi0NisoR=0.0d0
nqiso=1
!allocate(chi0(nq(1),nq(2),nq(3),nnu) 
do ind=1,el%npoints
   ei=el%energy(ind)
   if (abs(ei).gt.2*numax) cycle
   ist = el%istrand(ind)
   occi= el%kocc(ind)
   EEs=sqrt(ei**2 + el%Delta(ind)**2)
   do indp=1,el%npoints
      eip=el%energy(indp)
      if (abs(eip).gt.2*numax) cycle
      jst=el%istrand(eip)
      occj=el%kocc(eip)
      ww=el%weight(ind)*el%weight(indp)

      !!!!!!!!!!!!!!WARNING

      iq=el%ikmap(ind,indp,1)
      ik=el%ikmap(ind,indp,2)

      ! qrand = k'rand-krand
      qrand(:) = el%k(indp,:)-el%k(ind,:)

      iq1=modulo(nint(qrand(1)*qgr%nfinq(1)),qgr%nfinq(1))
      iq2=modulo(nint(qrand(2)*qgr%nfinq(2)),qgr%nfinq(2))
      iq3=modulo(nint(qrand(3)*qgr%nfinq(3)),qgr%nfinq(3))

      !iqdens=qgr%ivfinqiq(iq1,iq2,iq3)

      Mkkp=vzzkq(ig,jg,ist,jst,iq,ik)

      !WARNING END


      np2=np2+1 ! counter
      EEsp=sqrt(eip**2 + el%Delta(indp)**2)
      f12=fEN(ind)-feN(indp)
      e21=eip-ei
      ww=el%we(ind)*el%we(indp)
      !chi0iso(1,:,1) = chi0iso(1,:,1) + (EEs-EEsp)*(fE(ind)-fE(indp))/(nu2(:)+(EEs-EEsp)**2)*  & 
      !                                  ww *  uk2(ind)*uk2(indp)
      !chi0iso(1,:,2) = chi0iso(1,:,2) + (-EEs+EEsp)*(fE(-ind)-fE(-indp))/(nu2(:)+(-EEs+EEsp)**2)*  & 
      !                                  ww *  vk2(ind)*vk2(indp)
      !chi0iso(1,:,3) = chi0iso(1,:,3) + (-EEs-EEsp)*(fE(-ind)-fE(indp))/(nu2(:)+(-EEs-EEsp)**2)*  & 
      !                                  ww *  vk2(ind)*uk2(indp)
      !chi0iso(1,:,4) = chi0iso(1,:,4) + (EEs+EEsp)*(fE(ind)-fE(-indp))/(nu2(:)+(EEs+EEsp)**2)*  & 
      !                                  ww *  uk2(ind)*vk2(indp)
      indq=indqmat(ind,indp)
      if (indq.eq.0) cycle
      if (indq.gt.10) stop 'boh'
 
      chi0isoR(indq,:,1) = chi0isoR(indq,:,1) +Mkkp* (fE(ind)-fE(indp))/(om(:)+(EEs-EEsp))*  & 
                                        ww *  ( uk2(ind)*uk2(indp) )!  +   uk(ind)*vk(ind)*vk(indp)*uk(indp) ) 
      chi0isoR(indq,:,2) = chi0isoR(indq,:,2) +Mkkp* (fE(-ind)-fE(-indp))/(om(:)+(-EEs+EEsp))*  & 
                                        ww *  ( vk2(ind)*vk2(indp) )!  +   vk(ind)*uk(ind)*uk(indp)*vk(indp) )
      chi0isoR(indq,:,3) = chi0isoR(indq,:,3) +Mkkp* (fE(-ind)-fE(indp))/(om(:)+(-EEs-EEsp))*  & 
                                        ww *  ( vk2(ind)*uk2(indp) )!  -   vk(ind)*uk(ind)*uk(indp)*vk(indp) )
      chi0isoR(indq,:,4) = chi0isoR(indq,:,4) +Mkkp* (fE(ind)-fE(-indp))/(om(:)+(EEs+EEsp))*  & 
                                        ww *  ( uk2(ind)*vk2(indp) )!  -   uk(ind)*vk(ind)*uk(indp)*vk(indp) )
      chi0isoR(indq,:,5) = chi0isoR(indq,:,5) +Mkkp* (fE(ind)-fE(indp))/(om(:)+(EEs-EEsp))*  & 
                                        ww *  ( uk(ind)*vk(ind)*vk(indp)*uk(indp) ) 
      chi0isoR(indq,:,6) = chi0isoR(indq,:,6) +Mkkp* (fE(-ind)-fE(-indp))/(om(:)+(-EEs+EEsp))*  & 
                                        ww *  ( vk(ind)*uk(ind)*uk(indp)*vk(indp) )
      chi0isoR(indq,:,7) = chi0isoR(indq,:,7) -Mkkp* (fE(-ind)-fE(indp))/(om(:)+(-EEs-EEsp))*  & 
                                        ww *  ( vk(ind)*uk(ind)*uk(indp)*vk(indp) )
      chi0isoR(indq,:,8) = chi0isoR(indq,:,8) -Mkkp* (fE(ind)-fE(-indp))/(om(:)+(EEs+EEsp))*  & 
                                        ww *  ( uk(ind)*vk(ind)*uk(indp)*vk(indp) )

      chi0Niso(1,:) = chi0Niso(1,:)  -  f12 * e21 / ( nu2(:)+e21**2 ) *ww
      chi0NisoR(1,:)= chi0NisoR(1,:) -  f12       / ( om(:)+e21 )     *ww
   enddo
enddo    

chi0isoR(:,:,5)=real(chi0isoR(:,:,5))
chi0isoR(:,:,6)=real(chi0isoR(:,:,6))
chi0isoR(:,:,7)=real(chi0isoR(:,:,7))
chi0isoR(:,:,8)=real(chi0isoR(:,:,8))

write(*,*) np2,'np2'

open(5000,file='chi0iso.dat')
open(5001,file='chi0isoR.dat')
open(5002,file='chi0Niso.dat')
open(5003,file='chi0NisoR.dat')
open(5004,file='chi0isoR_parts.dat')
do inu=0,nnu
   !write(5000,'(10f18.6)') sqrt(nu2(inu)) , chi0iso(1,inu,1),chi0iso(1,inu,2),chi0iso(1,inu,3),chi0iso(1,inu,4)
   write(5001,'(30f18.6)') sqrt(nu2(inu)) , (sum(chi0isoR(indq,inu,:)),indq=1,10)
   write(5004,'(30f18.6)') sqrt(nu2(inu)) , chi0isoR(1,inu,1),chi0isoR(1,inu,2),chi0isoR(1,inu,3),chi0isoR(1,inu,4) &
                                          , chi0isoR(1,inu,5),chi0isoR(1,inu,6),chi0isoR(1,inu,7),chi0isoR(1,inu,8) 
   write(5002,'(10f18.6)') sqrt(nu2(inu)) , chi0Niso(1,inu)
   write(5003,'(10f18.6)') sqrt(nu2(inu)) , chi0NisoR(1,inu)
enddo
close(5000)
close(5001)
close(5002)
close(5003)
close(5004)


write(*,*) 'END'

end

! ---------------------------------------------------------

real(8) function fermiII(xxx)
use modmain, only : beta 

implicit none
real(8) :: xxx,arg

arg=beta*xxx
call adjust(arg) 

fermiII=1.d0/(dexp(arg)+1.d0)

end function fermiII

! ---------------------------------------------------------

real(8) function boseII(xxx)
use modmain, only : beta
implicit none
real(8) :: xxx,tol2,arg
parameter(tol2=1.d-14)

if (dabs(xxx).lt.tol2) xxx=dsign(1.d0,xxx)*tol2
arg=beta*xxx
call  adjust(arg)

boseII=1.d0/(dexp(arg)-1.d0)

!----------------------------------------------------

subroutine adjust(ene)
use modmain, only: small, big
real(8) ::  ene
if (dabs(ene).lt.small) ene = small * dsign(1d0,ene)
if (dabs(ene).gt.big)   ene = big   * dsign(1d0,ene)
end

!----------------------------------------------------


