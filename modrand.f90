        module modrand
        
! Data related to the k-resolved energy
        type tel
           ! Shift applied to the eigenvalues 
           real(8) :: fe_shift
           ! Number of bands for which the sampling is performed
           integer :: nob
           ! Index of the first band from which the sampling is done
           integer :: n1
           ! Total number of k-points that will be sampled 
           integer :: nktot
           ! Initial k point grid dimension (corresponds to ngridk)
           integer :: nk(3)
           ! Number of random k points generated on a given band
           integer, allocatable :: nkib(:)
           ! K point vectors that are generated 
           real(8), allocatable :: k(:,:)
           ! Initial eigenvalues of the bands of interest
           real(8), dimension(:,:,:,:), allocatable :: emat
           ! Eigenvalues of the random k points
           real(8), dimension(:), allocatable :: energy
           ! Weight of the random k points
           real(8), dimension(:), allocatable :: weight
           ! KS or MB gap of the random k points 
           real(8), dimension(:), allocatable :: Delta
           ! Band associated to a random k points
           integer, dimension(:), allocatable :: istrand
           ! Mapping of krand and k'rand with irreducible q and rotated k
           integer, dimension(:,:,:), allocatable :: ikmap
           ! Mapping of krand, k'rand and g vectors with the rotated g vectors
           integer, dimension(:,:,:), allocatable :: igmap
           ! Occupation using Fermi-Dirac distribution of the k states
           real(8), dimension(:), allocatable :: kocc
        end type tel
        type(tel) :: el

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Data related to the final q grid for the dielectric function
        type tqgr
           ! Dense q-point grid size
           integer nfinq(3)
           ! Number of all dense q-points
           integer nfinqpt
           ! Vector in lattice coordinates of dense q-points
           real(8), allocatable :: vfinql(:,:)
           ! Locations of q-points on integer grid
           integer, allocatable :: ivfinq(:,:)
           ! Map from integer grid to dense q-point index
           integer, allocatable :: ivfinqiq(:,:,:)
           ! Dielectric function computed on the q dense grid
           complex(8), dimension(:,:), allocatable :: eps
           ! Full epsgg(w,q) for istrogram loss function
           complex(8), allocatable :: totaleps(:,:,:)
           end type tqgr
        type(tqgr) :: qgr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Quantities used during the interpolation and control quantities
        type tinp
           ! final k point grid dimension after fft 
           integer :: nr(3) 
           ! random initialization seed of the kp sampling seed 
           integer :: irand 
           ! seed used during the kp sampling 
           integer, allocatable :: iseed(:)
           ! logical for assigning automatically nkib, if false read
           logical :: autokp
           ! parameter used in the distribution for assigning nkib
           real(8) :: ef_window
           ! logical if true the velocity at the random kp are computed
           logical :: fermi_velocity
        end type tinp
        type(tinp) :: inp 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Saxon-Woods function for the random sampling selection
! p(x)=Pmin+(1-Pmin)*(1+exp(-w/s))/(1+exp((|x|-w)/s))
        type tsamp  
           ! Pmin parameter of the Saxon-Woods function
           real(8) :: Pmin
           ! w parameter of the Saxon-Woods function
           real(8) :: width
           ! s parameter of the Saxon-Woods function
           real(8) :: skin 
        end type tsamp
        type(tsamp) :: samp

        end module
