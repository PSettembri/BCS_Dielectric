        module modmain
        implicit none

        ! General variables
        ! starting time
        real (4) tstart
        ! ending time
        real (4) tend

        ! Lattice parameters
        ! lattice vectors from ELK (atomic units) first index is the
        ! coordinate x,y,z , second index is vector index a1,a2,a3
        real (8) avec(3,3)
        ! reciprocal lattice vectors from ELK (atomic units)
        real(8) bvec(3,3)
        ! unit cell volume
        real(8) omega
        ! Brillouin zone volume
        real(8) omegabz
        ! Fermi energy of ELK
        real(8) elkfermi
        ! ELK gap estimate
        real(8) elkgap
        ! Fermi-Dirac smearing
        real(8) smear
        ! Eigenvalues from ELK
        real (8), allocatable :: elkeig(:,:)
        ! Occupancies from ELK
        real (8), allocatable :: elkocc(:,:)
        ! Scissor operator
        real (8) :: scissor

        ! K & Q points
        ! number of irreducible k-points
        integer nkpt
        ! number of reducible k-points
        integer nkptnr
        ! vector in lattice coordinates of reducible k-points
        real(8), allocatable :: vkl(:,:)
        ! number of irreducible q points
        integer nqpt 
        ! number of reducible q-points
        integer nqptnr        
        ! vector in lattice coordinates of reducible k-points
        real(8), allocatable :: vql(:,:)
        ! mapping from qirr to q 
        integer, allocatable :: mapqirrtoq(:)


        ! Response function parameters
        
        ! G vectors
        ! |G| cut-off for response functions
        real(8) gmaxrf
        ! number of G-vectors for response functions
        integer ngrf
        ! G vectors
        integer, allocatable :: ivg(:,:) 
        ! Map of rotated G vector according to symmetries
        integer, allocatable :: rotg(:,:)

        ! States
        ! minimum state index considered for response functions
        integer nstminrf
        ! maximum state index considered for response functions
        integer nstmaxrf
        ! total number of states used in ELK calculation
        integer nstsv
        ! number of considered states
        integer nstates
        
        ! Frequencies
        ! minimal frequency computed for response functions
        real (8) wminrf 
        ! maximal frequency computed for response functions
        real (8) wmaxrf
        ! step between different frequencies considered
        real (8) dwrf
        ! number of frequencies considered
        integer nw
        ! array of frequencies w
        real (8), allocatable :: w(:)
        ! array of frequencies w in eV
        real (8), allocatable :: wplt(:)
        ! complex smearing of the response function
        real (8) eta
        ! Temperature, used for response Fermi functions, Ideally should be consistent with Delta(T)
        real (8) T,beta,kb
        ! other useful parameters
        real(8) small,small2,big
        ! Parameters
        parameter ( kb    = 6.3337008d-6          )  ! Ry/K
        parameter ( small = 5d-15                 )
        parameter ( small2 = 1d-12                )
        parameter ( big   = 60d0                  )


        !Symmetries
        ! Total number of crystal symmetries
        integer nsymcrys
        ! Traslation vector of each symmetry operation
        real (8), allocatable :: vtlsymc(:,:)
        ! Spacial rotation of each symmetry operation
        integer, allocatable :: symlatspa(:,:,:)
        ! Global spin rotation of each symmetry operation
        integer, allocatable :: symlatspin(:,:,:)
        ! K-point grid sizes
        integer ngridk(3)
        ! Q-point grid sizes
        integer ngridq(3)
        ! K-point grid offset
        real(8) vkloff(3)
        ! any vector with length less than epslat is considered zero
        real(8) epslat
        ! locations of k-points on integer grid
        integer, allocatable :: ivk(:,:)
        ! locations of q-points on integer grid
        integer, allocatable :: ivq(:,:)        
        ! map from integer grid to reducible k-point index
        integer, allocatable :: ivkik(:,:,:)
        ! map from integer grid to reducible q-point index
        integer, allocatable :: ivqiq(:,:,:)
        ! map from integer grid to non-reduced k-point index
        integer, allocatable :: ivkiknr(:,:,:)
        ! map from integer grid to non-reduced q-point index
        integer, allocatable :: ivqiqnr(:,:,:)

        ! vzz matrix elements
        complex(8), allocatable :: vzzkq(:,:,:,:,:,:)
        ! complex dielectric function 
        complex(8), allocatable :: eps(:,:)
        ! loss function 
        real (8), allocatable :: wloss(:,:)
        ! full epsgg(w,q) for istrogram loss function
        complex(8), allocatable :: totaleps(:,:,:)
        ! Istogram maximum k value
        real (8) kisto
        ! Spacing of k in loss function istogram
        real (8) dbin
        ! Logical for naming output files
        logical :: invers
        ! Logical if the elk and code epsilon are combined
        logical :: combination

        ! Tasks management
        ! maximum number of tasks
        integer, parameter :: maxtasks=40
        ! number of tasks
        integer ntasks
        ! task index
        integer itask
        ! task array
        integer tasks(maxtasks)
        ! current task
        integer task


        ! Conversion parameters
        ! Hartree in eV (CODATA 2018)
        real(8), parameter :: ha_ev = 27.211386245988d0
        ! Recuded Planck constant times c speed of light [KeV*A] 
        real(8), parameter :: hbarc_ev = 1.973269631d0
        ! Bohr radius in Angstrom (CODATA 2018)
        real(8), parameter :: br_ang=0.529177210903d0
        ! keV to eV conversion
        real(8), parameter :: kev_ev=1000.0d0
        ! pi 
        real(8), parameter :: pi= 3.14159265358979323844d0
        ! imaginary unit
        complex(8), parameter :: icomp = (0d0,1d0) 

        end module
