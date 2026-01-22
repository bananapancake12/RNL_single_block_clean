module declaration

  complex                      im                      ! Imaginary unit
  real(8)                      pi                      ! Pi number

  real(8)                      t,Re,mpgx               ! Time, Reynolds number, pressure gradient
  real(8)                      alp,bet,Lx,Ly,Lz        ! Wavelengths in x and z (alpha, beta) and sizes of the box (L = 2pi/alpha)
  integer,      allocatable :: dummint(:)

  integer                      nsamp,isamp             ! Number of files and counter
  integer                      istat,iter              ! Amount of statistics accumulated in the file, iteration in the simulation
  integer                      istatW                  ! Amount of conditional statistics accumulated in file after firstgood
  integer                      firstgood               ! Summing statistics after this number of files/iterations (isamp>=firstgood) 
  integer                      ptstilex,ptstilez       ! Points per tile in x and z (Ngal(1,1)/ntilex) (Ngal(2,1)/ntilez)

  integer                      nband                   ! Number of bands (3)
  integer,      allocatable :: N(:,:),Ngal(:,:)        ! Size of the mesh in spectral and physical space (aliasing)
  integer                      ntilex,ntilez           ! Number of tiles in x and z
  integer                      dnx,dnz!,dny             ! Number of points of a tile in x and z, height of the first band 
  integer                      dnyu,dnyv
  integer                      nx,nz,ny0,nyF           ! Long penciles (second band): N(1,2)/2, N(2,2)/2, height of the box 
  integer                      nyu0,nyuF,nyv0,nyvF           ! Long penciles (second band): N(1,2)/2, N(2,2)/2, height of the box 

  !real(8),      allocatable :: y(:),dy2i(:,:),dthdy(:)
  real(8),      allocatable :: yu(:),dyu2i(:,:),dthdyu(:)
  real(8),      allocatable :: yv(:),dyv2i(:,:),dthdyv(:)
  real(8)                      dtheta,dthetai,ddthetai

  real*8,       allocatable :: spU(:,:,:),spXu(:,:,:)   ! Variable to store spectra, temporal variable
  real*8,       allocatable :: spV(:,:,:),spXv(:,:,:)   ! Variable to store spectra, temporal variable
  real*8,       allocatable :: spP(:,:,:),spXp(:,:,:)   ! Variable to store spectra, temporal variable

  character*80                 fnameima,fnameimb
  character*80                 fildir,filelist,filout  ! Directory of files, name of the file with the list of files
  character*4                  ext1,ext2,ext3          ! Ngalx, Ngaly, Ngalz
  character*3                  extC,extD               ! ntilex, ntilez (for naming files)
  character*5,  allocatable :: extV(:),extL(:)         ! Variables: U,V,W,UV,P. Last number in file names (time): Used to read files.

  real(8),      allocatable :: Xm (:),        XmC (:,:,:)
  real(8),      allocatable :: Um (:),U2m (:),UmC (:,:,:),U2mC (:,:,:)
  real(8),      allocatable :: Vm (:),V2m (:),VmC (:,:,:),V2mC (:,:,:)
  real(8),      allocatable :: Wm (:),W2m (:),WmC (:,:,:),W2mC (:,:,:)
  real(8),      allocatable :: UVm(:),        UVmC(:,:,:)
  real(8),      allocatable :: UWm(:),        UWmC(:,:,:)
  real(8),      allocatable :: VWm(:),        VWmC(:,:,:)
  real(8),      allocatable :: wxm(:),wx2m(:),wxmC(:,:,:),wx2mC(:,:,:)

end module
