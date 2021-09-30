module variables

integer i,j,k,l,n,ii,nn,m,itr,mode,nmode,nbs
integer jmax,kmax,lmax,ary,nc
integer lda,ldvt,ldu,lwork,info

integer :: nbinst, nr, nz, nb
integer :: nbg, interval

real(8) dx,dy,dz,eps,jlength,llength

real(8) dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
!real(4) re,fsmach,alpha


    real(8) :: thetagT0, deltat, r0, rho0, p0, a0, deltat_nondim, theta, pi,uj
    real(8) :: uu, vv, ww, uvw, e
    real(8), allocatable,dimension(:)    :: rg,zg
    real(8), allocatable,dimension(:,:)  ::rhoT0moy,urT0moy,utT0moy,uzT0moy,pT0moy
    real(8), allocatable,dimension(:,:)  :: vortT0moy,dilatT0moy
    real(8), allocatable,dimension(:,:)  :: rhoT0,urT0,utT0,uzT0,pT0,eT0
    real(8), allocatable,dimension(:,:,:)  :: vortT0,dilatT0

    real   :: ds

real(8),allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)
real(4),allocatable :: buf(:,:,:,:),buf_ave(:,:,:)

real(8),allocatable :: rho(:,:),u(:,:),v(:,:),w(:,:),p(:,:),aimp(:,:)

real(8),allocatable :: rho_ave(:,:),u_ave(:,:),v_ave(:,:),w_ave(:,:),e_ave(:,:),p_ave(:,:),vdvd_ave(:,:)

real(8),allocatable :: rhod(:,:),ud(:,:),vd(:,:),wd(:,:),ed(:,:),pd(:,:),vdvd(:,:),kinetic_e(:,:),p_total(:,:),tke(:,:)

real(8),allocatable :: ptotal_all(:,:),ptotal_lami(:,:),ptotal_turb(:,:),aimp_all(:,:),aimp_lami(:,:),aimp_turb(:,:)


real(8),allocatable :: dpdx(:,:), d2pdx2(:,:)

real(8),allocatable :: vortxd(:,:),vortyd(:,:),vortzd(:,:)

real(8),allocatable :: q(:,:,:,:)

   data rho0,a0,p0,theta,interval/1.205,343,1.0d5,0d0,1/

character filename*256,flowinst*128,flowave*128,gridname*128,outfile*128



end module
