! Calculate the steady-state of Langmuir waves considering the
! effects of electrostatic bremsstrahlung and collisional damping
! for Maxwellian electron velocity distribution function
! (see https://doi.org/10.1088/1361-6587/ab4aad).
! Written by Sabrina F. Tigik.

module common_params
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp) :: Qxi=1.e-4_wp, Qxf=0.6_wp
  real(wp) :: Qzi=1.e-4_wp, Qzf=0.6_wp
  real(wp) :: Qi= 1.e-4_wp, Qf= 10._wp ! For VQQ
  real(wp) :: RTeTi= 4._wp, G=5.e-3_wp
  real(wp) :: Ve2C2= 4.e-3_wp
  real(wp) :: RMeMi,RMiMe
  real(wp) :: Geff,AA
  integer :: nqx=128, nqz=128
  integer :: nqcd=2048,nrz=1
end module common_params

module common_arrays
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp), dimension(:), allocatable :: VQ,VQQ
  real(wp), dimension(:), allocatable :: VQx
  real(wp), dimension(:), allocatable :: VQz
  real(wp), dimension(:,:), allocatable :: IL,IL0
  real(wp), dimension(:), allocatable :: IL1D,IL01D
  real(wp), dimension(:,:), allocatable :: GcollLp
  real(wp), dimension(:,:), allocatable :: GcollLm
  ! real(wp), dimension(:,:), allocatable :: GcollL 
  real(wp), dimension(:), allocatable:: GcollL1D
  real(wp), dimension(:), allocatable:: GcollLp1D,GcollLm1D
  ! real(wp), dimension(:,:), allocatable :: BremL
  real(wp), dimension(:,:), allocatable :: BremLp,BremLm
  real(wp), dimension(:), allocatable :: BremL1D
  real(wp), dimension(:,:), allocatable :: BrGcL
  real(wp), dimension(:), allocatable :: BrGcL1D
  ! real(wp), dimension(:,:), allocatable :: GqlLp,GqlLm
  ! real(wp), dimension(:,:), allocatable :: GqlL1D
  ! real(wp), dimension(:,:), allocatable:: GcollSp,GqlSp
  ! real(wp), dimension(:,:), allocatable:: GcollSm,GqlSm
  ! real(wp), dimension(:), allocatable:: GcollS1D,GqlS1D,BremS1D
  ! real(wp), dimension(:), allocatable:: GcollSp1D,GqlSp1D
  ! real(wp), dimension(:), allocatable:: GcollSm1D,GqlSm1D
  ! real(wp), dimension(:,:), allocatable :: BremSp,BremSm,BremS
  real(wp), dimension(:), allocatable :: VRNeNs,VRTeTs,VRTiTs
  ! Auxiliary for collisional damping:
  real(wp), dimension(:), allocatable :: Aux1_Gcoll
  ! integer, dimension(:), allocatable :: Aux2_Gcoll
  ! Auxiliary for electrostatic bremsstrahlung:
  real(wp), dimension(:), allocatable :: Aux1_Bremss
  ! integer, dimension(:), allocatable :: Aux2_Bremss
end module common_arrays

module math_constants 
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp) :: Pi = 4._wp*atan(1._wp)
  real(wp) :: Sqtwo = sqrt(2._wp)
  ! real(wp) :: Pi= 3.1415926535897932384626433832795_wp
  ! real(wp) :: Sqtwo= 1.41421356237309504880168872420969_wp
  real(wp) :: Infinity = 1.0e+30_wp
  real(wp) :: Degree = 0.01745329251994329576923690768488_wp
  real(wp) :: EpsMin = 1.0e-16_wp
  real(wp) :: Xacc = 1.0e-6_wp
  real(wp) :: Qmin = 1.0e-4_wp
  real(wp) :: Dqaux = 1.0e-2_wp
  complex(wp) :: Zi = (0._wp,1._wp)
  complex(wp) :: Zzero = (0._wp,0._wp)
end module math_constants

module phys_constants
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp) :: Me = 9.10938356e-28_wp       ! electron mass (g)
  real(wp) :: Mi = 1.6726219e-24_wp        ! proton mass (g)
  real(wp) :: MeC2 = 510.998946e+3_wp      ! electron mass (eV)
  real(wp) :: MpC2 = 938.272081e+6_wp      ! proton mass (eV)
  real(wp) :: C_SI = 2.997925e+8_wp        ! speed of light (m/s)
  real(wp) :: C_cgs = 2.997925e+10_wp      ! speed of light (cm/s)  
end module phys_constants

module sub_prog
contains

  subroutine allocate_arrays
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use common_params
    use common_arrays
    implicit none
    allocate (VQ(nrz))
    allocate (VQx(nqx))
    allocate (VQz(nqz))
    allocate (VQQ(nqcd))
    allocate (IL0(nqx,nqz),IL(nqx,nqz))
    allocate (IL01D(nqcd),IL1D(nqcd))
    allocate (GcollLp(nqx,nqz),GcollLm(nqx,nqz))
    ! allocate (GcollL(nqx,nqz))
    allocate (GcollL1D(nqcd))
    allocate (GcollLp1D(nqcd),GcollLm1D(nqcd))
    ! allocate (BremL(nqx,nqz))
    allocate (BremL1D(nqcd))
    allocate (BremLp(nqx,nqz))
    allocate (BremLm(nqx,nqz))
    allocate (BrGcL(nqx,nqz))
    allocate (BrGcL1D(nqcd))
    ! allocate (GqlLp(nqx,nqz),GqlLm(nqx,nqz))
    ! allocate (GqlL1D(nqcd))
    allocate (VRNeNs(nrz),VRTeTs(nrz),VRTiTs(nrz))
    ! Auxiliary for collisional damping:
    allocate (Aux1_Gcoll(3))
    ! allocate (Aux2_Gcoll(2))
    ! Auxiliary for electrostatic bremsstrahlung:
    allocate(Aux1_bremss(3))
    ! allocate(Aux2_bremss(2))
  end subroutine allocate_arrays

  subroutine Definitions
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use common_params
    use common_arrays
    use math_constants
    use phys_constants
    implicit none
    real(wp) :: RMpMe,RMeMp
    integer :: m

    m=1
    RMpMe= MpC2/MeC2
    RMeMp= MeC2/MpC2
    RMiMe= RMpMe
    RMeMi= 1._wp/RMiMe
    AA= sqrt(1._wp+3._wp/RTeTi)/sqrt(RMiMe)/sqrt(2._wp)
    !Geff= G/(2._wp*sqrt(2._wp)*(4._wp*Pi)**2)
    VRNeNs(m)= 1._wp
    VRTeTs(m)= 1._wp
    VRTiTs(m)= 1._wp/RTeTi
  end subroutine Definitions

  subroutine init_vec
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use common_params
    use common_arrays
    implicit none
    real(wp) :: Dqx,Dqz,Dq
    integer :: i
    
    Dqx= (Qxf-Qxi)/(nqx-1)
    do i= 1,nqx
       VQx(i)= Qxi+(i-1)*Dqx
    end do
    VQx(nqx)= Qxf
    Dqz= (Qzf-Qzi)/(nqz-1)
    do i= 1,nqz
       VQz(i)= Qzi+(i-1)*Dqz
    end do
    VQz(nqz)= Qzf

    Dq= (Qf-Qi)/(nqcd-1)
    do i= 1,nqcd
       VQQ(i)= Qi+(i-1)*Dq
    end do
    VQQ(nqcd)= Qf
  end subroutine init_vec

  subroutine wave_init
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use common_arrays
    use common_params
    use math_constants
    implicit none
    real(wp) :: Q,Q2
    real(wp) :: Qx,Qz
    real(wp) :: Zlq,Zlq2,Aux
    integer :: i,j,m,res

    m = 1
    Geff= G/(2._wp*sqrt(2._wp)*(4._wp*pi)**2)
    call coll_damping
    call bremsstrahlung

    do i=1,nqcd
       Q= VQQ(i)
       Q2=Q*Q
       Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
       Zlq2= Zlq*Zlq
       IL01D(i)= Geff*VRNeNs(m)*VRTeTs(m)/2._wp/Zlq2
       IL1D(i)= (sqrt(Pi)*Geff*VRNeNs(m)**2/sqrt(VRTeTs(m)) &
            * exp(-Zlq2/Q2/VRTeTs(m))/Q2/Q + BremL1D(i)) &
            / (2._wp*sqrt(Pi)*VRNeNS(m)/sqrt(VRTeTs(m)**3)*(Zlq2) &
            * exp(-Zlq2/Q2/VRTeTs(m))/Q2/Q - 2._wp*GcollL1D(i))
       BrGcL1D(i)= - BremL1D(i)/GcollL1D(i)
    end do
    
    Q= 0._wp
    Q2= 0._wp
    Zlq= 0._wp
    Zlq2= 0._wp
    do i=1,nqx
       Qx=VQx(i)
       do j=1,nqz
          Qz=VQz(j)
          Q2=Qx**2+Qz**2
          Q=sqrt(Q2)
          ! if(Q<=5.e-3_wp) Q=5.e-3_wp
          ! call Locate(VQQ,nqcd,Q,ires)
          ! call Aitp1d2(nqcd,VQQ,IL01D,Q,Aux,ires)
          ! IL0(i,j)=Aux
          ! call Locate(VQQ,nqcd,Q,ires)
          ! call Aitp1d2(nqcd,VQQ,IL1D,Q,Aux,ires)
          ! IL(i,j)=Aux
          Zlq= ZL(Qx,Qz)
          Zlq2= Zlq**2
          IL0(i,j)=Geff*VRNeNs(m)*VRTeTs(m)/2._wp/Zlq2
          IL(i,j)= (sqrt(Pi)*Geff*VRNeNs(m)**2/sqrt(VRTeTs(m)) &
               * exp(-Zlq2/Q2/VRTeTs(m))/Q2/Q + BremLp(i,j)) &
               / (2._wp*sqrt(Pi)*VRNeNS(m)/sqrt(VRTeTs(m)**3)*(Zlq2) &
               * exp(-Zlq2/Q2/VRTeTs(m))/Q2/Q - 2._wp*GcollLp(i,j))
          BrGcL(i,j)= -BremLp(i,j)/GcollLp(i,j)
       end do
    end do
    call output("GcL")
    call output("BrL")
    call output("IL0")
    call output("BGL")
    return
  end subroutine wave_init

  real(wp) function ZL(Qx,Qz)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    implicit none
    real(wp), intent(in) :: Qx,Qz
    real(wp) :: Q2
    integer :: m

    m= 1
    Q2= Qx**2+Qz**2
    ZL= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
    if (Qz < 0._wp) then
       ZL= -ZL
    else
    end if
    return
  end function ZL
  
  subroutine coll_damping
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use common_params
    use common_arrays
    use math_constants
    use phys_constants
    use DQsimp_m
    use zetaFunc_m
    implicit none
    real(wp) :: Q,Q2
    real(wp) :: Qx,Qz
    real(wp) :: Zlq
    real(wp) :: AuxSig,Aux
    real(wp) :: Res,Res1L,Res2L
    integer :: i,j,m
    integer :: ires,sigpm
    
    m= 1
    Geff= G/(2._wp*sqrt(2._wp)*(4._wp*Pi)**2)
    GcollL1D= 0._wp
    GcollLp1D= 0._wp
    GcollLm1D= 0._wp
    ! AuxSig= 0._wp
    do i=1,nqcd
       Q= VQQ(i)
       Q2= Q*Q
       Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
       Aux1_Gcoll(1)= Q
       Aux1_Gcoll(2)= Q2
       Aux1_Gcoll(3)= Zlq
       call dqsimp(Aux_GcollL,1.e-4_wp,4._wp,res1L)
       call dqsimpb(Aux_GcollL,4._wp,1.e+30_wp,res2L)
       res= Res1L+Res2L
       GcollL1D(i)=-8._wp*sqrt(Pi**3)*VRNeNs(m)**4/sqrt(VRTeTs(m)**7) &
               * Geff/Zlq**2/Q2/Q2*Res
       ! do sigpm= -1,1,2
       !    AuxSig=(dfloat(sigpm)*Zlq)
       !    GcollL1D(i)=-8._wp*sqrt(Pi**3)*VRNeNs(m)**4/sqrt(VRTeTs(m)**5) &
       !         * Geff/Zlq2/Q2/Q2*Res
       !    if (AuxSig > 0._wp) then
       !       GcollLp1D(i)= GcollL1D(i)
       !    else
       !       GcollLm1D(i)= GcollL1D(i)
       !    end if
       ! end do
    end do
    
    do i= 1,nqx
       Qx= VQx(i)
       do j= 1,nqz
          Qz= VQz(j)
          Q2= Qx**2+Qz**2
          Q= sqrt(Q2)
          if(Q<=5.e-3_wp) Q=5.e-3_wp
          call Locate(VQQ,nqcd,Q,ires)
          call Aitp1d2(nqcd,VQQ,GcollL1D,Q,Aux,ires)
          GcollLp(i,j)= Aux
          GcollLm(i,j)= Aux
          ! call Aitp1d2(nqcd,VQQ,GcollL1D,Q,Aux,ires)
          ! GcollLm(i,j)= Aux
       end do
    end do
    return
  end subroutine coll_damping

  real(wp) function Aux_GcollL(Qp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use DQsimp_m
    use zetaFunc_m
    implicit none
    real(wp), intent(in) :: Qp
    real(wp) :: Q,Q2,Qp2
    real(wp) :: Zlq,Aux
    real(wp) :: Qp2EpsQpSq
    complex(wp) :: Qp2EpsQp
    complex(wp) :: Zeta
    integer :: m
    
    m= 1 
    Q= Aux1_Gcoll(1)
    Q2= Aux1_Gcoll(2)
    Zlq= Aux1_Gcoll(3) 
    if(Q > 2.e-4_wp) then
       if(Qp > 2.e-4_wp) then
          Qp2= Qp**2
          Zeta= Zlq/Q/Qp/sqrt(VRTeTs(m))
          Qp2EpsQp= 1._wp+2._wp*VRNeNs(m)/VRTeTs(m)/Q2/Qp2 &
               * (1._wp+Zeta*Zfn(Zeta))
          Qp2EpsQpSq= (abs(Qp2EpsQp))**2
          Aux= ((1._wp+Qp2*Qp2)/(1._wp-Qp2)**2+(1._wp+Qp2)/4._wp/Qp &
               * log((1._wp-Qp)**2/(1._wp+Qp)**2))/Qp2/Qp
          Aux_GcollL= exp(-Zlq**2/Qp2/Q2/VRTeTs(m))*Aux/Qp2EpsQpSq
          ! Aux= ((1._wp+Qp2*Qp2)/(1._wp-Qp2)**2+(1._wp+Qp2)/4._wp/Qp &
          !      * log((1._wp-Qp)**2/(1._wp+Qp)**2))/Qp2
          ! Aux_GcollL= exp(-Zlq**2/Qp2/Q2/VRTeTs(m))*Aux/Qp2EpsQpSq
       else
          Aux_GcollL= 0._wp
       end if
    else
       Aux_GcollL= 0._wp
    end if
    return
  end function Aux_GcollL

  subroutine Bremsstrahlung
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    use DQsimp_m
    implicit none
    real(wp) :: Qx,Qz
    real(wp) :: Q,Q2
    real(wp) :: Zlq
    real(wp) :: Aux,ResL
    real(wp) :: Res1L,Res2L
    integer :: i,j
    integer :: m,ires

    m= 1
    Geff= G/(2._wp*sqrt(2._wp)*(4._wp*Pi)**2)
    BremL1D= 0._wp
    do i= 1,nqcd
       Q= VQQ(i)
       Q2= Q**2
       Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
       Aux1_Bremss(1)= Q
       Aux1_Bremss(2)= Q2
       Aux1_Bremss(3)= Zlq
       call DQsimp(Aux_BremL,1.e-4_wp,4._wp,Res1L)
       call DQsimpb(Aux_BremL,4._wp,Infinity,Res2L)
       ResL=Res1L+Res2L
       BremL1D(i)= 24.0_wp*sqrt(Pi**3)/Zlq**2*VRNeNs(m)**4 &
            / sqrt(VRTeTs(m)**5)*Geff**2/Q2/Q2*ResL
    end do

    do i= 1,nqx
       Qx= VQx(i)
       do j= 1,nqz
          Qz= VQz(j)
          Q2= Qx**2+Qz**2
          Q= sqrt(Q2)
          if(Q<=5.e-3_wp) Q=5.e-3_wp
          call Locate(VQQ,nqcd,Q,ires)
          call Aitp1d2(nqcd,VQQ,BremL1D,Q,Aux,ires)
          BremLp(i,j)= Aux
          BremLm(i,j)= Aux
       end do
    end do
    return
  end subroutine Bremsstrahlung
  
  real(wp) function Aux_BremL(Qp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use DQsimp_m
    use zetaFunc_m
    implicit none
    real(wp), intent(in) :: Qp
    real(wp) :: Q,Q2,Qp2
    real(wp) :: Qp2EpsQpSq
    real(wp) :: Zlq,Aux
    complex(wp) :: Qp2EpsQp
    complex(wp) :: Zeta
    integer :: m

    m= 1
    Q= Aux1_Bremss(1)
    Q2= Aux1_Bremss(2)
    Zlq= Aux1_Bremss(3)
    if(Q > 1.e-4_wp) then
       if(Qp > 1.e-4_wp) then
          Qp2= Qp**2
          Zeta= Zlq/Q/Qp/sqrt(VRTeTs(m))
          Qp2EpsQp= 1._wp+2._wp*VRNeNs(m)/VRTeTs(m)/Qp2/Q2*(1._wp+Zeta*Zfn(Zeta))
          Qp2EpsQpSq= (abs(Qp2EpsQp))**2
          Aux= ((1._wp+Qp2*Qp2)/(1._wp-Qp2)**2+(1._wp+Qp2)/4._wp/Qp &
               * log((1._wp-Qp)**2/(1._wp+Qp)**2))/Qp2/Qp
          Aux_BremL= exp(-Zlq**2/Qp2/Q2/VRTeTs(m))*Aux/Qp2EpsQpSq
       else
          Aux_BremL= 0._wp
       end if
    else
       Aux_BremL= 0._wp
    end if    
    return
  end function Aux_BremL
  
  subroutine Aitp1d2(nx,Vx,Fx,Xp,Fp,i)
    ! Uses linear interpolation in order to obtain the value of a function
    ! F, at the point Xp. The function F is given as a set of points Fx(x).
    ! Uses subroutine Locate (Numerical Recipes, P. 96)
    ! Version Fortran 95, Aug, 2006.
    ! Version of Aitp1d, without the call to Locate, which is called outside.
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    integer :: nx,i
    real(wp), dimension(nx) :: Vx,Fx
    real(wp) :: Xp,Fp,Aux
    
    !CALL Locate(Vx,nx,Xp,i)
    
    if ( i<=0 .or. i>=nx ) then
       FP= 0._wp
    else
       Aux= ( Xp-Vx(i) )/( Vx(i+1)-Vx(i) )
       Fp= Fx(i) + ( Fx(i+1)-Fx(i) )*Aux
    end if
    return
  end subroutine Aitp1d2
  
  subroutine Locate(Xx,N,X,J)
    ! Given an array Xx of lenght N, and given a value X, returns a value 
    ! J such that X is between Xx(J) and Xx(J+1).
    ! Xx must be monotonic, either increasing or decreasing.
    ! J=0 or J=N is returned to indicate that X is out of range.
    ! See NUMERICAL RECIPES.
    ! Version Fortran 95, Aug 2006.
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    integer :: N,J,JL,JU,JM
    real(wp), dimension (N) :: Xx
    real(wp) :: X
    
    JL= 0
    JU= N+1
10  if ( JU-JL .gt. 1 ) then
       JM= ( JU+JL ) / 2
       if ( (Xx(N).gt.Xx(1)).eqv.(X.gt.Xx(JM)) ) then
          JL= JM
       else
          JU= JM
       end if
       GO TO 10
    end if
    J= JL
    return
  end subroutine Locate

  subroutine output(WriteThis)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use common_params
    use common_arrays
    implicit none
    real(wp) :: Aux,Dqx,Dqz
    real(wp) :: Qx,Qz,Q
    integer :: i,k
    integer, parameter :: Step=1
    character(LEN=3) :: WriteThis
    
    select case(WriteThis)
    case("IL0")
       open(1,FILE='IL0.wt')
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),IL0(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),VQz(k),IL0(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) VQx(i),-VQz(nqz+1-k),IL0(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) VQx(i),VQz(k),IL0(i,k)
          end do
       end do
       close(1)
       open(1,FILE= "IL01D.wt")
       do k= 1,nqcd
          write(1,*) VQQ(k),IL01D(k)
       end do
       close(1)
       open(1,FILE='IL.wt')
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),IL(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),VQz(k),IL(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) VQx(i),-VQz(nqz+1-k),IL(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) VQx(i),VQz(k),IL(i,k)
          end do
       end do
       close(1)
       open(1,FILE= "IL1D.wt")
       do k= 1,nqcd
          write(1,*) VQQ(k),IL1D(k)
       end do
       close(1)
    case("GcL")
       open(1,FILE='GcollL.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,GcollLm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,GcollLp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,GcollLm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,GcollLp(i,k)! , &
                  ! GcollLp(i,k)*Geff
          end do
       end do
       close(1)
       open(1,FILE='GcollL1D.wt')
       do i= 1,nqcd
          write(1,*) VQQ(i),GcollL1D(i)
          ! ,GcollL1D(i)*Geff
       end do
       close(1)
    case("BrL")
       open(1,FILE='BremL.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,BremLm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,BremLp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,BremLm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,BremLp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='BremL1D.wt')
       do i= 1,nqcd
          write(1,*) VQQ(i),BremL1D(i)
       end do
       close(1)
       case("BGL")
       open(1,FILE='BrGcL.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,BrGcL(nqx+1-i,nqz+1-k)
                  ! BremLm(nqx+1-i,nqz+1-k)/GcollLm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,BrGcL(nqx+1-i,k)
             ! -BremLp(nqx+1-i,k)/GcollLp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,BrGcL(i,nqz+1-k)
             ! BremLm(i,nqz+1-k)/GcollLm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,BrGcL(i,k)
             ! -BremLp(i,k)/GcollLp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='BrGcL1D.wt')
       do i= 1,nqcd
          write(1,*) VQQ(i),BrGcL1D(i)
          ! -BremL1D(i)/GcollLp1D(i)
       end do
       close(1)
    end select
500 format(1x,5(e13.5e3,1x))
    return
  end subroutine output
end module sub_prog
   
program new_eff
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  use common_params
  use common_arrays
  use math_constants
  use phys_constants
  use sub_prog
  implicit none
  
  call allocate_arrays
  call init_vec
  call definitions
  call wave_init
  
end program new_eff
