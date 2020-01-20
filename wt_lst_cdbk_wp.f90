! Comments:
! 110205: Program to solve the time evolution of the equations for
!       electrons, Langmuir waves, ion-sound waves, and transverse waves, in
!       homogeneous medium, as formulated in notes "B11", dated Feb. 05, 2011.
!       Based on "wt_2d1d.f90", version 100826.
!       Uses subroutine "Split", with splitting method to solve the
!       equation for "Fe". 
!       The initial version does not yet include transverse waves, and should
!       be equivalent to "wt_2d1d.f90".
! 110219: Includes terms related to transverse waves, and the equation for
!       T waves.
! 110326: Corrects some mistakes found in notes B11:
!       1) L waves: in the scattering term related to L and T waves, eq. (43)
!       2) S waves: resonance conditions of decay involving L and T waves
!               (definitions of q1 ... q6).
!       3) T waves: in the scattering term related to L and T waves, eq. (50)
! 110403: Modifies the formulation according to notes E11.
! 110412: Modifies the formulation according to notes F11. These notes modify
!       the evaluation of the scattering terms, and correct a mistake in the
!       decay term involving L and S waves, in the equation for S waves. This
!       mistake was present in all previous notes.
! 110413: Modifies the dimensions of many arrays; instead of "-1:1" we use
!       dimension "2", creating vector "Iaux(2)" containing values "1" and "-1".
!       Warning: this version never worked properly. Although it compiles
!       correctly, it gives "segmentation errors", not identified, when it
!       is utilized. The version has been abandoned, and version 110425
!       was developed.
! 110425: Based on notes F11, as version 110412. The difference regarding that
!	version is that some quantities which were evaluated in subroutine
!	"Res_Cond" in version 110412 are evaluated in subroutines
!	"Coef_Lwave", "Coef_Swave" and "Coef_Twave". This is less demanding
!	in terms of allocated memory. Some tests have shown that the running
!	time of the code is not very different.
!	Introduced the option to include or not the scattering by electrons
!	(only spontaneous, see notes F-11).
!	This version corrects a mistake in the spontaneous terms for S waves
!	(which caused a small increase of the spectrum near the line Qx=0).
! 110630: Introduces a small modification in the use of boundary conditions.
!       The subroutine 'Boundary_Condition' is called after each of the
!       operators, Lx and Lz, in subroutine 'Split', instead of only at the
!       end, as in the previous versions. The subroutine is also modified
!       in the case of 'ConstantDF': When the value at the edge is negative,
!       keeps the previous value instead of attributing zero value.
!       The use of allocatable arrays is introduced, so that the dimensions
!       nux, nuz, nqx, nqz and nrz are input parameters of the code.
! 110723: Introduces the Runge-Kutta method for the wave equations, using
!       subroutine "RK4" (as in "wt_2d"). The present version considers
!       only the case of fixed time step. 
! 111115: Based on notex J11. Some mistakes found in previous versions of
! 	the notes have been corrected, and some features have been improved.
!	Particularly, some resonance conditions have been improved.
! 111211: Introduces some improvements regarding version 111115. One of them
!       is that the resonance values of "Qxp" or "Qzp" and of "Uxp" or "Uzp"
!       are evaluated only in subroutine "Res_Cond" and preserved for use
!       along the time evolution.
! 120302: Introduces a new approach for the scattering terms, following
!	notes C12 and "12b_comm", by Rudi Gaelzer. 
! 120325: Introduces an improvement in the new approach for the scattering 
!	terms, adding the contribution of induced scattering due to ions, 
!       following notes C12 dated March 11, 2012. 
! 120514: Introduces improvements on the resonance conditions for decay terms,
!       following notes D12 dated May 16, 2012.
! 120810: Further improvements on the resonance conditions for decay terms 
!       (particularly TdecayLL) and on scattering terms, following notes E12 
!       dated August 08, 2012.
! 120818: Introduces finite difference equations for the extreme points in
!       velocity space (except Ux=0, where the condition of zero derivative
!       is utilized). The subroutines for boundary conditions are not
!       necessary anymore, as well as the variable "BoundaryCondition" 
!       (see notes E12, dated August 18, 2012).
! 121229: Introduces improvements on the resonance conditions for decay terms,
!       following notes K12 dated Oct 22, 2012.
! 130121: Introduces improvements on the resonance conditions for decay terms,
!       following notes A13 dated Jan 21, 2013. Some resonance conditions
!       are squared, which can introduce spurious roots. A procedure to verify
!       which are the real roots is introduced. Preliminary results are not
!       encouraging, and the version has been abandoned, at least temporarily.
! 130126: Introduces different procedures for verification of the real roots.
!       Preliminary results are not very encouragin, and the version has been
!       abandoned, at least temporarily.
! 130128: The roots of the resonance conditions are found by a numerical 
!       procedure, which does not require squaring and generation of spurious
!       roots. Subroutine 'RTSAFE', from Numerical Recipes, is utilized to
!       locate the roots. 
! 130301: Some modifications have been made to the interpolation procedures
!       in subroutines "Funcx_RcdXXX" and "Funcz_RcdXXX". 
!	Another modification has been in the indexes of resonant values of
!	Qxdif and Qzdif, in the decay terms. The index "iresdif2" has been
!	introduced.
!       A mistake in the coefficients of the term related to decay involving
!       S and T waves ("CoefLSTdA" and "CoefLSTdB") in the subroutine
!       "Coef_Lwave": the coefficients are now multiplied by "AA", as they 
!       should be.
! 130306: Part of the coefficients of decay terms are evaluated in the
!	subroutine "Aux_Res_Cond" and stored as arrays "CxpLLSd", "CzpLLSd",
!	etc., instead of evaluated at each iteration.
!	Another modification is that the initial spectra are evaluated by
!	separated subroutines, "ILwave_Init", "ISwave_Init" and "ITwave_Init",
!	instead of a single subroutine as in previous versions.
!	This version was intended to be faster than version 130301. However,
!	the tests have shown that version 130301 is faster (about 75% of the
!	running time). Therefore, this version will not be used. We will
!	keep using version 130301.
! 130320: This version starts from version 130301; however, the roots of the
!       resonance conditions in the decay terms are now searched using
!       subroutine "RTBIS" (NUmerical Recipes), which uses the method of
!       bisection instead of the Newton-Raphson and bisection as in "RTSAFE".
!       The auxiliary routines ZBRAC2 and ZBRAC1 have been introduced. The
!       subroutines "Funcx_RcdXXX" and "Funcz_RcdXXX" have been transformed
!       into "functions". 
! 130322: Version based on version 130320. It eliminates all files related to
!       the decision based on "Qz>Qx", for decay coefficients. For decay terms
!       we look for the resonant value of "Qxp", using routines ZBRAC2 and
!       RTBIS. This reduces the amount of memory allocated to the resonance
!       conditions.
! 130920: Version based on version 130322. It gives options for the starting
!       condition of the T waves. In addition to the options existing before,
!       the new option assumes that the initial spectrum of T waves is given
!       by a condition of "turbulent equilibrium" (see notes F-12).
! 140404: The terms associated to spontaneous and induced emission for L waves
!       are evaluated using a grid with (nqx2,nqz2) points. The values for the
!       grid with (nqx,nqz) points, used for all other coefficients, are
!       obtained using interpolation. Experimentation has shown that the
!       scattering terms must be evaluated with large number of points, like
!       (71,71), while "Lemis" does not behave well for grids larger than
!       (51,51).
! 140423: The term associated to decay of L waves, involving LS waves, is also
!       evaluated using the (nqx2,nqz2) grid, and then added to the QL term and 
!       interpolated into the (nqx,nqz) grid. The results obtained are good,
!       but the code runs in a time nearly twice the time of the previous
!       version, probably due to repeated interpolation.
! 140424: All decay terms, as well as the QL terms, are evaluated using the
!       (nqx2,nqz2) grid and then interpolated into the (nqx,nqz) grid. The
!       locator indexes for the interpolation are evaluated as the code
!       starts, in subroutine "Res_Cond", stored in the arrays "IQresx2" and
!       "IQresz2", and used along the time evolution. Uses subroutine
!       "Aitp2d2", who utilizes the locator indexes obtained using "Locate".
!       The points of the (nqx,nqz) grid are interpolated into the (nqx2,nqz2)
!       grid, with locators stored into "IQresx", and "IQresz".
!       Another modification:
!       In order to increase the speed of the code, interpolations appearing
!       in the scattering terms are now made with "Aitp2d2", with locating
!       indexes evaluated in subroutine "Res_Cond" and stored to be used along
!       time evolution.
!       Correction of mistakes:
!       1) In the scattering terms, the quantity "Qs" was used before being 
!       defined. The mistake is present in all the most recent versions of the
!       code, in the segments corresponding to induced scattering by the ions. 
!       2) In the scattering for T waves, induced contribution of the ions,
!       the quantity Qstar was defined in previous versions without the SQRT.
!       3) In subroutine "Coef_Twave", the quantity "Phi" was used in the term
!       of scattering for T waves, but was not defined.
! 140505: Incorporates the correction of mistakes in the scattering terms, 
!       which has been made in version 140424. Also incorporates the evaluation
!       of locator indexes for the scattering terms and subsequent use of
!       subroutine "Aitp2d2", as in version 140424.
!       Quasilinear and decay terms evaluated using the (nqx,nqz) grid, also
!       used for the wave spectra appearing in the particle equation (as in
!       version 130920. Scattering terms, which need more resolution, are
!       evaluated using another grid (nqx2,nqz2), and their effect is then
!       interpolated for the points of grid (nqx,nqz).
! 140731: Quasilinear and decay terms for L and S waves continue to be
!       evaluated using the (nqx,nqz) grid. Scattering terms for all types of
!       waves and also the decay terms for T waves are evaluated using the grid
!       with (nqx2,nqz2), and their effect is then interpolated for the points
!       of the grid (nqx,nqz). 
! 140824: Quasilinear and decay terms for L, S, and T waves are evaluated using
!       the (nqx,nqz) grid.
!       Scattering terms for L waves are evaluated using the (nqx,nqz) grid.
!       Scattering terms for T waves are evaluated using the (nqx2,nqz2) grid,
!       and their effect is then interpolated for the points of grid (nqx,nqz).
!       Scattering terms for S waves are not taken into account.
!       This version is in may regards similar to version 140505, but it
!       incorporates some improvements adopted in version 140731, particularly
!       regarding the interpolation procedure near "Q=0".
!       Subroutine "Rebuild" is used only for "Fe", "IL", and "IS". The
!       spectrum of T waves, which can be very peaked near "Q=0" in the case
!       with a beam, is not smoothed out in this version.
!       Introduces parameter "nph", for the angular variable "Phip". It is
!       evaluates as the maximum between "nqz" and "nqz2".
! 140828: This version is basically the same as version 140731, with some of
!       the changes incorporated to version 140824:
!       Quasilinear and decay terms for L and S waves are evaluated using
!       the (nqx,nqz) grid.
!       Scattering terms for L waves are also evaluated using the (nqx,nqz) 
!!!       grid. Scattering terms for S waves are not taken into account.
!!!       Scattering and decay terms for T waves are evaluated using the
!!!       (nqx2,nqz2) grid, and their effect is then interpolated for the 
!!!       points of grid (nqx,nqz).
!!! 141015: This version introduces subroutine "Energy2", which evaluates the
!!!       energy density of waves L, S, and T, along the time evolution. It
!!!       also evaluates the energy associated to the fundamental and the
!!!       harmonics 2 and 3 of the T wave emission. These results appear in
!!!       file "Ewave.wt".
!!!       In subroutine "Output", it is now possible to generate file "ITQ1",
!!!       which contains the integral of the T wave intensity over angle dPhi,
!!!       and file "ITQ2", which contains the integral over "Q*dPhi".
!!! 150225: Corrects a sign mistake in components "xz", in "Coef_Coll" (this
!!!       subroutine introduces collisional effects, although it has not yet
!!!       been sufficiently tested).
!!!       For the collision coefficients, see page 7 of notes B15, dated 
!!!       Apr. 28, 2015 (equation (16), section 4). See also
!!!       page 27 of notes D08, dated Nov. 15, 2013.
!!!       Introduces character parameters "TimeEvol", "OneDfiles", "TwoDfiles",
!!!       "Onecolumnfiles", and "ADcoefficients", which allow to decide which
!!!       sets of output files to be generated, and if the time evolution will
!!!       be performed or if only files will be generated, from previous output
!!!       files. These parameters must be provided in file "Start.wt"
!!! 151126: This is a special version of program "wt_lst.f90", adapted to 
!!!       produce results for the collisional damping of L and S waves, for the 
!!!       case of a Maxwellian distribution, using the formalism presented in 
!!!       "Weak turbulence theory for collisional plasmas", by P. H. Yoon, L. F.
!!!       Ziebell, and E. P. Kontar, submitted for publication in September 2015.
!!!       Uses subroutine "Coll_Damping" to generate the collisional damping and
!!!       the Landau damping as well, for Maxwellian distribution
!!!       After printing the damping coefficients, there is a "STOP", to finish
!!!       the execution of the code, before the time evolution of the weak 
!!!       turbulence equations.
!!! 190306: Sabrina started to implement dynamic precision adding the following
!!!       line, which uses an intrinsic function, to the modules and subroutines:
!!!       use,intrinsic:: iso_fortran_env, only: wp=>real64
!!!       In order to change the precision, one has only to change the value
!!!       assingned to wp. The variable wp can have three values 32 (single
!!!       precision, 64 (double precision) and 128 (crazy precision).
!!!       Please, do not use d, or c before any intrinsic function, like "exp".
!!!       Do not use "dexp", for instance. The intrinsiic function will use the
!!!       precision assigned to the variable that is is being appied. Also, do
!!!       not use d to write double precision numbers either. Any real or complex
!!!       number should be written like this: 1.e-4_wp, 5._wp, 3.141592653589_wp.
!!!       Floating point variables should be declared as real(wp) or complex(wp).
!!!       Such practice assure that we can change the precision of the code
!!!       according to our needs, without having a headache.
!!!---------------------------------------------------------------------
!!! List of modules, subroutines, and functions:

!!! MODULE Common_Params
!!! MODULE Common_Arrays
!!! MODULE Allocate_Arrays
!!! MODULE Math_Constants
!!! MODULE Phys_Constants
!!! MODULE Sub_Prog

!!! SUBROUTINE Definitions
!!! SUBROUTINE Space_Profiles
!!! SUBROUTINE Init_Wave(Tau,Dqx,Dqz,Dux,Duz,Anorm,Epart,Ewave,EppEw,&
!!!     EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
!!! SUBROUTINE Iwave_Init(Qx,Qz,ILinit,ISinit,ITinit)
!!! SUBROUTINE Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
!!! SUBROUTINE Save_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
!!!     Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
!!! SUBROUTINE Read_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
!!!     Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
!!! SUBROUTINE Res_Cond
!!! SUBROUTINE Coef_A
!!! SUBROUTINE Coef_D
!!! SUBROUTINE Coef_Lwave(sigma,Dfdux,Dfduz,CoefA,CoefB)
!!! SUBROUTINE Coef_Swave(sigma,Dfdux,Dfduz,CoefA,CoefB)
!!! SUBROUTINE Coef_Twave(sigma,CoefA,CoefB)
!!! SUBROUTINE Aux_Coef_Decay_z(nqx,nqz,ires,iresdif,&
!!!  Qxp,Qzp,Qxdif,Qzdif,Wavep,Wavedif,VQx,VQz,Vauxqxp,Vauxqxdif,&
!!!  Iqp,Iqdif)
!!! SUBROUTINE Fnorm(Anorm)
!!! SUBROUTINE Energy(Epart,Ewave,EppEw)
!!! SUBROUTINE Energy2(EwaveL,EWaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
!!! SUBROUTINE Output_Coef(CoefChoice,CoefA,CoefB)
!!! SUBROUTINE Output(WriteChoice)
!!! SUBROUTINE Output2(WriteChoice)
!!! SUBROUTINE Split(Dux,Duz,DTau)
!!! SUBROUTINE Evol_Iwave(DTau,WaveType)
!!! SUBROUTINE Tridag(If,k,A,B,C,D,V,Beta,Gamma)
!!! SUBROUTINE Cor_Ampli(Vf,n)
!!! SUBROUTINE Cor_Ampli_2(Vf,Vfnew,n)
!!! SUBROUTINE Aitp1d2(nx,Vx,Fx,Xp,Fp,i)
!!! SUBROUTINE Locate(Xx,N,X,J)
!!! SUBROUTINE Aitp2d2(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp,i,j)
!!! SUBROUTINE Simpson(Vx,F,N,Res)
!!! SUBROUTINE Derivxy5p2d(nx,ny,Vx,Vy,Fxy,Dfdx,Dfdy)
!!! SUBROUTINE Derivxy5pln2d(nx,ny,Vx,Vy,Fxy,Dfdx,Dfdy)
!!! SUBROUTINE Derivx5p2d(nx,ny,Vx,Fxy,Dfdx)
!!! SUBROUTINE Derivy5p2d(nx,ny,Vy,Fxy,Dfdy)
!!! SUBROUTINE Coef_Coll
!!! SUBROUTINE Rebuild
!!! SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT)
!!! SUBROUTINE DERIVS(Tau,Ft,Dfdt,N)
!!! SUBROUTINE ZBRAC2(FUNC,X1,X2,SUCCES)

!!! REAL FUNCTION Funcx_RcdLLS(Qxp)
!!! REAL FUNCTION Funcx_RcdLLT(Qxp)
!!! REAL FUNCTION Funcx_RcdLST(Qxp)
!!! REAL FUNCTION Funcx_RcdLTT(Qxp)
!!! REAL FUNCTION Funcx_RcdSLL(Qxp)
!!! REAL FUNCTION Funcx_RcdSLT(Qxp)
!!! REAL FUNCTION Funcx_RcdTLL(Qxp)
!!! REAL FUNCTION Funcx_RcdTLS(Qxp)
!!! REAL FUNCTION Funcx_RcdTTL(Qxp)

!!! REAL FUNCTION ZL(Qx,Qz)
!!! REAL FUNCTION ZS(Qx,Qz)
!!! REAL FUNCTION ZT(Qx,Qz)
!!! REAL FUNCTION RTBIS(FUNC,X1,X2,XACC) 
!!! REAL FUNCTION UgL(Qx,Qz)
!!! REAL FUNCTION UgS(Qx,Qz)
!!! REAL FUNCTION PERF(X)
!!! REAL FUNCTION PERFC(X)
!!! REAL FUNCTION PERFCE(X)

!!!---------------------------------------------------------------------

!!!---------------------------------------------------------------------
!!! Modules with definitions:
!!!---------------------------------------------------------------------

module Common_Params
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  ! Read parameters:
  integer :: nqx,nqz,nux,nuz,nrz  ! nqx,nqz,nux,nuz,nrz must be odd numbers.
  integer :: nqx2,nqz2  ! nqx2,nqz2 must be odd numbers.
  integer :: AsympT
  real(wp) :: Kappae, Kappai
  real(wp) :: Ulim,Ucrit
  real(wp) :: Qxi,Qxf,Qzi,Qzf 
  real(wp) :: Iw0 
  real(wp) :: RatioNf,Uf,RTfTe 
  real(wp) :: RatioNb,Ub,RTbTe 
  real(wp) :: RTeTi,G
  real(wp) :: Ve2C2
  real(wp) :: ScaleDensInhom
  character(LEN=6) :: InitialLevel
  character(LEN=3) :: Lemis,LdecayLS,LdecayLT,LdecayST,LdecayTT,LscatLL,LscatLT
  character(LEN=3) :: Semis,SdecayLL,SdecayLT,Sscat
  character(LEN=3) :: TdecayLL,TdecayLS,TdecayTL,TscatLT
  character(LEN=3) :: CollTerm
  character(LEN=8) :: CollTermForm
  character(LEN=3) :: SpontEmis
  character(LEN=3) :: ScatElSpo
  character(LEN=3) :: Gcoll
  character(LEN=3) :: GcollEvol
  character(LEN=3) :: Bremss
  character(LEN=3) :: BremssEvol
  character(LEN=3) :: NewEff_Init
  !CHARACTER(LEN=3) :: ScatElInd
  character(LEN=3) :: RenormFe
  character(LEN=3) :: DerivLn

  ! Parameters generated in the code:
  real(wp) :: RMiMe
  real(wp) :: AA,Geff 
  real(wp) :: U0
  real(wp) :: AuxInitialLevel
  real(wp) :: AuxAlpha
  integer :: Auxinterp
  integer :: nph
  integer, parameter :: nqcd= 1024 !512 ! 1536  ! 2048 ! 1024  ! number of points in Q, for Coll_Damp !odd
  !JPi
  real(wp), parameter :: QxSqr=0.4_wp,QzSqr=0.12_wp
  !JPf

end module Common_Params

module Common_Arrays
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp), dimension(:), allocatable :: VRz
  real(wp), dimension(:), allocatable :: VQx
  real(wp), dimension(:), allocatable :: VQz
  real(wp), dimension(:), allocatable :: VQx2
  real(wp), dimension(:), allocatable :: VQz2
  real(wp), dimension(:), allocatable :: VQQ
  real(wp), dimension(:,:), allocatable :: ILp,ILm,ISp,ISm,ITp,ITm
  real(wp), dimension(:,:), allocatable :: GcollLp,GqlLp
  real(wp), dimension(:,:), allocatable :: GcollLm,GqlLm
  real(wp), dimension(:,:), allocatable :: GcollL 
  real(wp), dimension(:), allocatable:: GcollL1D,GqlL1D
  real(wp), dimension(:), allocatable:: GcollLp1D,GcollLm1D
  real(wp), dimension(:,:), allocatable:: GcollSp,GqlSp
  real(wp), dimension(:,:), allocatable:: GcollSm,GqlSm
  real(wp), dimension(:), allocatable:: GcollS1D,GqlS1D
  ! real(wp), dimension(:), allocatable:: GcollSp1D,GqlSp1D
  ! real(wp), dimension(:), allocatable:: GcollSm1D,GqlSm1D
  real(wp), dimension(:,:), allocatable :: BremL,BremS
  real(wp), dimension(:,:), allocatable :: BremLp,BremSp
  real(wp), dimension(:,:), allocatable :: BremLm,BremSm
  !real(wp), dimension(:,:), allocatable ::  BremLpiL0m,BremSpIS0m
  !real(wp), dimension(:,:), allocatable ::  BremLpIL0p,BremSpIS0p
  real(wp), dimension(:), allocatable :: BremL1D,BremS1D
  real(wp), dimension(:), allocatable :: VPhip
  real(wp), dimension(:), allocatable :: VUx
  real(wp), dimension(:), allocatable :: VUz
  real(wp), dimension(:,:), allocatable :: Fe0,Fe,Fi
  real(wp), dimension(:,:), allocatable :: Ax,Az,Dxx,Dxz,Dzx,Dzz
  real(wp), dimension(:,:), allocatable :: ColAx,ColAz,ColDxx,ColDxz,ColDzz
  real(wp), dimension(:,:), allocatable :: CoefARK,CoefBRK
  real(wp), dimension(:), allocatable :: Fepzm,Fepzp
  real(wp), dimension(:), allocatable :: Fepxm,Fepxp
  real(wp), dimension(:,:), allocatable :: Fepsm,Fepsp
  real(wp), dimension(:), allocatable :: VRNeNs,VRTeTs,VRTfTs,VRTbTs,VRTiTs,VDRNeNs
  ! Resonance conditions, and interpolation indexes:
  ! For the values of Qx,Qz in the second grid of the wave spectra:
  integer, dimension(:), allocatable :: IQresx2,IQresz2
  integer, dimension(:), allocatable :: IQresx,IQresz
  ! For the values of Qx,Qz in the wave spectra, scattering terms:
  integer, dimension(:,:,:), allocatable :: IQxLLL1,IQzLLL1
  integer, dimension(:,:,:,:,:,:), allocatable :: IQxLLL2,IQzLLL2
  integer, dimension(:,:,:), allocatable :: IQxLLT1,IQzLLT1
  integer, dimension(:,:,:,:,:,:), allocatable :: IQxLLT2,IQzLLT2
  integer, dimension(:,:,:), allocatable :: IQxTLT1,IQzTLT1
  integer, dimension(:,:,:,:,:,:), allocatable :: IQxTLT2,IQzTLT2
  ! For the coefficients Dij:
  real(wp), dimension(:,:,:,:), allocatable :: Qzr1D,Qzr2D
  real(wp), dimension(:,:,:,:), allocatable :: Qxr1D,Qxr2D
  integer, dimension(:,:,:,:), allocatable :: IQzr1D,IQzr2D
  integer, dimension(:,:,:,:), allocatable :: IQxr1D,IQxr2D
  ! L waves, spontaneous and induced emission:
  real(wp), dimension(:,:,:,:), allocatable :: UzrpLql,UzrmLql
  real(wp), dimension(:,:,:,:), allocatable :: UxrpLql,UxrmLql
  integer, dimension(:,:,:,:), allocatable :: IUzrpLql,IUzrmLql
  integer, dimension(:,:,:,:), allocatable :: IUxrpLql,IUxrmLql
  ! L waves, decay, L and S waves:
  real(wp), dimension(:,:,:,:,:,:,:), allocatable :: QxpLLSd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpLLSd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpdifLLSd
  ! L waves, decay, L and T waves:
  real(wp), dimension(:,:,:,:,:,:,:), allocatable :: QxpLLTd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpLLTd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpdifLLTd
  ! L waves, decay, S and T waves:
  real(wp), dimension(:,:,:,:,:,:,:), allocatable :: QxpLSTd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpLSTd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpdifLSTd
  ! L waves, decay, T and T waves:
  real(wp), dimension(:,:,:,:,:,:,:), allocatable :: QxpLTTd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpLTTd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpdifLTTd
  ! S waves, spontaneous and induced emission:
  real(wp), dimension(:,:,:,:), allocatable :: UzrpSql,UzrmSql
  real(wp), dimension(:,:,:,:), allocatable :: UxrpSql,UxrmSql
  integer, dimension(:,:,:,:), allocatable :: IUzrpSql,IUzrmSql
  integer, dimension(:,:,:,:), allocatable :: IUxrpSql,IUxrmSql
  ! S waves, decay, L and L waves:
  real(wp), dimension(:,:,:,:,:,:,:), allocatable :: QxpSLLd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpSLLd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpdifSLLd
  ! S waves, decay, L and T waves:
  real(wp), dimension(:,:,:,:,:,:,:), allocatable :: QxpSLTd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpSLTd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpdifSLTd
  ! T waves, decay, L and L waves:
  real(wp), dimension(:,:,:,:,:,:,:), allocatable :: QxpTLLd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpTLLd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpdifTLLd
  ! T waves, decay, L and S waves:
  real(wp), dimension(:,:,:,:,:,:,:), allocatable :: QxpTLSd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpTLSd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpdifTLSd
  ! T waves, decay, T and T waves:
  real(wp), dimension(:,:,:,:,:,:,:), allocatable :: QxpTTLd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpTTLd
  integer, dimension(:,:,:,:,:,:,:), allocatable :: IQxpdifTTLd
  ! Auxiliary for resonance conditions:
  real(wp), dimension(:), allocatable :: Aux1_Rcd
  integer, dimension(:), allocatable :: Aux2_Rcd
  ! Auxiliary for collisional damping:
  real(wp), dimension(:), allocatable :: Aux1_Gcoll
  integer, dimension(:), allocatable :: Aux2_Gcoll
  ! Auxiliary for electrostatic bremsstrahlung:
  real(wp), dimension(:), allocatable :: Aux1_Bremss
  integer, dimension(:), allocatable :: Aux2_Bremss
    
end module Common_Arrays

module Math_Constants
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none

  real(wp) :: Pi= 3.1415926535897932384626433832795_wp
  real(wp) :: Sqtwo= 1.41421356237309504880168872420969_wp
  real(wp) :: Infinity= 1.0E+30_wp
  real(wp) :: Degree= 0.01745329251994329576923690768488_wp
  real(wp) :: EpsMin= 1.0E-16_wp
  real(wp) :: Xacc= 1.0E-6_wp
  real(wp) :: Qmin= 1.0E-4_wp
  real(wp) :: Dqaux= 1.0E-2_wp

  complex(wp) :: Zi= (0._wp,1._wp)
  complex(wp) :: Zzero= (0._wp,0._wp)

end module Math_Constants

module Phys_Constants
  use,intrinsic:: iso_fortran_env, only: wp=>real64
  implicit none
  real(wp) :: Me= 9.1093856E-28_wp        ! electron mass (g)
  real(wp) :: Mi= 1.6726219E-24_wp        ! proton mass (g)
  real(wp) :: MeC2= 510.998946E+3_wp      ! electron mass (eV)
  real(wp) :: MpC2= 938.272081E+6_wp      ! proton mass (eV)
  real(wp) :: C_SI= 2.997925E+8_wp        ! speed of light (m/s)
  real(wp) :: C_cgs= 2.997925E+10_wp      ! speed of light (cm/s)

end module Phys_Constants

!---------------------------------------------------------------------
! Module with Functions and Subroutines:
!---------------------------------------------------------------------

module Sub_Prog
  use,intrinsic:: iso_fortran_env, only: wp=>real64
contains
  
  subroutine Allocate_Arrays
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    implicit none

    allocate (VRz(nrz))
    allocate (VQx(nqx))
    allocate (VQz(nqz))
    allocate (VQx2(nqx2))
    allocate (VQz2(nqz2))
    allocate (VQQ(nqcd))
    allocate (ILp(nqx,nqz),ILm(nqx,nqz))
    allocate (ISp(nqx,nqz),ISm(nqx,nqz))
    allocate (ITp(nqx,nqz),ITm(nqx,nqz))
    allocate (GcollLp(nqx,nqz),GqlLp(nqx,nqz))
    allocate (GcollLm(nqx,nqz),GqlLm(nqx,nqz))
    allocate (GcollL(nqx,nqz))
    allocate (GcollL1D(nqcd),GqlL1D(nqcd))
    allocate (GcollLp1D(nqcd),GcollLm1D(nqcd))
    allocate (GcollSp(nqx,nqz),GqlSp(nqx,nqz))
    allocate (GcollSm(nqx,nqz),GqlSm(nqx,nqz))
    allocate (GcollS1D(nqcd),GqlS1D(nqcd))
    allocate (BremL(nqx,nqz),BremS(nqx,nqz))
    allocate (BremL1D(nqcd),BremS1D(nqcd))
    allocate (BremLp(nqx,nqz),BremSp(nqx,nqz))
    allocate (BremLm(nqx,nqz),BremSm(nqx,nqz))
    !allocate (BremLpIL0m(nqx,nqz),BremSpIS0m(nqx,nqz))
    !allocate (BremLpIL0p(nqx,nqz),BremSpIS0p(nqx,nqz))
    allocate (VPhip(nph))
    allocate (VUx(nux))
    allocate (VUz(nuz))
    allocate (Fe0(nux,nuz),Fe(nux,nuz),Fi(nux,nuz))
    allocate (Ax(nux,nuz),Az(nux,nuz))
    allocate (Dxx(nux,nuz),Dxz(nux,nuz),Dzx(nux,nuz),Dzz(nux,nuz))
    allocate (ColAx(nux,nuz),ColAz(nux,nuz))
    allocate (CoefARK(nqx,nqz),CoefBRK(nqx,nqz))
    allocate (ColDxx(nux,nuz),ColDxz(nux,nuz),ColDzz(nux,nuz))
    allocate (Fepzm(nux),Fepzp(nux))
    allocate (Fepxm(nuz),Fepxp(nuz))
    allocate (Fepsm(nux,nuz),Fepsp(nux,nuz))
    allocate (VRNeNs(nrz),VRTeTs(nrz),VRTfTs(nrz),VRTbTs(nrz))
    allocate (VRTiTs(nrz),VDRNeNs(nrz))
    ! Resonance conditions:
    ! For the coefficients Dij:
    allocate (Qzr1D(nux,nuz,nqx,-1:1),Qzr2D(nux,nuz,nqx,-1:1))
    allocate (Qxr1D(nux,nuz,nqz,-1:1),Qxr2D(nux,nuz,nqz,-1:1))
    allocate (IQzr1D(nux,nuz,nqx,-1:1),IQzr2D(nux,nuz,nqx,-1:1))
    allocate (IQxr1D(nux,nuz,nqz,-1:1),IQxr2D(nux,nuz,nqz,-1:1))
    ! Auxiliary for resonance conditions:
    allocate (Aux1_Rcd(3))
    allocate (Aux2_Rcd(4))
    ! Auxiliary for collisional damping:
    allocate (Aux1_Gcoll(6))
    allocate (Aux2_Gcoll(4))
    ! Auxiliary for electrostatic bremsstrahlung:
    allocate(Aux1_bremss(6))
    allocate(Aux2_bremss(4))
  end subroutine Allocate_Arrays

  subroutine Definitions
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: RMpMe,RMeMp

    RMpMe= MpC2/MeC2
    RMeMp= MeC2/MpC2
    RMiMe= RMpMe
    AA= sqrt(1._wp+3._wp/RTeTi)/sqrt(RMiMe)/sqrt(2._wp)
    Geff= G/(2._wp*sqrt(2._wp)*(4._wp*Pi)**2)

  end subroutine Definitions

  subroutine Space_Profiles
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Rz,Drz,R
    real(wp) :: Rzf,Rzi
    integer :: m

    Rzf= 101._wp
    Rzi= 100._wp

    if (nrz==1) then
       VRz(1)= Rzi
    else
       Drz= (Rzf-Rzi)/(nrz-1)
       do m= 1,nrz
          VRz(m)= Rzi+(m-1)*Drz
       end do
       VRz(nrz)= Rzf
    end if

    do m= 1,nrz
       Rz= VRz(m)
       R= Rz

       VRNeNs(m)= 1._wp
       VDRNeNs(m)= 0._wp
       VRTeTs(m)= 1._wp   ! In the present version we assume homogeneous temperatures.
       VRTfTs(m)= RTfTe
       VRTbTs(m)= RTbTe
       VRTiTs(m)= 1._wp/RTeTi
    end do

    return
  end subroutine Space_Profiles


  !---------------------------------------------------------------------
  ! Functions and Subroutines:
  !---------------------------------------------------------------------

  subroutine Init_Wave(Tau,Dqx,Dqz,Dux,Duz,Anorm,Epart,Ewave,EppEw,&
       EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Fesum,Femax,Fef,Feb
    ! real(wp) :: Fekappa
    real(wp) :: ILinit,ISinit,ITinit
    real(wp), intent(out) :: Tau
    real(wp) :: Dqx,Dqz,Dux,Duz,Dqx2,Dqz2 
    real(wp) :: DPhip,Phip
    real(wp) :: Ux,Uz,U,U2,Qx,Qz 
    real(wp) :: Anorm,Epart,Ewave,EppEw
    real(wp) :: EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
    real(wp) :: BremL0,BremS0,GcollL0,GcollS0
    real(wp), dimension(nux,nuz) :: Dfdux,Dfduz
    integer :: i,j,k,nu2

    Tau= 0._wp      ! Initializes the time counter.
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
    Dqx2= (Qxf-0.1_wp*Qmin-Qxi-0.1_wp*Qmin)/(nqx2-1)
    do i= 1,nqx2
       VQx2(i)= Qxi+0.1_wp*Qmin+(i-1)*Dqx2
    end do
    Dqz2= (Qzf-0.1_wp*Qmin-Qzi-0.1_wp*Qmin)/(nqz2-1)
    do i= 1,nqz2
       VQz2(i)= Qzi+0.1_wp*Qmin+(i-1)*Dqz2
    end do

    ! Initializes an auxiliary array, to be used in the scattering terms:
    DPhip= (Pi/2._wp-0._wp)/(nph-1)
    do i= 1,nph
       Phip= 0._wp+(i-1)*DPhip
       VPhip(i)= Phip
    end do
    VPhip(nph)= Pi/2._wp

    call Coll_Damping
    call Bremsstrahlung
    
    ! Initializes the spectrum of L waves and S waves:
    select case(InitialLevel)

    case("Choose")
       ILp= Iw0
       ILm= Iw0
       ISp= Iw0
       ISm= Iw0
       ITp= Iw0
       ITm= Iw0
       AuxInitialLevel= 0.E0_wp

    case("Auto  ")
       AuxInitialLevel= 1.E0_wp
       do i= 1,nqx
          Qx= VQx(i)
          do j= 1,nqz
             Qz= VQz(j)
             GcollL0= GcollLp(i,j)!-GcollLm(i,j)
             GcollS0= GcollSp(i,j)
             BremL0= BremLp(i,j)
             BremS0= BremSp(i,j)
             call Iwave_Init(Qx,Qz,BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
             ILp(i,j)= ILinit
             ILm(i,j)= ILinit
             ISp(i,j)= ISinit
             ISm(i,j)= ISinit
             ITp(i,j)= ITinit
             ITm(i,j)= ITinit
          end do
       end do

    case("Null  ")
       Iw0= 0.E0_wp
       ILp= Iw0
       ILm= Iw0
       ISp= Iw0
       ISm= Iw0
       ITp= Iw0
       ITm= Iw0
       AuxInitialLevel= 0.E0_wp

    case DEFAULT
       open(98,FILE='Warning_Init_Wave.wt')
       write(98,*) ' InitialLevel= ',InitialLevel
       write(98,*) ' InitialLevel must be (Choose) or (Auto  ) or (Null  ) !!'
       close(98)
       stop

    end select

    Dux= (Ulim-0._wp)/(nux-1)
    do i= 1,nux
       VUx(i)= 0._wp+(i-1)*Dux
    end do
    Duz= (Ulim-(-Ulim))/(nuz-1)
    do i= 1,nuz
       VUz(i)= -Ulim+(i-1)*Duz
    end do

    !! Symmetrization of the array for Uz:
    if(mod(nuz,2) .ne. 0._wp) then
       nu2= (nuz-1)/2+1
       do i= 1,nu2
          VUz(nuz+1-i)= -VUz(i)
       end do
       VUz(nu2)= 0._wp
    end if

    ! ! Symmetrization of the array for Uz:
    ! nu2= (nuz-1)/2+1
    ! do i= 1,nu2
    !    VUz(nuz+1-i)= -VUz(i)
    ! end do
    ! VUz(nu2)= 0._wp

    ! Initializes the electron and ion distribution functions:
    do k= 1,nuz
       Uz= VUz(k)
       do i= 1,nux
          Ux= VUx(i)
          U2= Ux**2+Uz**2
          U= sqrt(U2)
          call Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
          Fe(i,k)= Fesum
          ! Fi(i,k)= (RMiMe*RTeTi)*Kappae_i/(Kappae_i-1._wp)/(Pi) &
          !      * (1._wp+RMiMe*RTeTi*U2/(Kappae_i-1._wp))**(-(Kappae_i+1._wp))
          Fi(i,k)= (RMiMe*RTeTi)*exp(-RMiMe*RTeTi*U2)/(Pi)
       end do
    end do
    ! write(13,*) RMiMe, RTeTi, RMiMe*RTeTi
    call Fnorm(Anorm)
    call Energy(Epart,Ewave,EppEw)
    call Energy2(EwaveL,EWaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
    Fe0= Fe

    ! Initializes the boundary conditions for the electron distribution: 
    select case(DerivLn)
    case("Yes")
       call Derivxy5pln2d(nux,nuz,VUx,VUz,Fe0,Dfdux,Dfduz)
    case("No ")
       call Derivxy5p2d(nux,nuz,VUx,VUz,Fe0,Dfdux,Dfduz)
    end select
    do k= 1,nuz
       Fepxm(k)= Dfdux(2,k)
       Fepxp(k)= Dfdux(nux-1,k)
    end do
    do i= 1,nux
       Fepzp(i)= Dfduz(i,nuz-1)
       Fepzm(i)= Dfduz(i,2)
    end do

    return
  end subroutine Init_Wave

  subroutine Iwave_Init(Qx,Qz,BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp) :: Qx,Qz,Q,Q2
    real(wp) :: BremL0,BremS0
    real(wp) :: GcollL0,GcollS0
    real(wp) :: Zlq,Zsq,Zsoq,Ztq,Aux
    real(wp) :: Zlq2,Zsq2,Muq
    real(wp) :: ILinit,ISinit,ITinit
    integer :: m,i,k

    m= 1   ! Auxiliary quantity, for the spatial profile.

    ! AA= SQRT(1._wp+3._wp/RTeTi)/SQRT(RMiMe)/SQRT(2._wp)
    ! Geff= G/(2._wp*SQRT(2._wp)*(4._wp*Pi)**2)

    if (AuxInitialLevel==0._wp) then
       ILinit= Iw0
       ISinit= Iw0
       ITinit= EpsMin
    else
       Q= sqrt(Qx**2+Qz**2)
       Q2= Q**2
       Zlq= ZL(Qx,Qz)
       Zlq2= Zlq**2
       Zsq= ZS(Qx,Qz)
       Zsq2= Zsq**2
       Ztq= ZT(Qx,Qz)
       Muq= Q2*Q*AA/2._wp
       Zsoq= AA*sqrt(VRTeTs(m))/sqrt(1._wp+Q**2/2._wp*VRTeTs(m)/VRNeNs(m))
       if(NewEff_Init=="Yes")then
!!!!!!!!!!!!!!!!!!!!!!!! Maxwellian !!!!!!!!!!!!!!!!!!!!!!!
          ILinit= (sqrt(Pi)*Geff*VRNeNs(m)**2/sqrt(VRTeTs(m)) &
            * exp(-Zlq2/Q2/VRTeTs(m))/Q2/Q + BremL0 )&
            / (2._wp*sqrt(Pi)*VRNeNS(m)/sqrt(VRTeTs(m)**3)*(Zlq2) &
            * exp(-Zlq2/Q2/VRTeTs(m))/Q2/Q - 2._wp*GcollL0)
          ISinit= (sqrt(Pi)*(AA/2._wp)*Geff &
               * sqrt(VRTeTs(m)/VRNeNs(m))*VRNeNs(m)*VRTeTs(m) &
               * (exp(-Zsq2/Q2/VRTeTs(m))/sqrt(VRTeTS(m)) &
               + exp(-Zsq2/Q2*RMiMe/VRTiTs(m)) &
               * sqrt(RMiMe/VRTiTs(m))) + BremS0) &
               / (2._wp*sqrt(Pi)*(AA/2._wp)*Zlq*Zsq &
               * sqrt(VRTeTs(m)/VRNeNs(m))*VRTeTs(m) &
               * (exp(-Zsq2/Q2/VRTeTs(m))/sqrt(VRTeTs(m)**3) &
               + exp(-Zsq2/Q2*RMiMe/VRTiTs(m)) &
               * sqrt(RMiMe/VRTiTs(m))/VRTiTs(m)) - 2._wp*GcollS0)
          ! ILinit= (SQRT(Pi)*Geff*VRNeNs(m)**2/SQRT(VRTeTs(m)*(Kappae-1.5_wp)) &
          !      * GAMMA(Kappae)/GAMMA(Kappae-0.5_wp) &
          !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(Kappae-1.5_wp))**(-Kappae)+Q*Q2*BremL0) &
          !      / (2._wp*SQRT(Pi)*VRNeNS(m)/(VRTeTs(m)*(Kappae-1.5_wp))**(1.5_wp)*(Zlq2) &
          !      * GAMMA(Kappae+1._wp)/GAMMA(Kappae-0.5_wp) &
          !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(Kappae-1.5_wp))**(-(Kappae+1._wp))-2._wp*Q*Q2*GcollL0)
          ! ISinit= (SQRT(pi)*Muq/Q/Q2*Geff &
          !      * SQRT(VRTeTs(m)/VRNeNs(m))*VRTeTs(m)*VRNeNs(m) &
          !      * (GAMMA(Kappae)/GAMMA(Kappae-0.5_wp)/SQRT((Kappae-1.5_wp)/VRTeTs(m)) &
          !      * (1._wp+Zsq2/Q2/(Kappae-1.5_wp)/VRTeTs(m))**(-Kappae) &
          !      + SQRT(RMiMe/VRTiTs(m)) &
          !      * EXP(-RMiMe/VRTiTs(m)*Zsq2/Q2))+BremS0) &
          !      / (SQRT(Pi)*2._wp*Muq/Q/Q2*Zsq*Zlq*SQRT(VRTeTs(m)/VRNeNs(m))*VRTeTs(m) &
          !      * (GAMMA(Kappae+1._wp)/GAMMA(Kappae-0.5_wp)/((Kappae-1.5_wp)*VRTeTs(m))**(1.5_wp) &
          !      * (1._wp+Zsq2/Q2/(Kappae-1.5_wp)/VRTeTs(m))**(-(Kappae+1._wp)) &
          !      + SQRT(RMiMe/VRTiTs(m))/VRTiTs(m) &
          !      * EXP(-RMiMe/VRTiTs(m)*Zsq2/Q2))-2._wp*GcollS0)

!!!!!!!!!!!!!!!!!!!!! alpha= 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! ILinit= (SQRT(Pi)*Geff*VRNeNs(m)**2/SQRT(VRTeTs(m)*(Kappae)) &
          !      * GAMMA(Kappae-1._wp)/GAMMA(Kappae-1.5_wp) &
          !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(Kappae))**(-(Kappae-1._wp))+Q*Q2*BremL0) &
          !      / (2._wp*SQRT(Pi)*VRNeNS(m)/SQRT(VRTeTs(m)**3*(Kappae)**3)*(Zlq2) &
          !      * GAMMA(Kappae)/GAMMA(Kappae-1.5_wp) &
          !      * (1._wp+Zlq2/Q2/VRTeTs(m)/(Kappae))**(-(Kappae))-2._wp*Q*Q2*GcollL0)

          ! ISinit= (SQRT(pi)*AA/2._wp*Geff &
          !      * SQRT(VRTeTs(m)/VRNeNs(m))*VRTeTs(m)*VRNeNs(m) &
          !      * (GAMMA(Kappae-1._wp)/GAMMA(Kappae-1.5_wp)/SQRT((Kappae)/VRTeTs(m)) &
          !      * (1._wp+Zsq2/Q2/(Kappae)/VRTeTs(m))**(-(Kappae-1._wp)) &
          !      + SQRT(RMiMe/VRTiTs(m)) &
          !      * EXP(-RMiMe/VRTiTs(m)*Zsq2/Q2))+BremS0) &
          !      / (SQRT(Pi)*AA*Zsq*Zlq*SQRT(VRTeTs(m)/VRNeNs(m))*VRTeTs(m) &
          !      * (GAMMA(Kappae)/GAMMA(Kappae-1.5_wp)/((Kappae)*VRTeTs(m))**(1.5_wp) &
          !      * (1._wp+Zsq2/Q2/(Kappae)/VRTeTs(m))**(-(Kappae)) &
          !      + SQRT(RMiMe/VRTiTs(m))/VRTiTs(m) &
          !      * EXP(-RMiMe/VRTiTs(m)*Zsq2/Q2))-2._wp*GcollS0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Ions and electrons with kappa distribution !!!

          ! ISinit=(SQRT(Pi)*(AA/2._wp)*Geff &
          !      * SQRT(VRTeTs(m)/VRNeNs(m))*VRNeNs(m)*VRTeTs(m) &
          !      * ((1._wp+Zsq2/Q2/VRTeTs(m)/(Kappae-1.5_wp))**(-Kappae)/SQRT(VRTeTS(m)*Kappae) &
          !      * GAMMA(Kappae)/GAMMA(Kappae-0.5_wp) &
          !      + (1._wp+Zsq2/Q2/(Kappai-1.5_wp)*RMiMe/VRTiTs(m))**(-(Kappai)) &
          !      * GAMMA(Kappai)/GAMMA(Kappai-0.5_wp) &
          !      * SQRT(RMiMe/VRTiTs(m)/(Kappai-1.5_wp))) + BremS0) &
          !      / (2._wp*SQRT(Pi)*(AA/2._wp)*Zlq*Zsq &
          !      * SQRT(VRTeTs(m)/VRNeNs(m))*VRTeTs(m) &
          !      * ((1._wp+Zsq2/Q2/VRTeTs(m)/(Kappae-1.5_wp))**(-(Kappae+1._wp)) &
          !      / SQRT(VRTeTs(m)**3*(Kappae-1.5_wp)**3) &
          !      * GAMMA(Kappae+1._wp)/GAMMA(Kappae-0.5_wp) &
          !      + (1._wp+Zsq2/Q2/(Kappai-1.5_wp)*RMiMe/VRTiTs(m))**(-(Kappai+1)) &
          !      * GAMMA(Kappai+1._wp)/GAMMA(Kappai-0.5_wp) &
          !      * SQRT(RMiMe/VRTiTs(m)/(Kappai-1.5_wp)**3)/VRTiTs(m)) - 2._wp*GcollS0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       else
          ILinit= Geff*VRNeNs(m)*VRTeTs(m)/2._wp/Zlq**2
          Aux= (exp(-(Zsoq**2)/VRTeTs(m))+sqrt(RMiMe*RTeTi) &
               * exp(-RMiMe/VRTiTs(m)*(Zsoq)**2)) &
               / (exp(-(Zsoq)**2/VRTeTs(m))+RTeTi*sqrt(RMiMe*RTeTi) &
               * exp(-RMiMe/VRTiTs(m)*(Zsoq)**2))
          ISinit= Geff*VRNeNs(m)*VRTeTs(m)/2._wp/Zlq/Zsq * Aux
       end if
       !   WRITE(13,400)Q,ILinit,ISinit
       !400 FORMAT(1x,5(e13.5e3,1x))
       if (AsympT==1) then
          ITinit= Geff*VRNeNs(m)*VRTeTs(m)/Ztq**2
       else
          !   ITinit= EpsMin
          ITinit= 0._wp
       end if
    end if
    return
  end subroutine Iwave_Init

  subroutine Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp) :: Ux,Uz
    real(wp) :: Fesum,Femax,Fef,Feb
    ! real(wp) :: Kappa,Fekappa
    integer :: m

    m= 1   ! Auxiliary quantity, for the spatial profile.

    U0= -(RatioNf*Uf+RatioNb*Ub)/(1._wp-RatioNf-RatioNb)

    !! ! Kappa distribution function for electrons !!!!!!!!!!!!!!!!!!!!!!!!!!
    !! ! Using alpha=1 !!!
    ! Fekappa= (1._wp-RatioNf-RatioNb)/(Pi*VRTeTs(m))*Kappae/(Kappae-1._wp) &
    !      * (1._wp+(Ux**2+(Uz-U0)**2)/(Kappae-1._wp)/VRTeTs(m))**(-(Kappae+1._wp))
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Femax= (1._wp-RatioNf-RatioNb)/(Pi*VRTeTs(m))*exp(-(Ux**2+(Uz-U0)**2)/VRTeTs(m))
    Fef= RatioNf/(Pi*VRTfTs(m))*exp(-(Ux**2+(Uz-Uf)**2)/VRTfTs(m)) 
    Feb= RatioNb/(Pi*VRTbTs(m))*exp(-(Ux**2+(Uz-Ub)**2)/VRTbTs(m))
    !Fesum= Fekappa+Fef+Feb
    Fesum= Femax+Fef+Feb

    ! Adding a kappa distribution to avoid values too small for large 'U':
    ! Kappa= 2._wp
    ! Fesum= Fesum+EpsMin*((1._wp+(Ux**2+Uz**2)/VRTeTs(m)/Kappa)**(-Kappa) &
    !      /(Pi*Kappa*VRTeTs(m))/SQRT(Pi))

    !! !! Kappa beam !!!!
    ! Fef= RatioNf/(Pi*VRTfTs(m))*Kappae/(Kappae-1._wp) &
    !      * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m)/(Kappae-1._wp))**(-(Kappae+1._wp))
    ! Feb= RatioNb/(Pi*VRTbTs(m))*Kappae/(Kappae-1._wp) &
    !      * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m)/(Kappae-1._wp))**(-(Kappae+1._wp))
    ! Fesum= Femax+Fef+Feb

    return
  end subroutine Fe_Init

  subroutine Save_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
       Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Tau 
    real(wp) :: Dqx,Dqz,Dux,Duz
    real(wp) :: Anorm0,Ewave0,Epart0,EppEw0
    real(wp) :: Rn,Rp,Rw,Rs
    real(wp) :: EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
    integer :: i,k,iflag

    Iflag= 1
    open(1,FILE='Out.wt')
    write(1,*) Iflag
    write(1,*) nux,nuz
    write(1,*) nqx,nqz
    write(1,*) nqx2,nqz2
    write(1,*) nrz
    write(1,*) Ulim
    write(1,*) Ucrit
    write(1,*) Qxi 
    write(1,*) Qxf 
    write(1,*) Qzi 
    write(1,*) Qzf 
    write(1,*) InitialLevel
    write(1,*) Iw0 
    write(1,*) AsympT
    write(1,*) Kappae
    write(1,*) Kappai
    write(1,*) RatioNf 
    write(1,*) Uf 
    write(1,*) RTfTe 
    write(1,*) RatioNb 
    write(1,*) Ub 
    write(1,*) RTbTe 
    write(1,*) RTeTi 
    write(1,*) G
    write(1,*) Ve2C2 
    write(1,*) Lemis
    write(1,*) LdecayLS
    write(1,*) LdecayLT
    write(1,*) LdecayST
    write(1,*) LdecayTT
    write(1,*) LscatLL
    write(1,*) LscatLT
    write(1,*) Semis
    write(1,*) SdecayLL
    write(1,*) SdecayLT
    write(1,*) Sscat
    write(1,*) TdecayLL
    write(1,*) TdecayLS
    write(1,*) TdecayTL
    write(1,*) TscatLT
    write(1,*) CollTerm
    write(1,*) CollTermForm
    write(1,*) Gcoll
    write(1,*) GcollEvol
    write(1,*) Bremss
    write(1,*) BremssEvol
    write(1,*) NewEff_Init
    write(1,*) SpontEmis
    write(1,*) ScatElSpo
    write(1,*) RenormFe
    write(1,*) DerivLn
    ! Parameters generated in the code:
    write(1,*) RMiMe
    write(1,*) AA,Geff 
    write(1,*) U0
    write(1,*) AuxInitialLevel
    write(1,*) Auxinterp

    do i= 1,nux
       write(1,*) VUx(i)
    end do
    do k= 1,nuz
       write(1,*) VUz(k)
    end do
    do i= 1,nux
       do k= 1,nuz
          write(1,*) Fe0(i,k)
       end do
    end do
    do i= 1,nux
       do k= 1,nuz
          write(1,*) Fe(i,k)
       end do
    end do
    do i= 1,nux
       do k= 1,nuz
          write(1,*) Fi(i,k)
       end do
    end do
    do i= 1,nqx
       write(1,*) VQx(i)
    end do
    do k= 1,nqz
       write(1,*) VQz(k)
    end do
    do i= 1,nqx2
       write(1,*) VQx2(i)
    end do
    do k= 1,nqz2
       write(1,*) VQz2(k)
    end do
    do i= 1,nph
       write(1,*) VPhip(i)
    end do
    do i= 1,nqx
       do k= 1,nqz
          write(1,*) ILp(i,k)
       end do
    end do
    do i= 1,nqx
       do k= 1,nqz
          write(1,*) ILm(i,k)
       end do
    end do
    do i= 1,nqx
       do k= 1,nqz
          write(1,*) ISp(i,k)
       end do
    end do
    do i= 1,nqx
       do k= 1,nqz
          write(1,*) ISm(i,k)
       end do
    end do
    do i= 1,nqx
       do k= 1,nqz
          write(1,*) ITp(i,k)
       end do
    end do
    do i= 1,nqx
       do k= 1,nqz
          write(1,*) ITm(i,k)
       end do
    end do
    do i= 1,nux
       write(1,*) Fepzp(i),Fepzm(i)
    end do
    do k= 1,nuz
       write(1,*) Fepxp(k),Fepxm(k)
    end do
    write(1,*) nux
    write(1,*) nuz
    write(1,*) nqx
    write(1,*) nqz
    write(1,*) nqx2
    write(1,*) nqz2
    write(1,*) nrz
    write(1,*) Dqx
    write(1,*) Dqz
    write(1,*) Dux
    write(1,*) Duz
    write(1,*) U0
    write(1,*) Anorm0
    write(1,*) Epart0
    write(1,*) Ewave0
    write(1,*) EppEw0
    write(1,*) Rn
    write(1,*) Rp
    write(1,*) Rw
    write(1,*) Rs
    write(1,*) EwaveL
    write(1,*) EwaveS
    write(1,*) EwaveT
    write(1,*) EwaveTF
    write(1,*) EwaveTH
    write(1,*) EwaveT3
    write(1,*) Tau
    close(1)

    return
  end subroutine Save_Results

  subroutine Read_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
       Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Tau
    real(wp) :: Dqx,Dqz,Dux,Duz
    real(wp) :: Anorm0,Ewave0,Epart0,EppEw0
    real(wp) :: Rn,Rp,Rw,Rs
    real(wp) :: EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
    integer :: nuxback,nuzback,nqxback,nqzback,nqx2back,nqz2back,nrzback
    integer :: i,k

    ! Parameters generated in the code:
    read(1,*) RMiMe
    read(1,*) AA,Geff 
    read(1,*) U0
    read(1,*) AuxInitialLevel
    read(1,*) Auxinterp
    !
    do i= 1,nux
       read(1,*) VUx(i)
    end do
    do k= 1,nuz
       read(1,*) VUz(k)
    end do
    do i= 1,nux
       do k= 1,nuz
          read(1,*) Fe0(i,k)
       end do
    end do
    do i= 1,nux
       do k= 1,nuz
          read(1,*) Fe(i,k)
       end do
    end do
    do i= 1,nux
       do k= 1,nuz
          read(1,*) Fi(i,k)
       end do
    end do
    do i= 1,nqx
       read(1,*) VQx(i)
    end do
    do k= 1,nqz
       read(1,*) VQz(k)
    end do
    do i= 1,nqx2
       read(1,*) VQx2(i)
    end do
    do k= 1,nqz2
       read(1,*) VQz2(k)
    end do
    do i= 1,nph
       read(1,*) VPhip(i)
    end do
    do i= 1,nqx
       do k= 1,nqz
          read(1,*) ILp(i,k)
       end do
    end do
    do i= 1,nqx
       do k= 1,nqz
          read(1,*) ILm(i,k)
       end do
    end do
    do i= 1,nqx
       do k= 1,nqz
          read(1,*) ISp(i,k)
       end do
    end do
    do i= 1,nqx
       do k= 1,nqz
          read(1,*) ISm(i,k)
       end do
    end do
    do i= 1,nqx
       do k= 1,nqz
          read(1,*) ITp(i,k)
       end do
    end do
    do i= 1,nqx
       do k= 1,nqz
          read(1,*) ITm(i,k)
       end do
    end do
    do i= 1,nux
       read(1,*) Fepzp(i),Fepzm(i)
    end do
    do k= 1,nuz
       read(1,*) Fepxp(k),Fepxm(k)
    end do
    read(1,*) nuxback
    read(1,*) nuzback
    read(1,*) nqxback
    read(1,*) nqzback
    read(1,*) nqx2back
    read(1,*) nqz2back
    read(1,*) nrzback
    read(1,*) Dqx
    read(1,*) Dqz
    read(1,*) Dux
    read(1,*) Duz
    read(1,*) U0
    read(1,*) Anorm0
    read(1,*) Epart0
    read(1,*) Ewave0
    read(1,*) EppEw0
    read(1,*) Rn
    read(1,*) Rp
    read(1,*) Rw
    read(1,*) Rs
    read(1,*) EwaveL
    read(1,*) EwaveS
    read(1,*) EwaveT
    read(1,*) EwaveTF
    read(1,*) EwaveTH
    read(1,*) EwaveT3
    read(1,*) Tau
    close(1)

    if(nuxback/=nux .or. nuzback/=nuz .or. nqxback/=nqx .or. nqzback/=nqz &
         .or. nrzback/=nrz .or. nqx2back/=nqx2 .or. nqz2back/=nqz2) then
       open(98,FILE='Warning_Read_Results.wt')
       write(98,*) ' nuxback= ',nuxback,' nux= ',nux
       write(98,*) ' nuzback= ',nuzback,' nuz= ',nuz
       write(98,*) ' nqxback= ',nqxback,' nqx= ',nqx
       write(98,*) ' nqzback= ',nqzback,' nqz= ',nqz
       write(98,*) ' nqx2back= ',nqx2back,' nqx2= ',nqx2
       write(98,*) ' nqz2back= ',nqz2back,' nqz2= ',nqz2
       write(98,*) ' nrzback= ',nrzback,' nrz= ',nrz
       write(98,*) ' The quantities nux, nuz, nqx, nqz, nqx2, nqz2 and nrz, can not be modified !'
       close(98)
       stop
    else
    end if

    return
  end subroutine Read_Results

  subroutine Res_Cond
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Ux,Uz,Qx,Qz
    real(wp) :: AbsUz,AbsUx
    real(wp) :: Q1,Q2,Qxp,Qzp
    real(wp) :: Qx2,Qz2,Q,Uresp,Uresm
    real(wp) :: Qxdif
    real(wp) :: Zlq,Zsq,Ztq
    real(wp) :: Muq
    real(wp) :: X1,X2
    real(wp) :: Qstar,Qs,Phi,Phip,Aux1,Aalpha
    integer :: i,l,m,sigma,S1,S2,Asp,Aspp
    integer :: iqx,kqz,kqzp
    integer :: ires
    integer :: iresq,kresq
    integer :: j
    integer :: SS
    logical :: Succes 

    j= 1   ! Auxiliary quantity, for the spatial profile.

    Qzr1D= 0._wp
    Qzr2D= 0._wp
    Qxr1D= 0._wp
    Qxr2D= 0._wp
    IQzr1D= 0._wp
    IQzr2D= 0._wp
    IQxr1D= 0._wp
    IQxr2D= 0._wp

    ! For the values of Qx,Qz in the second grid of the wave spectra:
    allocate (IQresx2(nqx2),IQresz2(nqz2))
    ! The points of the grid (nqx2,nqz2) are "inside" the grid (nqx,nqz), and
    ! therefore are correctly interpolated.
    do iqx= 1,nqx2
       Qx= VQx2(iqx)
       call Locate(VQx,nqx,Qx,iresq)
       IQresx2(iqx)= iresq
    end do
    do kqz= 1,nqz2
       Qz= VQz2(kqz)
       call Locate(VQz,nqz,Qz,kresq)
       IQresz2(kqz)= kresq
    end do
    ! For the values of Qx,Qz in the first grid of the wave spectra:
    allocate (IQresx(nqx),IQresz(nqz))
    ! The extreme points of the grid (nqx,nqz) are "outside" of the grid
    ! (nqx2,nqz2), and therefore need a correction for correct interpolation.
    ! This is made with routine "Aitp2d2b"
    do iqx= 1,nqx
       Qx= VQx(iqx)
       call Locate(VQx2,nqx2,Qx,iresq)
       IQresx(iqx)= iresq
    end do
    do kqz= 1,nqz
       Qz= VQz(kqz)
       call Locate(VQz2,nqz2,Qz,kresq)
       IQresz(kqz)= kresq
    end do

    ! For the values of Qx,Qz in the wave spectra, scattering terms:
    ! L waves, scattering involving L and L waves
    allocate (IQxLLL1(nqx,nqz,nph),IQzLLL1(nqx,nqz,nph))
    allocate (IQxLLL2(nqx,nqz,nph,-1:1,-1:1,-1:1),&
         IQzLLL2(nqx,nqz,nph,-1:1,-1:1,-1:1))
    do iqx= 1,nqx
       Qx= VQx(iqx)
       do kqz= 1,nqz
          Qz= VQz(kqz)
          Qx2= Qx**2
          Qz2= Qz**2
          Q2= Qx2+Qz2
          Q= sqrt(Q2)
          Phi= acos(Qz/Q)
          Zlq= ZL(Qx,Qz)
          Qstar= Q
          do i= 1,nph
             Phip= VPhip(i)
             Qxp= Qstar*sin(Phip)
             Qzp= Qstar*cos(Phip)
             call Locate(VQx,nqx,Qxp,ires)
             IQxLLL1(iqx,kqz,i)= ires
             call Locate(VQz,nqz,Qzp,ires)
             IQzLLL1(iqx,kqz,i)= ires
             do S2= -1,1,2
                do S1= -1,1,2
                   Aux1= (Q2+Qstar**2-2._wp*S1*Q*Qstar*cos(Phi-S1*S2*Phip))
                   if (Aux1<1.E-6_wp) then
                      Aalpha= RMiMe/VRTiTs(j)*(9._wp/4._wp)*(Qstar/Zlq)**2/1.E-6_wp 
                   else
                      Aalpha= RMiMe/VRTiTs(j)*(9._wp/4._wp)*(Qstar/Zlq)**2/Aux1 
                   end if
                   do SS= -1,1,2
                      Qs= Qstar+SS/sqrt(2._wp*Aalpha)
                      Qxp= Qs*sin(Phip)
                      Qzp= Qs*cos(Phip)
                      call Locate(VQx,nqx,Qxp,ires)
                      IQxLLL2(iqx,kqz,i,S2,S1,SS)= ires
                      call Locate(VQz,nqz,Qzp,ires)
                      IQzLLL2(iqx,kqz,i,S2,S1,SS)= ires
                   end do
                end do
             end do
          end do
       end do
    end do

    ! L waves, scattering involving L and T waves
    allocate (IQxLLT1(nqx,nqz,nph),IQzLLT1(nqx,nqz,nph))
    allocate (IQxLLT2(nqx,nqz,nph,-1:1,-1:1,-1:1),&
         IQzLLT2(nqx,nqz,nph,-1:1,-1:1,-1:1))
    do iqx= 1,nqx
       Qx= VQx(iqx)
       do kqz= 1,nqz
          Qz= VQz(kqz)
          Qx2= Qx**2
          Qz2= Qz**2
          Q2= Qx2+Qz2
          Q= sqrt(Q2)
          Phi= acos(Qz/Q)
          Zlq= ZL(Qx,Qz)
          Qstar= sqrt(3._wp/2._wp*Ve2C2)*Q
          do i= 1,nph
             Phip= VPhip(i)
             Qxp= Qstar*sin(Phip)
             Qzp= Qstar*cos(Phip)
             call Locate(VQx,nqx,Qxp,ires)
             IQxLLT1(iqx,kqz,i)= ires
             call Locate(VQz,nqz,Qzp,ires)
             IQzLLT1(iqx,kqz,i)= ires
             do S2= -1,1,2
                do S1= -1,1,2
                   Aux1= (Q2+Qstar**2-2._wp*S1*Q*Qstar*cos(Phi-S1*S2*Phip))
                   if (Aux1<1.E-6_wp) then
                      Aalpha= RMiMe/VRTiTs(j)*(1._wp/Ve2C2**2)*(Qstar/Zlq)**2/1.E-6_wp
                   else
                      Aalpha= RMiMe/VRTiTs(j)*(1._wp/Ve2C2**2)*(Qstar/Zlq)**2/Aux1 
                   end if
                   do SS= -1,1,2
                      Qs= Qstar+SS/sqrt(2._wp*Aalpha)
                      Qxp= Qs*sin(Phip)
                      Qzp= Qs*cos(Phip)
                      call Locate(VQx,nqx,Qxp,ires)
                      IQxLLT2(iqx,kqz,i,S2,S1,SS)= ires
                      call Locate(VQz,nqz,Qzp,ires)
                      IQzLLT2(iqx,kqz,i,S2,S1,SS)= ires
                   end do
                end do
             end do
          end do
       end do
    end do

    ! T waves, scattering involving L and T waves
    allocate (IQxTLT1(nqx2,nqz2,nph),IQzTLT1(nqx2,nqz2,nph))
    allocate (IQxTLT2(nqx2,nqz2,nph,-1:1,-1:1,-1:1),&
         IQzTLT2(nqx2,nqz2,nph,-1:1,-1:1,-1:1))
    do iqx= 1,nqx2
       Qx= VQx2(iqx)
       do kqz= 1,nqz2
          Qz= VQz2(kqz)
          Qx2= Qx**2
          Qz2= Qz**2
          Q2= Qx2+Qz2
          Q= sqrt(Q2)
          Phi= acos(Qz/Q)
          Zlq= ZL(Qx,Qz)
          Ztq= ZT(Qx,Qz)
          Qstar= sqrt(2._wp/3._wp/Ve2C2)*Q
          do i= 1,nph
             Phip= VPhip(i)
             Qxp= Qstar*sin(Phip)
             Qzp= Qstar*cos(Phip)
             call Locate(VQx,nqx,Qxp,ires)
             IQxTLT1(iqx,kqz,i)= ires
             call Locate(VQz,nqz,Qzp,ires)
             IQzTLT1(iqx,kqz,i)= ires
             do S2= -1,1,2
                do S1= -1,1,2
                   Aux1= (Q2+Qstar**2-2._wp*S1*Q*Qstar*cos(Phi-S1*S2*Phip))
                   if (Aux1<1.E-6_wp) then
                      Aalpha= RMiMe/VRTiTs(j)*(9._wp/4._wp)*(Qstar/Ztq)**2/1.E-6_wp
                   else
                      Aalpha= RMiMe/VRTiTs(j)*(9._wp/4._wp)*(Qstar/Ztq)**2/Aux1 
                   end if
                   do SS= -1,1,2
                      Qs= Qstar+SS/sqrt(2._wp*Aalpha)
                      Qxp= Qs*sin(Phip)
                      Qzp= Qs*cos(Phip)
                      call Locate(VQx,nqx,Qxp,ires)
                      IQxTLT2(iqx,kqz,i,S2,S1,SS)= ires
                      call Locate(VQz,nqz,Qzp,ires)
                      IQzTLT2(iqx,kqz,i,S2,S1,SS)= ires
                   end do
                end do
             end do
          end do
       end do
    end do

    ! For the coefficients Dij:
    do m= 1,nuz
       Uz= VUz(m)
       if (Uz==0._wp) Uz= sign(1._wp,Uz)*EpsMin
       AbsUz= abs(Uz)
       do l= 1,nux
          Ux= VUx(l)
          if (Ux==0._wp) Ux= EpsMin
          AbsUx= abs(Ux)
          if (AbsUz>=Ux) then
             do sigma= -1,1,2
                do i= 1,nqx
                   Qx= VQx(i)
                   Q1= (sigma*sqrt(VRNeNs(j))-Qx*Ux)/Uz
                   if (Q1>0._wp) then
                      call Locate(VQz,nqz,Q1,ires)
                      IQzr1D(l,m,i,sigma)= ires
                      Qzr1D(l,m,i,sigma)= Q1
                   else
                      IQzr1D(l,m,i,sigma)= nqz+1
                      Qzr1D(l,m,i,sigma)= 0._wp
                   end if
                   Q2= (sigma*sqrt(VRNeNs(j))+Qx*Ux)/Uz
                   if (Q2>0._wp) then
                      call Locate(VQz,nqz,Q2,ires)
                      IQzr2D(l,m,i,sigma)= ires
                      Qzr2D(l,m,i,sigma)= Q2
                   else
                      IQzr2D(l,m,i,sigma)= nqz+1
                      Qzr2D(l,m,i,sigma)= 0._wp
                   end if
                end do
             end do
          else
             do sigma= -1,1,2
                do i= 1,nqz
                   Qz= VQz(i)
                   Q1= (sigma*sqrt(VRNeNs(j))-Qz*Uz)/Ux
                   if (Q1>0._wp) then
                      call Locate(VQx,nqx,Q1,ires)
                      IQxr1D(l,m,i,sigma)= ires
                      Qxr1D(l,m,i,sigma)= Q1
                   else
                      IQxr1D(l,m,i,sigma)= nqx+1
                      Qxr1D(l,m,i,sigma)= 0._wp
                   end if
                   Q2= (-sigma*sqrt(VRNeNs(j))+Qz*Uz)/Ux
                   if (Q2>0._wp) then
                      call Locate(VQx,nqx,Q2,ires)
                      IQxr2D(l,m,i,sigma)= ires
                      Qxr2D(l,m,i,sigma)= Q2
                   else
                      IQxr2D(l,m,i,sigma)= nqx+1
                      Qxr2D(l,m,i,sigma)= 0._wp
                   end if
                end do
             end do
          end if
       end do
    end do

    if(Lemis=="Yes") then
       ! L waves, contribution due to spontaneous and induced emission:
       ! L waves, spontaneous and induced emission:
       allocate (UzrpLql(nqx,nqz,nux,-1:1),UzrmLql(nqx,nqz,nux,-1:1))
       allocate (UxrpLql(nqx,nqz,nuz,-1:1),UxrmLql(nqx,nqz,nuz,-1:1))
       allocate (IUzrpLql(nqx,nqz,nux,-1:1),IUzrmLql(nqx,nqz,nux,-1:1))
       allocate (IUxrpLql(nqx,nqz,nuz,-1:1),IUxrmLql(nqx,nqz,nuz,-1:1))
       UzrpLql= 0._wp
       UzrmLql= 0._wp
       UxrpLql= 0._wp
       UxrmLql= 0._wp
       IUzrpLql= 0
       IUzrmLql= 0
       IUxrpLql= 0
       IUxrmLql= 0
       do iqx= 1,nqx
          Qx= VQx(iqx)
          do kqz= 1,nqz
             Qz= VQz(kqz)
             Qx2= Qx**2
             Qz2= Qz**2
             Q2= Qx2+Qz2
             Q= sqrt(Q2)
             Zlq= ZL(Qx,abs(Qz))
             if (Qz>Qx) then
                do l= 1,nux
                   Ux= VUx(l)
                   do sigma= -1,1,2
                      Uresp= (sigma*Zlq+Qx*Ux)/Qz
                      Uresm= (sigma*Zlq-Qx*Ux)/Qz
                      UzrpLql(iqx,kqz,l,sigma)= Uresp
                      UzrmLql(iqx,kqz,l,sigma)= Uresm
                      call Locate(VUz,nuz,Uresp,ires)
                      IUzrpLql(iqx,kqz,l,sigma)= ires
                      call Locate(VUz,nuz,Uresm,ires)
                      IUzrmLql(iqx,kqz,l,sigma)= ires
                   end do
                end do
             else
                do l= 1,nuz
                   Uz= VUz(l)
                   do sigma= -1,1,2
                      Uresp= (sigma*Zlq-Qz*Uz)/Qx
                      Uresm= -(sigma*Zlq-Qz*Uz)/Qx
                      call Locate(VUx,nux,Uresp,ires)
                      IUxrpLql(iqx,kqz,l,sigma)= ires
                      UxrpLql(iqx,kqz,l,sigma)= Uresp
                      call Locate(VUx,nux,Uresm,ires)
                      IUxrmLql(iqx,kqz,l,sigma)= ires
                      UxrmLql(iqx,kqz,l,sigma)= Uresm
                   end do
                end do
             end if
          end do
       end do
    else
    end if

    if(LdecayLS=="Yes") then
       ! L waves, decay associated to L and S waves:
       ! L waves, decay, L and S waves:
       allocate (QxpLLSd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpLLSd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpdifLLSd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       QxpLLSd= 0._wp
       IQxpLLSd= 0
       IQxpdifLLSd= 0
       do iqx= 1,nqx
          Qx= VQx(iqx)
          do kqz= 1,nqz
             Qz= VQz(kqz)
             Qx2= Qx**2
             Qz2= Qz**2
             Q2= Qx2+Qz2
             do Asp= -1,1,2
                do Aspp= -1,1,2
                   do S1= -1,1,2
                      do S2= -1,1,2
                         Aux1_Rcd(1)= Qx
                         Aux1_Rcd(2)= Qz
                         Aux2_Rcd(1)= Asp
                         Aux2_Rcd(2)= Aspp
                         Aux2_Rcd(3)= S1
                         Aux2_Rcd(4)= S2
                         !      IF (Qz>Qx) THEN
                         !      IF (Qz<0._wp) THEN
                         !      ELSE
                         do kqzp= 1,nqz
                            Qzp= VQz(kqzp)
                            Aux1_Rcd(3)= Qzp
                            X1= 0._wp
                            X2= 1.0_wp
                            call ZBRAC2(Funcx_RcdLLS,X1,X2,Succes) 
                            if (Succes .eqv. .true.) then 
                               Qxp= RTBIS(Funcx_RcdLLS,X1,X2,Xacc)
                            else
                               Qxp= 0._wp
                            end if
                            if(Qxp>0._wp)then
                               Qxdif= Qx-S1*Qxp
                               call Locate(VQx,nqx,Qxp,ires)
                               IQxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                               QxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
                               call Locate(VQx,nqx,abs(Qxdif),ires)
                               IQxpdifLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                            else
                            end if
                         end do
                         !      END IF
                      end do
                   end do
                end do
             end do
          end do
       end do
    else
    end if

    if(LdecayLT=="Yes") then
       ! L waves, decay associated to L and T waves:
       ! L waves, decay, L and T waves:
       allocate (QxpLLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpLLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpdifLLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       QxpLLTd= 0._wp
       IQxpLLTd= 0
       IQxpdifLLTd= 0
       do iqx= 1,nqx
          Qx= VQx(iqx)
          do kqz= 1,nqz
             Qz= VQz(kqz)
             Qx2= Qx**2
             Qz2= Qz**2
             Q2= Qx2+Qz2
             do Asp= -1,1,2
                do Aspp= -1,1,2
                   do S1= -1,1,2
                      do S2= -1,1,2
                         Aux1_Rcd(1)= Qx
                         Aux1_Rcd(2)= Qz
                         Aux2_Rcd(1)= Asp
                         Aux2_Rcd(2)= Aspp
                         Aux2_Rcd(3)= S1
                         Aux2_Rcd(4)= S2
                         !      IF (Qz>Qx) THEN
                         !      IF (Qz<0._wp) THEN
                         !      ELSE
                         do kqzp= 1,nqz
                            Qzp= VQz(kqzp)
                            Aux1_Rcd(3)= Qzp
                            X1= 0._wp
                            X2= 1.0_wp
                            call ZBRAC2(Funcx_RcdLLT,X1,X2,Succes) 
                            if (Succes .eqv. .true.) then 
                               Qxp= RTBIS(Funcx_RcdLLT,X1,X2,Xacc)
                            else
                               Qxp= 0._wp
                            end if
                            if(Qxp>0._wp)then
                               Qxdif= Qx-S1*Qxp
                               call Locate(VQx,nqx,Qxp,ires)
                               IQxpLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                               QxpLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
                               call Locate(VQx,nqx,abs(Qxdif),ires)
                               IQxpdifLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                            else
                            end if
                         end do
                         !      END IF
                      end do
                   end do
                end do
             end do
          end do
       end do
    else
    end if

    if(LdecayST=="Yes") then
       ! L waves, decay associated to S and T waves:
       ! L waves, decay, S and T waves:
       allocate (QxpLSTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpLSTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpdifLSTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       QxpLSTd= 0._wp
       IQxpLSTd= 0
       IQxpdifLSTd= 0
       do iqx= 1,nqx
          Qx= VQx(iqx)
          do kqz= 1,nqz
             Qz= VQz(kqz)
             Qx2= Qx**2
             Qz2= Qz**2
             Q2= Qx2+Qz2
             do Asp= -1,1,2
                do Aspp= -1,1,2
                   do S1= -1,1,2
                      do S2= -1,1,2
                         Aux1_Rcd(1)= Qx
                         Aux1_Rcd(2)= Qz
                         Aux2_Rcd(1)= Asp
                         Aux2_Rcd(2)= Aspp
                         Aux2_Rcd(3)= S1
                         Aux2_Rcd(4)= S2
                         !      IF (Qz>Qx) THEN
                         !      IF (Qz<0._wp) THEN
                         !      ELSE
                         do kqzp= 1,nqz
                            Qzp= VQz(kqzp)
                            Aux1_Rcd(3)= Qzp
                            X1= 0._wp
                            X2= 1.0_wp
                            call ZBRAC2(Funcx_RcdLST,X1,X2,Succes) 
                            if (Succes .eqv. .true.) then 
                               Qxp= RTBIS(Funcx_RcdLST,X1,X2,Xacc)
                            else
                               Qxp= 0._wp
                            end if
                            if(Qxp>0._wp)then
                               Qxdif= Qx-S1*Qxp
                               call Locate(VQx,nqx,Qxp,ires)
                               IQxpLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                               QxpLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
                               call Locate(VQx,nqx,abs(Qxdif),ires)
                               IQxpdifLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                            else
                            end if
                         end do
                         !      END IF
                      end do
                   end do
                end do
             end do
          end do
       end do
    else
    end if

    if(LdecayTT=="Yes") then
       ! L waves, decay associated to T and T waves:
       ! L waves, decay, T and T waves:
       allocate (QxpLTTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpLTTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpdifLTTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       QxpLTTd= 0._wp
       IQxpLTTd= 0
       IQxpdifLTTd= 0
       do iqx= 1,nqx
          Qx= VQx(iqx)
          do kqz= 1,nqz
             Qz= VQz(kqz)
             Qx2= Qx**2
             Qz2= Qz**2
             Q2= Qx2+Qz2
             do Asp= -1,1,2
                do Aspp= -1,1,2
                   do S1= -1,1,2
                      do S2= -1,1,2
                         Aux1_Rcd(1)= Qx
                         Aux1_Rcd(2)= Qz
                         Aux2_Rcd(1)= Asp
                         Aux2_Rcd(2)= Aspp
                         Aux2_Rcd(3)= S1
                         Aux2_Rcd(4)= S2
                         !      IF (Qz>Qx) THEN
                         !      IF (Qz<0._wp) THEN
                         !      ELSE
                         do kqzp= 1,nqz
                            Qzp= VQz(kqzp)
                            Aux1_Rcd(3)= Qzp
                            X1= 0._wp
                            X2= 1.0_wp
                            call ZBRAC2(Funcx_RcdLTT,X1,X2,Succes) 
                            if (Succes .eqv. .true.) then 
                               Qxp= RTBIS(Funcx_RcdLTT,X1,X2,Xacc)
                            else
                               Qxp= 0._wp
                            end if
                            if(Qxp>0._wp)then
                               Qxdif= Qx-S1*Qxp
                               call Locate(VQx,nqx,Qxp,ires)
                               IQxpLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                               QxpLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
                               call Locate(VQx,nqx,abs(Qxdif),ires)
                               IQxpdifLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                            else
                            end if
                         end do
                         !      END IF
                      end do
                   end do
                end do
             end do
          end do
       end do
    else
    end if

    if(Semis=="Yes") then
       ! S waves, contribution due to spontaneous and induced emission:
       ! S waves, spontaneous and induced emission:
       allocate (UzrpSql(nqx,nqz,nux,-1:1),UzrmSql(nqx,nqz,nux,-1:1))
       allocate (UxrpSql(nqx,nqz,nuz,-1:1),UxrmSql(nqx,nqz,nuz,-1:1))
       allocate (IUzrpSql(nqx,nqz,nux,-1:1),IUzrmSql(nqx,nqz,nux,-1:1))
       allocate (IUxrpSql(nqx,nqz,nuz,-1:1),IUxrmSql(nqx,nqz,nuz,-1:1))
       UzrpSql= 0._wp
       UzrmSql= 0._wp
       UxrpSql= 0._wp
       UxrmSql= 0._wp
       IUzrpSql= 0
       IUzrmSql= 0
       IUxrpSql= 0
       IUxrmSql= 0
       do iqx= 1,nqx
          Qx= VQx(iqx)
          do kqz= 1,nqz
             Qz= VQz(kqz)
             Q2= Qx**2+Qz**2
             Q= sqrt(Q2)
             Muq= Q**3*AA/2._wp
             Zlq= ZL(Qx,Qz)
             Zsq= ZS(Qx,Qz)
             if (Qz>Qx) then 
                do l= 1,nux
                   Ux= VUx(l)
                   do sigma= -1,1,2
                      Uresp= (sigma*Zsq+Qx*Ux)/Qz
                      Uresm= (sigma*Zsq-Qx*Ux)/Qz
                      call Locate(VUz,nuz,Uresp,ires)
                      IUzrpSql(iqx,kqz,l,sigma)= ires
                      UzrpSql(iqx,kqz,l,sigma)= Uresp
                      call Locate(VUz,nuz,Uresm,ires)
                      IUzrmSql(iqx,kqz,l,sigma)= ires
                      UzrmSql(iqx,kqz,l,sigma)= Uresm
                   end do
                end do
             else
                do l= 1,nuz
                   Uz= VUz(l)
                   do sigma= -1,1,2
                      Uresp= (sigma*Zsq-Qz*Uz)/Qx
                      Uresm= -(sigma*Zsq-Qz*Uz)/Qx
                      call Locate(VUx,nux,Uresp,ires)
                      IUxrpSql(iqx,kqz,l,sigma)= ires
                      UxrpSql(iqx,kqz,l,sigma)= Uresp
                      call Locate(VUx,nux,Uresm,ires)
                      IUxrmSql(iqx,kqz,l,sigma)= ires
                      UxrmSql(iqx,kqz,l,sigma)= Uresm
                   end do
                end do
             end if
          end do
       end do
    else
    end if

    if(SdecayLL=="Yes") then
       ! S waves, decay associated to L and L waves:
       ! S waves, decay, L and L waves:
       allocate (QxpSLLd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpSLLd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpdifSLLd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       QxpSLLd= 0._wp
       IQxpSLLd= 0
       IQxpdifSLLd= 0
       do iqx= 1,nqx
          Qx= VQx(iqx)
          do kqz= 1,nqz
             Qz= VQz(kqz)
             Qx2= Qx**2
             Qz2= Qz**2
             Q2= Qx2+Qz2
             do Asp= -1,1,2
                do Aspp= -1,1,2
                   do S1= -1,1,2
                      do S2= -1,1,2
                         Aux1_Rcd(1)= Qx
                         Aux1_Rcd(2)= Qz
                         Aux2_Rcd(1)= Asp
                         Aux2_Rcd(2)= Aspp
                         Aux2_Rcd(3)= S1
                         Aux2_Rcd(4)= S2
                         !      IF (Qz>Qx) THEN
                         !      IF (Qz<0._wp) THEN
                         !      ELSE
                         do kqzp= 1,nqz
                            Qzp= VQz(kqzp)
                            Aux1_Rcd(3)= Qzp
                            X1= 0._wp
                            X2= 1.0_wp
                            call ZBRAC2(Funcx_RcdSLL,X1,X2,Succes) 
                            if (Succes .eqv. .true.) then 
                               Qxp= RTBIS(Funcx_RcdSLL,X1,X2,Xacc)
                            else
                               Qxp= 0._wp
                            end if
                            if(Qxp>0._wp)then
                               Qxdif= Qx-S1*Qxp
                               call Locate(VQx,nqx,Qxp,ires)
                               IQxpSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                               QxpSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
                               call Locate(VQx,nqx,abs(Qxdif),ires)
                               IQxpdifSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                            else
                            end if
                         end do
                         !      END IF
                      end do
                   end do
                end do
             end do
          end do
       end do
    else
    end if

    if(SdecayLT=="Yes") then
       ! S waves, decay associated to L and T waves:
       ! S waves, decay, L and T waves:
       allocate (QxpSLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpSLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpdifSLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
       QxpSLTd= 0._wp
       IQxpSLTd= 0
       IQxpdifSLTd= 0
       do iqx= 1,nqx
          Qx= VQx(iqx)
          do kqz= 1,nqz
             Qz= VQz(kqz)
             Qx2= Qx**2
             Qz2= Qz**2
             Q2= Qx2+Qz2
             do Asp= -1,1,2
                do Aspp= -1,1,2
                   do S1= -1,1,2
                      do S2= -1,1,2
                         Aux1_Rcd(1)= Qx
                         Aux1_Rcd(2)= Qz
                         Aux2_Rcd(1)= Asp
                         Aux2_Rcd(2)= Aspp
                         Aux2_Rcd(3)= S1
                         Aux2_Rcd(4)= S2
                         !      IF (Qz>Qx) THEN
                         !      IF (Qz<0._wp) THEN
                         !      ELSE
                         do kqzp= 1,nqz
                            Qzp= VQz(kqzp)
                            Aux1_Rcd(3)= Qzp
                            X1= 0._wp
                            X2= 1.0_wp
                            call ZBRAC2(Funcx_RcdSLT,X1,X2,Succes) 
                            if (Succes .eqv. .true.) then 
                               Qxp= RTBIS(Funcx_RcdSLT,X1,X2,Xacc)
                            else
                               Qxp= 0._wp
                            end if
                            if(Qxp>0._wp)then
                               Qxdif= Qx-S1*Qxp
                               call Locate(VQx,nqx,Qxp,ires)
                               IQxpSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                               QxpSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
                               call Locate(VQx,nqx,abs(Qxdif),ires)
                               IQxpdifSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                            else
                            end if
                         end do
                         !      END IF
                      end do
                   end do
                end do
             end do
          end do
       end do
    else
    end if

    if(TdecayLL=="Yes") then
       ! T waves, decay associated to L and L waves:
       ! T waves, decay, L and L waves:
       allocate (QxpTLLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpTLLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpdifTLLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
       QxpTLLd= 0._wp
       IQxpTLLd= 0
       IQxpdifTLLd= 0
       do iqx= 1,nqx2
          Qx= VQx2(iqx)
          do kqz= 1,nqz2
             Qz= VQz2(kqz)
             Qx2= Qx**2
             Qz2= Qz**2
             Q2= Qx2+Qz2
             do Asp= -1,1,2
                do Aspp= -1,1,2
                   do S1= -1,1,2
                      do S2= -1,1,2
                         Aux1_Rcd(1)= Qx
                         Aux1_Rcd(2)= Qz
                         Aux2_Rcd(1)= Asp
                         Aux2_Rcd(2)= Aspp
                         Aux2_Rcd(3)= S1
                         Aux2_Rcd(4)= S2
                         !      IF (Qz>Qx) THEN
                         !      IF (Qz<0._wp) THEN
                         !      ELSE
                         do kqzp= 1,nqz2
                            Qzp= VQz2(kqzp)
                            Aux1_Rcd(3)= Qzp
                            X1= 0._wp
                            X2= 1.0_wp
                            call ZBRAC2(Funcx_RcdTLL,X1,X2,Succes) 
                            if (Succes .eqv. .true.) then 
                               Qxp= RTBIS(Funcx_RcdTLL,X1,X2,Xacc)
                            else
                               Qxp= 0._wp
                            end if
                            if(Qxp>0._wp)then
                               Qxdif= Qx-S1*Qxp
                               call Locate(VQx2,nqx2,Qxp,ires)
                               IQxpTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                               QxpTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
                               call Locate(VQx2,nqx2,abs(Qxdif),ires)
                               IQxpdifTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                            else
                            end if
                         end do
                         !      END IF
                      end do
                   end do
                end do
             end do
          end do
       end do
    else
    end if

    if(TdecayLS=="Yes") then
       ! T waves, decay associated to L and S waves:
       ! T waves, decay, L and S waves:
       allocate (QxpTLSd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpTLSd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpdifTLSd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
       QxpTLSd= 0._wp
       IQxpTLSd= 0
       IQxpdifTLSd= 0
       do iqx= 1,nqx2
          Qx= VQx2(iqx)
          do kqz= 1,nqz2
             Qz= VQz2(kqz)
             Qx2= Qx**2
             Qz2= Qz**2
             Q2= Qx2+Qz2
             do Asp= -1,1,2
                do Aspp= -1,1,2
                   do S1= -1,1,2
                      do S2= -1,1,2
                         Aux1_Rcd(1)= Qx
                         Aux1_Rcd(2)= Qz
                         Aux2_Rcd(1)= Asp
                         Aux2_Rcd(2)= Aspp
                         Aux2_Rcd(3)= S1
                         Aux2_Rcd(4)= S2
                         !      IF (Qz>Qx) THEN
                         !      IF (Qz<0._wp) THEN
                         !      ELSE
                         do kqzp= 1,nqz2
                            Qzp= VQz2(kqzp)
                            Aux1_Rcd(3)= Qzp
                            X1= 0._wp
                            X2= 1.0_wp
                            call ZBRAC2(Funcx_RcdTLS,X1,X2,Succes) 
                            if (Succes .eqv. .true.) then 
                               Qxp= RTBIS(Funcx_RcdTLS,X1,X2,Xacc)
                            else
                               Qxp= 0._wp
                            end if
                            if(Qxp>0._wp)then
                               Qxdif= Qx-S1*Qxp
                               call Locate(VQx2,nqx2,Qxp,ires)
                               IQxpTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                               QxpTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
                               call Locate(VQx2,nqx2,abs(Qxdif),ires)
                               IQxpdifTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                            else
                            end if
                         end do
                         !      END IF
                      end do
                   end do
                end do
             end do
          end do
       end do
    else
    end if

    if(TdecayTL=="Yes") then
       ! T waves, decay associated to T and L waves:
       ! T waves, decay, T and L waves:
       allocate (QxpTTLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpTTLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
       allocate (IQxpdifTTLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
       QxpTTLd= 0._wp
       IQxpTTLd= 0
       IQxpdifTTLd= 0
       do iqx= 1,nqx2
          Qx= VQx2(iqx)
          do kqz= 1,nqz2
             Qz= VQz2(kqz)
             Qx2= Qx**2
             Qz2= Qz**2
             Q2= Qx2+Qz2
             do Asp= -1,1,2
                do Aspp= -1,1,2
                   do S1= -1,1,2
                      do S2= -1,1,2
                         Aux1_Rcd(1)= Qx
                         Aux1_Rcd(2)= Qz
                         Aux2_Rcd(1)= Asp
                         Aux2_Rcd(2)= Aspp
                         Aux2_Rcd(3)= S1
                         Aux2_Rcd(4)= S2
                         !      IF (Qz>Qx) THEN
                         !      IF (Qz<0._wp) THEN
                         !      ELSE
                         do kqzp= 1,nqz2
                            Qzp= VQz2(kqzp)
                            Aux1_Rcd(3)= Qzp
                            X1= 0._wp
                            X2= 1.0_wp
                            call ZBRAC2(Funcx_RcdTTL,X1,X2,Succes) 
                            if (Succes .eqv. .true.) then 
                               Qxp= RTBIS(Funcx_RcdTTL,X1,X2,Xacc)
                            else
                               Qxp= 0._wp
                            end if
                            if(Qxp>0._wp)then
                               Qxdif= Qx-S1*Qxp
                               call Locate(VQx2,nqx2,Qxp,ires)
                               IQxpTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                               QxpTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
                               call Locate(VQx2,nqx2,abs(Qxdif),ires)
                               IQxpdifTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
                            else
                            end if
                         end do
                         !      END IF
                      end do
                   end do
                end do
             end do
          end do
       end do
    else
    end if

    return
  end subroutine Res_Cond

  subroutine Coef_A
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Ux,Uz,Qx,Qx2,Qz,Qz2
    real(wp) :: AbsUz,AbsUx
    real(wp) :: Q1,Q12
    real(wp) :: Q2,Q22
    real(wp) :: Auxx,Auxz,AuxSE
    real(wp), dimension(nqx) :: VxAx,VxAz
    real(wp), dimension(nqz) :: VzAx,VzAz
    integer :: i,j,l,m,sigma

    !   Initialization of the Ai:
    Ax= 0._wp
    Az= 0._wp

    if(SpontEmis=="Yes") then
       AuxSE= 1.E0_wp
    else
       AuxSE= 0.E0_wp
    end if

    j= 1   ! Auxiliary quantity, for the spatial profile.
    do m= 1,nuz
       Uz= VUz(m)
       if (Uz==0._wp) Uz= sign(1._wp,Uz)*EpsMin
       AbsUz= abs(Uz)
       do l= 1,nux
          Ux= VUx(l)
          if (Ux==0._wp) Ux= EpsMin
          AbsUx= abs(Ux)
          if (AbsUz>=Ux) then
             do sigma= -1,1,2
                VxAx= 0._wp
                VxAz= 0._wp
                do i= 1,nqx
                   Qx= VQx(i)
                   Qx2= Qx**2
                   Q1= (sigma*sqrt(VRNeNs(j))-Qx*Ux)/Uz
                   if (Q1>0._wp) then
                      Q12= Q1**2
                      VxAx(i)= VxAx(i)+(Qx*Ux+Q1*Uz)/(Qx2+Q12) * Qx
                      VxAz(i)= VxAz(i)+(Qx*Ux+Q1*Uz)/(Qx2+Q12) * Q1
                   else
                      VxAx(i)= VxAx(i)+0._wp
                      VxAz(i)= VxAz(i)+0._wp
                   end if
                   Q2= (sigma*sqrt(VRNeNs(j))+Qx*Ux)/Uz
                   if (Q2>0._wp) then
                      Q22= Q2**2
                      VxAx(i)= VxAx(i)+(-Qx*Ux+Q2*Uz)/(Qx2+Q22) * (-Qx)
                      VxAz(i)= VxAz(i)+(-Qx*Ux+Q2*Uz)/(Qx2+Q22) * Q2
                   else
                      VxAx(i)= VxAx(i)+0._wp
                      VxAz(i)= VxAz(i)+0._wp
                   end if
                end do
                call Simpson(VQx,VxAx,nqx,Auxx)
                call Simpson(VQx,VxAz,nqx,Auxz)
                Ax(l,m)= Ax(l,m)+(2._wp*Geff/AbsUz)*Auxx * AuxSE
                Az(l,m)= Az(l,m)+(2._wp*Geff/AbsUz)*Auxz * AuxSE
             end do
          else
             do sigma= -1,1,2
                VzAx= 0._wp
                VzAz= 0._wp
                do i= 1,nqz
                   Qz= VQz(i)
                   Qz2= Qz**2
                   Q1= (sigma*sqrt(VRNeNs(j))-Qz*Uz)/Ux
                   if (Q1>0._wp) then
                      Q12= Q1**2
                      VzAx(i)= VzAx(i)+(Q1*Ux+Qz*Uz)/(Q12+Qz2) * Q1
                      VzAz(i)= VzAz(i)+(Q1*Ux+Qz*Uz)/(Q12+Qz2) * Qz
                   else
                      VzAx(i)= VzAx(i)+0._wp
                      VzAz(i)= VzAz(i)+0._wp
                   end if
                   Q2= (-sigma*sqrt(VRNeNs(j))+Qz*Uz)/Ux
                   if (Q2>0._wp) then
                      Q22= Q2**2
                      VzAx(i)= VzAx(i)+(-Q2*Ux+Qz*Uz)/(Q22+Qz2) * (-Q2)
                      VzAz(i)= VzAz(i)+(-Q2*Ux+Qz*Uz)/(Q22+Qz2) * Qz
                   else
                      VzAx(i)= VzAx(i)+0._wp
                      VzAz(i)= VzAz(i)+0._wp
                   end if
                end do
                call Simpson(VQz,VzAx,nqz,Auxx)
                call Simpson(VQz,VzAz,nqz,Auxz)
                Ax(l,m)= Ax(l,m)+(2._wp*Geff/AbsUx)*Auxx * AuxSE
                Az(l,m)= Az(l,m)+(2._wp*Geff/AbsUx)*Auxz * AuxSE
             end do
          end if
       end do
    end do

    return
  end subroutine Coef_A

  subroutine Coef_D !! restart from here
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Ux,Uz,Qx,Qx2,Qz,Qz2
    real(wp) :: AbsUz,AbsUx
    real(wp) :: EsLq1,EsLq2
    real(wp) :: Q1,Q12
    real(wp) :: Q2,Q22
    real(wp) :: Auxxx,Auxxz,Auxzz
    real(wp), dimension(nqx,nqz) :: Iwave
    real(wp), dimension(nqz) :: Vxq
    real(wp), dimension(nqx) :: Vzq
    real(wp), dimension(nqx) :: VxDxx,VxDxz,VxDzz
    real(wp), dimension(nqz) :: VzDxx,VzDxz,VzDzz
    integer :: i,l,m,sigma
    integer :: ires1,ires2
    integer :: j

    j= 1   ! Auxiliary quantity, for the spatial profile.
    !   Initialization of the Dij:
    Dxx= 0._wp
    Dxz= 0._wp
    Dzx= 0._wp
    Dzz= 0._wp

    do m= 1,nuz
       Uz= VUz(m)
       if (Uz==0._wp) Uz= sign(1._wp,Uz)*EpsMin
       AbsUz= abs(Uz)
       do l= 1,nux
          Ux= VUx(l)
          if (Ux==0._wp) Ux= EpsMin
          AbsUx= abs(Ux)
          if (AbsUz>=Ux) then
             do sigma= -1,1,2
                VxDxx= 0._wp
                VxDxz= 0._wp
                VxDzz= 0._wp
                if (sigma==1) then
                   Iwave= ILp
                else
                   Iwave= ILm
                end if
                do i= 1,nqx
                   Qx= VQx(i)
                   Qx2= Qx**2
                   Q1= Qzr1D(l,m,i,sigma)
                   if (Q1>0._wp) then
                      ires1= IQzr1D(l,m,i,sigma)
                      Q12= Q1**2
                      Vxq(:)= Iwave(i,:)
                      call Aitp1d2(nqz,VQz,Vxq,Q1,EsLq1,ires1)
                      VxDxx(i)= VxDxx(i)+Qx2/(Qx2+Q12) * EsLq1
                      VxDxz(i)= VxDxz(i)+Qx*Q1/(Qx2+Q12) * EsLq1
                      VxDzz(i)= VxDzz(i)+Q12/(Qx2+Q12) * EsLq1
                   else
                      VxDxx(i)= VxDxx(i)+0._wp
                      VxDxz(i)= VxDxz(i)+0._wp
                      VxDzz(i)= VxDzz(i)+0._wp
                   end if
                   Q2= Qzr2D(l,m,i,sigma)
                   if (Q2>0._wp) then
                      ires2= IQzr2D(l,m,i,sigma)
                      Q22= Q2**2
                      Vxq(:)= Iwave(i,:)
                      call Aitp1d2(nqz,VQz,Vxq,Q2,EsLq2,ires2)
                      VxDxx(i)= VxDxx(i)+Qx2/(Qx2+Q22) * EsLq2
                      VxDxz(i)= VxDxz(i)+(-Qx*Q2)/(Qx2+Q22) * EsLq2
                      VxDzz(i)= VxDzz(i)+Q22/(Qx2+Q22) * EsLq2
                   else
                      VxDxx(i)= VxDxx(i)+0._wp
                      VxDxz(i)= VxDxz(i)+0._wp
                      VxDzz(i)= VxDzz(i)+0._wp
                   end if
                end do
                call Simpson(VQx,VxDxx,nqx,Auxxx)
                call Simpson(VQx,VxDxz,nqx,Auxxz)
                call Simpson(VQx,VxDzz,nqx,Auxzz)
                Dxx(l,m)= Dxx(l,m)+(2._wp/AbsUz)*Auxxx
                Dxz(l,m)= Dxz(l,m)+(2._wp/AbsUz)*Auxxz
                Dzx(l,m)= Dzx(l,m)+(2._wp/AbsUz)*Auxxz
                Dzz(l,m)= Dzz(l,m)+(2._wp/AbsUz)*Auxzz
             end do
          else
             do sigma= -1,1,2
                VzDxx= 0._wp
                VzDxz= 0._wp
                VzDzz= 0._wp
                if (sigma==1) then
                   Iwave= ILp
                else
                   Iwave= ILm
                end if
                do i= 1,nqz
                   Qz= VQz(i)
                   Qz2= Qz**2
                   Q1= Qxr1D(l,m,i,sigma)
                   if (Q1>0._wp) then
                      ires1= IQxr1D(l,m,i,sigma)
                      Q12= Q1**2
                      Vzq(:)= Iwave(:,i)
                      call Aitp1d2(nqx,VQx,Vzq,Q1,EsLq1,ires1)
                      VzDxx(i)= VzDxx(i)+Q12/(Q12+Qz2) * EsLq1
                      VzDxz(i)= VzDxz(i)+Q1*Qz/(Q12+Qz2) * EsLq1
                      VzDzz(i)= VzDzz(i)+Qz2/(Q12+Qz2) * EsLq1
                   else
                      VzDxx(i)= VzDxx(i)+0._wp
                      VzDxz(i)= VzDxz(i)+0._wp
                      VzDzz(i)= VzDzz(i)+0._wp
                   end if
                   Q2= Qxr2D(l,m,i,sigma)
                   if (Q2>0._wp) then
                      ires2= IQxr2D(l,m,i,sigma)
                      Q22= Q2**2
                      Vzq(:)= Iwave(:,i)
                      call Aitp1d2(nqx,VQx,Vzq,Q2,EsLq2,ires2)
                      VzDxx(i)= VzDxx(i)+Q22/(Q22+Qz2) * EsLq2
                      VzDxz(i)= VzDxz(i)+(-Q2*Qz)/(Q22+Qz2) * EsLq2
                      VzDzz(i)= VzDzz(i)+Qz2/(Q22+Qz2) * EsLq2
                   else
                      VzDxx(i)= VzDxx(i)+0._wp
                      VzDxz(i)= VzDxz(i)+0._wp
                      VzDzz(i)= VzDzz(i)+0._wp
                   end if
                end do
                call Simpson(VQz,VzDxx,nqz,Auxxx)
                call Simpson(VQz,VzDxz,nqz,Auxxz)
                call Simpson(VQz,VzDzz,nqz,Auxzz)
                Dxx(l,m)= Dxx(l,m)+(2._wp/AbsUx)*Auxxx
                Dxz(l,m)= Dxz(l,m)+(2._wp/AbsUx)*Auxxz
                Dzx(l,m)= Dzx(l,m)+(2._wp/AbsUx)*Auxxz
                Dzz(l,m)= Dzz(l,m)+(2._wp/AbsUx)*Auxzz
             end do
          end if
       end do
    end do

    return
  end subroutine Coef_D

  subroutine Coef_Lwave(sigma,Dfdux,Dfduz,CoefA,CoefB)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Ux,Qx,Qz
    real(wp) :: Qx2,Qz2,Q2,Q,Uresp,Uresm
    real(wp) :: Zlq,Zlqp,Zlqdif
    real(wp) :: Ztqp,Ztqdif
    real(wp) :: Aux,AuxSE,Aux1,D1,AuxCoef
    real(wp) :: AuxScatElSpo
    !REAL(WP) :: AuxScatElInd
    real(wp) :: CoefEA,CoefLLSdA,CoefLLTdA,CoefLSTdA,CoefLTTdA,CoefLLLsA,CoefLLTsA
    real(wp) :: CoefEB,CoefLLSdB,CoefLLTdB,CoefLSTdB,CoefLTTdB,CoefLLLsB,CoefLLTsB
    real(wp) :: BremLA,GcollLB
    real(wp) :: Qxp,Qzp,Qxp2,Qzp2
    real(wp) :: Qxdif,Qzdif
    real(wp) :: Iqp,Iqdif,Iqpp,Iqpm
    real(wp) :: Fesum,Femax,Fef,Feb
    ! real(wp) :: Fekappa
    real(wp) :: Uz
    real(wp) :: Auxe,Aux0
    real(wp) :: ILinit,ISinit,ITinit
    real(wp) :: Phip,Phi,Sphip,Cphip
    real(wp) :: Aalpha,AuxNum
    real(wp) :: Qstar,Qs
    real(wp) :: Beta
    real(wp) :: Qp,Qp2
    real(wp) :: BremL0,BremS0,GcollL0,GcollS0
    real(wp), dimension(nux) :: VintuxA,VauxuxA
    real(wp), dimension(nux) :: VintuxB,VauxuxB
    real(wp), dimension(nuz) :: VauxuzA,VintuzA
    real(wp), dimension(nuz) :: VauxuzB,VintuzB
    real(wp), dimension(nux,nuz) :: Dfdux,Dfduz
    real(wp), dimension(nqz) :: VauxqzA,VauxqzB
    real(wp), dimension(nqx,nqz) :: CoefA,CoefB
    !REAL(WP), DIMENSION(nqx,nqz) :: BremssL
    real(wp), dimension(nqx,nqz,-1:1) :: EcalL
    real(wp), dimension(nqx,nqz,-1:1) :: EcalS
    real(wp), dimension(nqx,nqz,-1:1) :: EcalT
    real(wp), dimension(nqx) :: Vauxqxp,Vauxqxdif
    real(wp), dimension(nqx) :: VintqxA,VintqxB
    real(wp), dimension(nqz) :: VintqzA,VintqzB
    real(wp), dimension(nph) :: VintPhipA,VintPhipB
    integer :: i,k,l,sigma
    integer :: m
    integer :: sigmap,sigmapp
    integer :: iqx,kqz,kqzp
    integer :: iresp,iresm
    integer :: ires,iresdif,iresdif2,kres
    integer :: S1,S2,S3,SS
    integer :: Asp,Aspp

    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    if(SpontEmis=="Yes") then
       AuxSE= 1.E0_wp
    else
       AuxSE= 0.E0_wp
    end if
    if(ScatElSpo=="Yes") then
       AuxScatElSpo= 1.E0_wp
    else
       AuxScatElSpo= 0.E0_wp
    end if
    !IF(ScatElInd=="Yes") THEN
    ! AuxScatElInd= 1.E0_wp
    !ELSE
    ! AuxScatElInd= 0.E0_wp
    !END IF
    BremLA= 0._wp
    GcollLB= 0._wp
    CoefA= 0._wp
    CoefB= 0._wp
    do iqx= 1,nqx
       Qx= VQx(iqx)
       do kqz= 1,nqz
          Qz= VQz(kqz)
          Qx2= Qx**2
          Qz2= Qz**2
          Q2= Qx2+Qz2
          Q= sqrt(Q2)
          Zlq= ZL(Qx,Qz)
          Phi= acos(Qz/Q) 
          if(Lemis=="Yes") then
             ! Contribution due to spontaneous and induced emission:
             if (Qz>Qx) then
                do l= 1,nux
                   Ux= VUx(l)
                   Uresp= UzrpLql(iqx,kqz,l,sigma)
                   Uresm= UzrmLql(iqx,kqz,l,sigma)
                   iresp= IUzrpLql(iqx,kqz,l,sigma)
                   iresm= IUzrmLql(iqx,kqz,l,sigma)
                   if (iresp==0 .or. iresp==nuz) then
                      Uz= Uresp
                      call Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
                      VintuxA(l)= AuxSE * VRNeNs(m)*Geff*Fesum
                      ! VintuxB(l)= -(sigma*Zlq) &
                      !      * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) & 
                      !      - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                      !      * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m) &
                      !      / (Kappae-1._wp))**(-1._wp)*Fef &
                      !      - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                      !      * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m) &
                      !      / (Kappae-1._wp))**(-1._wp)*Feb &
                      !      + Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                      !      + Qz*2._wp*(Uz-Uf)*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                      !      * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m) &
                      !      / (Kappae-1._wp))**(-1._wp)*Fef &
                      !      + Qz*2._wp*(Uz-Ub)*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                      !      * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m) &
                      !      / (Kappae-1._wp))**(-1._wp)*Feb )
                      VintuxB(l)= -(sigma*Zlq) &
                           * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) &
                           - Qx*2._wp*Ux*Fef/VRTfTs(m) &
                           - Qx*2._wp*Ux*Feb/VRTbTs(m) &
                           + Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                           + Qz*2._wp*(Uz-Uf)*Fef/VRTfTs(m) &
                           + Qz*2._wp*(Uz-Ub)*Feb/VRTbTs(m) )
                   else
                      VauxuzA(:)= AuxSE * VRNeNs(m)*Geff*Fe(l,:)
                      VauxuzB(:)= (sigma*Zlq)*(-Qx*Dfdux(l,:) &
                           +Qz*Dfduz(l,:))
                      call Aitp1d2(nuz,VUz,VauxuzA,Uresp,Aux,iresp)
                      VintuxA(l)= Aux
                      call Aitp1d2(nuz,VUz,VauxuzB,Uresp,Aux,iresp)
                      VintuxB(l)= Aux
                   end if
                   if (iresm==0 .or. iresm==nuz) then
                      Uz= Uresm
                      call Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
                      VintuxA(l)= VintuxA(l) + AuxSE*VRNeNs(m)*Geff*Fesum
                      ! VintuxB(l)= VintuxB(l)+ (sigma*Zlq) &
                      !      * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) & 
                      !      - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                      !      * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m) &
                      !      / (Kappae-1._wp))**(-1._wp)*Fef &
                      !      - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                      !      * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m) &
                      !      / (Kappae-1._wp))**(-1._wp)*Feb &
                      !      - Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                      !      - Qz*2._wp*(Uz-Uf)*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                      !      * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m) &
                      !      / (Kappae-1._wp))**(-1._wp)*Fef &
                      !      - Qz*2._wp*(Uz-Ub)*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                      !      * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m) &
                      !      / (Kappae-1._wp))**(-1._wp)*Feb )
                      VintuxB(l)= VintuxB(l)+ (sigma*Zlq) &
                           * ( -Qx*2._wp*Ux*Femax/VRTeTs(m)-Qx*2. &
                           * Ux*Fef/VRTfTs(m) &
                           - Qx*2._wp*Ux*Feb/VRTbTs(m) &
                           - Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m)-Qz*2. &
                           * (Uz-Uf)*Fef/VRTfTs(m) &
                           - Qz*2._wp*(Uz-Ub)*Feb/VRTbTs(m) )
                   else
                      VauxuzA(:)= AuxSE * VRNeNs(m)*Geff*Fe(l,:)
                      VauxuzB(:)= (sigma*Zlq)*(Qx*Dfdux(l,:) &
                           +Qz*Dfduz(l,:))
                      call Aitp1d2(nuz,VUz,VauxuzA,Uresm,Aux,iresm)
                      VintuxA(l)= VintuxA(l) + Aux
                      call Aitp1d2(nuz,VUz,VauxuzB,Uresm,Aux,iresm)
                      VintuxB(l)= VintuxB(l) + Aux
                   end if
                end do
                call Simpson(VUx,VintuxA,nux,Aux)
                CoefEA= (Pi/Q2/(abs(Qz))) * VRNeNs(m)*Aux
                call Simpson(VUx,VintuxB,nux,Aux)
                CoefEB= (Pi/Q2/(abs(Qz))) * VRNeNs(m)*Aux
             else
                do l= 1,nuz
                   Uz= VUz(l)
                   Uresp= UxrpLql(iqx,kqz,l,sigma)
                   Uresm= UxrmLql(iqx,kqz,l,sigma)
                   iresp= IUxrpLql(iqx,kqz,l,sigma)
                   iresm= IUxrmLql(iqx,kqz,l,sigma)
                   if (iresm==0 .or. iresm==nux) then
                      if (iresm==nux) then
                         Ux= Uresm
                         call Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
                         VintuzA(l)= AuxSE * VRNeNs(m)*Geff*Fesum
                         ! VintuzB(l)= -(sigma*Zlq) &
                         !      * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) & 
                         !      - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                         !      * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m) &
                         !      / (Kappae-1._wp))**(-1._wp)*Fef &
                         !      - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                         !      * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m) &
                         !      / (Kappae-1._wp))**(-1._wp)*Feb &
                         !      + Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                         !      + Qz*2._wp*(Uz-Uf)*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                         !      * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m) &
                         !      / (Kappae-1._wp))**(-1._wp)*Fef &
                         !      + Qz*2._wp*(Uz-Ub)*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                         !      * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m) &
                         !      / (Kappae-1._wp))**(-1._wp)*Feb )
                         VintuzB(l)= -(sigma*Zlq) &
                              * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) &
                              - Qx*2._wp*Ux*Fef/VRTfTs(m) &
                              - Qx*2._wp*Ux*Feb/VRTbTs(m) &
                              + Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                              + Qz*2._wp*(Uz-Uf)*Fef/VRTfTs(m) &
                              + Qz*2._wp*(Uz-Ub)*Feb/VRTbTs(m) )
                      else
                         VintuzA(l)= 0._wp
                         VintuzB(l)= 0._wp
                      end if
                   else
                      VauxuxA(:)= AuxSE * VRNeNs(m)*Geff*Fe(:,l)
                      VauxuxB(:)= (sigma*Zlq)*(-Qx*Dfdux(:,l) &
                           +Qz*Dfduz(:,l))
                      call Aitp1d2(nux,VUx,VauxuxA,Uresm,Aux,iresm)
                      VintuzA(l)= Aux
                      call Aitp1d2(nux,VUx,VauxuxB,Uresm,Aux,iresm)
                      VintuzB(l)= Aux
                   end if
                   if (iresp==0 .or. iresp==nux) then
                      if (iresp==nux) then
                         Ux= Uresp
                         call Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
                         VintuzA(l)= VintuzA(l) + AuxSE * VRNeNs(m)*Geff*Fesum
                         ! VintuzB(l)= VintuzB(l)+ (sigma*Zlq) &
                         !      * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) & 
                         !      - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                         !      * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m) &
                         !      / (Kappae-1._wp))**(-1._wp)*Fef &
                         !      - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                         !      * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m) &
                         !      / (Kappae-1._wp))**(-1._wp)*Feb &
                         !      - Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                         !      - Qz*2._wp*(Uz-Uf)*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                         !      * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m) &
                         !      / (Kappae-1._wp))**(-1._wp)*Fef &
                         !      - Qz*2._wp*(Uz-Ub)*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                         !      * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m) &
                         !      / (Kappae-1._wp))**(-1._wp)*Feb )
                         VintuzB(l)= VintuzB(l) + (sigma*Zlq) &
                              * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) &
                              - Qx*2._wp*Ux*Fef/VRTfTs(m) &
                              - Qx*2._wp*Ux*Feb/VRTbTs(m) &
                              - Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                              - Qz*2._wp*(Uz-Uf)*Fef/VRTfTs(m) &
                              - Qz*2._wp*(Uz-Ub)*Feb/VRTbTs(m) )
                      else
                         VintuzA(l)= VintuzA(l)+ 0._wp
                         VintuzB(l)= VintuzB(l)+ 0._wp
                      end if
                   else
                      VauxuxA(:)= AuxSE * VRNeNs(m)*Geff*Fe(:,l)
                      VauxuxB(:)= (sigma*Zlq)*(Qx*Dfdux(:,l) &
                           +Qz*Dfduz(:,l))
                      call Aitp1d2(nux,VUx,VauxuxA,Uresp,Aux,iresp)
                      VintuzA(l)= VintuzA(l) + Aux
                      call Aitp1d2(nux,VUx,VauxuxB,Uresp,Aux,iresp)
                      VintuzB(l)= VintuzB(l) + Aux
                   end if
                end do
                call Simpson(VUz,VintuzA,nuz,Aux)
                CoefEA= (Pi/Q2/(abs(Qx))) * VRNeNs(m)*Aux
                call Simpson(VUz,VintuzB,nuz,Aux)
                CoefEB= (Pi/Q2/(abs(Qx))) * VRNeNs(m)*Aux
             end if
          else
             CoefEA= 0._wp
             CoefEB= 0._wp
          end if

          if(LdecayLS=="Yes") then
             ! Contribution due to spontaneous and induced decay involving L and S waves:
             EcalL(:,:,-1)= ILm(:,:)
             EcalL(:,:,+1)= ILp(:,:)
             EcalS(:,:,-1)= ISm(:,:)
             EcalS(:,:,+1)= ISp(:,:)
             !   IF (Qz>Qx) THEN
             !   IF (Qz<0._wp) THEN
             !   ELSE
             VauxqzA= 0._wp   ! Initializes the integrand of the Qzp integrals
             VauxqzB= 0._wp
             do kqzp= 1,nqz
                Qzp= VQz(kqzp)
                Qzp2= Qzp**2
                do sigmap= -1,1,2
                   Asp= sigmap/sigma
                   Vauxqxp(:)= EcalL(:,kqzp,sigmap)
                   do sigmapp= -1,1,2
                      Aspp= sigmapp/sigma
                      do S1= -1,1,2
                         do S2= -1,1,2
                            Qxp= QxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            ires= IQxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            iresdif= IQxpdifLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            Qxp2= Qxp**2
                            Qxdif= Qx-S1*Qxp
                            Qzdif= Qz-S2*Qzp
                            S3= sign(1._wp,Qzdif)
                            Zlqp= ZL(Qxp,Qzp)
                            Zlqdif= ZL(Qxdif,abs(Qzdif))
                            D1= -S2*sigmap*(1.5_wp*VRTeTs(m)/sqrt(VRNeNs(m))*Qxp &
                                 /sqrt(1._wp+1.5_wp*Beta*(Qxp2+Qzp2))) &
                                 + S1*S3*sigmapp*AA*sqrt(VRTeTs(m))*Qxdif &
                                 /sqrt(Qxdif**2+Qzdif**2+EpsMin) &
                                 /(1._wp+Beta*(Qxdif**2+Qzdif**2)/2._wp)**(1.5_wp)
                            AuxCoef= sqrt((Qxdif)**2+(Qzdif)**2) &
                                 * (S1*Qx*Qxp+S2*Qz*Qzp)**2 &
                                 /(Qxp2+Qzp2+EpsMin)/(abs(D1)+EpsMin)
                            iresdif2= kqz-S2*kqzp
                            Vauxqxdif(:)= EcalS(:,abs(iresdif2+S2 &
                                 * sign(1,iresdif2)),sigmapp)
                            call Aux_Coef_Decay_z(nqx,nqz,ires,iresdif,&
                                 Qxp,Qzp,Qxdif,Qzdif,"L","S",VQx,VQz, &
                                 Vauxqxp,Vauxqxdif,Iqp,Iqdif,nqcd,VQQ, &
                                 BremL1D,BremS1D,GcollL1D,GcollS1D)
                            VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma &
                                 * Zlq*Iqp*Iqdif
                            VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef &
                                 * (S2*sigmap*Zlqp*Iqdif &
                                 +S3*sigmapp*Zlqdif*Iqp)
                         end do
                      end do
                   end do
                end do
             end do
             call Simpson(VQz,VauxqzA,nqz,Aux)
             CoefLLSdA= Aux
             call Simpson(VQz,VauxqzB,nqz,Aux)
             CoefLLSdB= Aux
             !   END IF
             CoefLLSdA= AA/VRNeNs(m)/sqrt(VRNeNs(m)*VRTeTs(m)) &
                  * (sigma*Zlq/Q2) * CoefLLSdA
             CoefLLSdB= AA/VRNeNs(m)/sqrt(VRNeNs(m)*VRTeTs(m)) &
                  * (sigma*Zlq/Q2) * CoefLLSdB
          else
             CoefLLSdA= 0._wp
             CoefLLSdB= 0._wp
          end if

          if(LdecayLT=="Yes") then
             ! Contribution due to spontaneous and induced decay involving L and T waves:
             EcalL(:,:,-1)= ILm(:,:)
             EcalL(:,:,+1)= ILp(:,:)
             EcalT(:,:,-1)= ITm(:,:)
             EcalT(:,:,+1)= ITp(:,:)
             !   IF (Qz>Qx) THEN
             !   IF (Qz<0._wp) THEN
             !   ELSE
             VauxqzA= 0._wp   ! Initializes the integrand of the Qzp integrals
             VauxqzB= 0._wp
             do kqzp= 1,nqz
                Qzp= VQz(kqzp)
                Qzp2= Qzp**2
                do sigmap= -1,1,2
                   Asp= sigmap/sigma
                   Vauxqxp(:)= EcalL(:,kqzp,sigmap)
                   do sigmapp= -1,1,2
                      Aspp= sigmapp/sigma
                      do S1= -1,1,2
                         do S2= -1,1,2
                            Qxp= QxpLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            ires= IQxpLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            iresdif= IQxpdifLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            Qxp2= Qxp**2
                            Qxdif= Qx-S1*Qxp
                            Qzdif= Qz-S2*Qzp
                            S3= sign(1._wp,Qzdif)
                            Zlqp= ZL(Qxp,Qzp)
                            Ztqdif= ZT(Qxdif,abs(Qzdif))
                            D1= -S2*sigmap*(1.5_wp*Qxp/sqrt(1._wp+1.5_wp*(Qxp2+Qzp2)) &
                                 +S1*S3*sigmapp*Qxdif/Ve2C2 &
                                 / sqrt(1._wp+(Qxdif**2+Qzdif**2)/Ve2C2))
                            AuxCoef= (Q2*(Qxp2+Qzp2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
                                 /(Qxp2+Qzp2+EpsMin)/(Qxdif**2+Qzdif**2+EpsMin) &
                                 *((Qxp2+Qzp2)/S2/sigmap/Zlqp &
                                 + Q2/sigma/Zlq)**2/(abs(D1)+EpsMin)
                            iresdif2= kqz-S2*kqzp
                            Vauxqxdif(:)= EcalT(:,abs(iresdif2 &
                                 + S2*sign(1,iresdif2)),sigmapp)
                            call Aux_Coef_Decay_z(nqx,nqz,ires,iresdif,&
                                 Qxp,Qzp,Qxdif,Qzdif,"L","T",VQx,VQz,Vauxqxp, &
                                 Vauxqxdif,Iqp,Iqdif,nqcd,VQQ, &
                                 BremL1D,BremS1D,GcollL1D,GcollS1D)
                            VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma &
                                 * Zlq*Iqp*Iqdif
                            VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef &
                                 * (S2*sigmap*Zlqp*Iqdif &
                                 +S3*sigmapp*2._wp*Ztqdif*Iqp)
                         end do
                      end do
                   end do
                end do
             end do
             call Simpson(VQz,VauxqzA,nqz,Aux)
             CoefLLTdA= Aux
             call Simpson(VQz,VauxqzB,nqz,Aux)
             CoefLLTdB= Aux
             !   END IF
             CoefLLTdA= (1._wp/16._wp) * (sigma*Zlq/Q2) * CoefLLTdA
             CoefLLTdB= (1._wp/16._wp) * (sigma*Zlq/Q2) * CoefLLTdB
          else
             CoefLLTdA= 0._wp
             CoefLLTdB= 0._wp
          end if

          if(LdecayST=="Yes") then
             ! Contribution due to spontaneous and induced decay involving S and T waves:
             EcalS(:,:,-1)= ISm(:,:)
             EcalS(:,:,+1)= ISp(:,:)
             EcalT(:,:,-1)= ITm(:,:)
             EcalT(:,:,+1)= ITp(:,:)
             !   IF (Qz>Qx) THEN
             !   IF (Qz<0._wp) THEN
             !   ELSE
             VauxqzA= 0._wp   ! Initializes the integrand of the Qzp integrals
             VauxqzB= 0._wp
             do kqzp= 1,nqz
                Qzp= VQz(kqzp)
                Qzp2= Qzp**2
                do sigmap= -1,1,2
                   Asp= sigmap/sigma
                   Vauxqxp(:)= EcalS(:,kqzp,sigmap)
                   do sigmapp= -1,1,2
                      Aspp= sigmapp/sigma
                      do S1= -1,1,2
                         do S2= -1,1,2
                            Qxp= QxpLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            ires= IQxpLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            iresdif= IQxpdifLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            Qxp2= Qxp**2
                            Qxdif= Qx-S1*Qxp
                            Qzdif= Qz-S2*Qzp
                            S3= sign(1._wp,Qzdif)
                            Zlqp= ZL(Qxp,Qzp)
                            Ztqdif= ZT(Qxdif,abs(Qzdif))
                            D1= -S2*sigmap*AA*Qxp/sqrt(Qxp2+Qzp2+EpsMin) &
                                 / (1._wp+(Qxp2+Qzp2)/2._wp)**(1.5_wp) &
                                 + S1*S3*sigmapp*Qxdif/Ve2C2 &
                                 / sqrt(1._wp+(Qxdif**2+Qzdif**2)/Ve2C2)
                            AuxCoef= sqrt(Qxp2+Qzp2)*(Q2*(Qxp2+Qzp2) &
                                 - (S1*Qx*Qxp+S2*Qz*Qzp)**2) &
                                 /(Qxdif**2+Qzdif**2+EpsMin)/(abs(D1)+EpsMin)
                            iresdif2= kqz-S2*kqzp
                            Vauxqxdif(:)= EcalT(:,abs(iresdif2 &
                                 + S2*sign(1,iresdif2)),sigmapp)
                            call Aux_Coef_Decay_z(nqx,nqz,ires,iresdif,&
                                 Qxp,Qzp,Qxdif,Qzdif,"S","T",VQx,VQz, &
                                 Vauxqxp,Vauxqxdif,Iqp,Iqdif,nqcd,VQQ, &
                                 BremL1D,BremS1D,GcollL1D,GcollS1D)
                            VauxqzA(kqzp)= VauxqzA(kqzp) &
                                 + AuxCoef*sigma*Zlq*Iqp*Iqdif
                            VauxqzB(kqzp)= VauxqzB(kqzp) &
                                 -AuxCoef*(S2*sigmap*Zlqp*Iqdif &
                                 +S3*sigmapp*2._wp*Ztqdif*Iqp)
                         end do
                      end do
                   end do
                end do
             end do
             call Simpson(VQz,VauxqzA,nqz,Aux)
             CoefLSTdA= Aux
             call Simpson(VQz,VauxqzB,nqz,Aux)
             CoefLSTdB= Aux
             !   END IF
             CoefLSTdA= (AA/2._wp) * (sigma*Zlq/Q2) * CoefLSTdA
             CoefLSTdB= (AA/2._wp) * (sigma*Zlq/Q2) * CoefLSTdB
          else
             CoefLSTdA= 0._wp
             CoefLSTdB= 0._wp
          end if

          if(LdecayTT=="Yes") then
             ! Contribution due to spontaneous and induced decay involving two T waves:
             EcalT(:,:,-1)= ITm(:,:)
             EcalT(:,:,+1)= ITp(:,:)
             !   IF (Qz>Qx) THEN
             !   IF (Qz<0._wp) THEN
             !   ELSE
             VauxqzA= 0._wp   ! Initializes the integrand of the Qzp integrals
             VauxqzB= 0._wp
             do kqzp= 1,nqz
                Qzp= VQz(kqzp)
                Qzp2= Qzp**2
                do sigmap= -1,1,2
                   Asp= sigmap/sigma
                   Vauxqxp(:)= EcalT(:,kqzp,sigmap)
                   do sigmapp= -1,1,2
                      Aspp= sigmapp/sigma
                      do S1= -1,1,2
                         do S2= -1,1,2
                            Qxp= QxpLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            ires= IQxpLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            iresdif= IQxpdifLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            Qxp2= Qxp**2
                            Qxdif= Qx-S1*Qxp
                            Qzdif= Qz-S2*Qzp
                            S3= sign(1._wp,Qzdif)
                            Ztqp= ZT(Qxp,Qzp)
                            Ztqdif= ZT(Qxdif,abs(Qzdif))
                            D1= -S2*sigmap*Qxp/Ve2C2/sqrt(1._wp+(Qxp2+Qzp2)/Ve2C2) &
                                 + S1*S3*sigmapp*Qxdif/Ve2C2 &
                                 / sqrt(1._wp+(Qxdif**2+Qzdif**2)/Ve2C2)
                            AuxCoef= (Q2+(S1*Qxp*Qxdif+S2*Qzp*Qzdif)**2 &
                                 /(Qxdif**2+Qzdif**2+EpsMin)) &
                                 / (abs(D1)+EpsMin)/Ztqp**2/Ztqdif**2
                            iresdif2= kqz-S2*kqzp
                            Vauxqxdif(:)= EcalT(:,abs(iresdif2 &
                                 + S2*sign(1,iresdif2)),sigmapp)
                            call Aux_Coef_Decay_z(nqx,nqz,ires,iresdif,&
                                 Qxp,Qzp,Qxdif,Qzdif,"T","T",VQx,VQz, &
                                 Vauxqxp,Vauxqxdif,Iqp,Iqdif,nqcd,VQQ, &
                                 BremL1D,BremS1D,GcollL1D,GcollS1D)
                            VauxqzA(kqzp)= VauxqzA(kqzp) &
                                 + AuxCoef*sigma*Zlq*Iqp*Iqdif
                            VauxqzB(kqzp)= VauxqzB(kqzp) &
                                 - AuxCoef*(S2*sigmap*2._wp*Ztqp*Iqdif &
                                 +S3*sigmapp*2._wp*Ztqdif*Iqp)
                         end do
                      end do
                   end do
                end do
             end do
             call Simpson(VQz,VauxqzA,nqz,Aux)
             CoefLTTdA= Aux
             call Simpson(VQz,VauxqzB,nqz,Aux)
             CoefLTTdB= Aux
             !   END IF
             CoefLTTdA= (sigma*Zlq/16._wp) * CoefLTTdA
             CoefLTTdB= (sigma*Zlq/16._wp) * CoefLTTdB
          else
             CoefLTTdA= 0._wp
             CoefLTTdB= 0._wp
          end if

          if(LscatLL=="Yes") then
             ! Contribution due to spontaneous and induced scattering involving L waves:
             EcalL(:,:,-1)= ILm(:,:)
             EcalL(:,:,+1)= ILp(:,:)
             ! Electron contribution:
             VintqxA= 0._wp
             VintqxB= 0._wp
             do i= 1,nqx
                Qxp= VQx(i)
                VintqzA= 0._wp
                VintqzB= 0._wp
                do k= 1,nqz
                   Qzp= VQz(k)
                   Zlqp= ZL(Qxp,Qzp)
                   do sigmap= -1,1,2
                      do S1= -1,1,2
                         Qxdif= Qx-S1*Qxp
                         S2= sigma*sigmap !The other value of S2 gives negligible contribution.
                         Qzdif= Qz-S2*Qzp
                         if((Qxdif**2+Qzdif**2)<=1.e-6_wp) then
                            Auxe= 0._wp
                         else
                            Aux0= (S1*Qx*Qxp+S2*Qz*Qzp)**2 &
                                 /(Qxp**2+Qzp**2+EpsMin) &
                                 /sqrt(Qxdif**2+Qzdif**2)
                            Aux= (sigma*Zlq-S2*sigmap*Zlqp)**2 &
                                 /(Qxdif**2+Qzdif**2)
                            Auxe= Aux0 * exp(-Aux/VRTeTs(m))*sqrt(1._wp/VRTeTs(m))
                         end if
                         VintqzA(k)= VintqzA(k) + (Auxe*AuxScatElSpo) &
                              * (AuxSE*Geff*(sigma*Zlq)*EcalL(i,k,sigmap))
                         VintqzB(k)= VintqzB(k) &
                              + Auxe * (AuxSE*AuxScatElSpo) &
                              * (AuxSE*Geff*(-S2*sigmap*Zlqp)) 
                      end do
                   end do
                end do
                call Simpson(VQz,VintqzA,nqz,Aux)
                VintqxA(i)= Aux
                call Simpson(VQz,VintqzB,nqz,Aux)
                VintqxB(i)= Aux
             end do
             call Simpson(VQx,VintqxA,nqx,Aux)
             CoefLLLsA= 1._wp/sqrt(Pi)/VRNeNs(m) * (sigma*Zlq/Q2) * Aux
             call Simpson(VQx,VintqxB,nqx,Aux)
             CoefLLLsB= 1._wp/sqrt(Pi)/VRNeNs(m) * (sigma*Zlq/Q2) * Aux
             ! Ion contribution, spontaneous scattering:
             Qstar= Q
             VintPhipA= 0._wp
             VintPhipB= 0._wp
             do i= 1,nph
                Phip= VPhip(i)
                Qxp= Qstar*sin(Phip)
                Qzp= Qstar*cos(Phip)
                if( (Qxp>=VQx(nqx) .or. Qxp<=VQx(1)) .or. &
                     (Qzp>=VQz(nqz) .or. Qzp<=VQz(1)) ) then
                   Qp2= Qxp**2+Qzp**2
                   Qp= sqrt(Qp2)
                   if(Qp<=5.e-3_wp) Qp=5.e-3_wp
                   call Locate(VQQ,nqcd,Qp,ires)
                   call Aitp1d2(nqcd,VQQ,BremL1D,Qp,Aux,ires)
                   BremL0= Aux
                   call Aitp1d2(nqcd,VQQ,BremS1D,Qp,Aux,ires)
                   BremS0= Aux
                   call Aitp1d2(nqcd,VQQ,GcollL1D,Qp,Aux,ires)
                   GcollL0= Aux
                   call Aitp1d2(nqcd,VQQ,GcollS1D,Qp,Aux,ires)
                   GcollS0= Aux
                   call Iwave_Init(Qx,Qz,BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
                   Iqpp= ILinit
                   Iqpm= ILinit
                else
                   ires= IQxLLL1(iqx,kqz,i)
                   kres= IQzLLL1(iqx,kqz,i)
                   call Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxp,Qzp,Iqpp,ires,kres)
                   call Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxp,Qzp,Iqpm,ires,kres)
                end if
                do S1= -1,1,2
                   sigmap= 1._wp
                   Aux= (cos(Phi-S1*sigma*sigmap*Phip))**2
                   VintPhipA(i)= VintPhipA(i)+Aux*Iqpp
                   VintPhipB(i)= VintPhipB(i)-Aux
                   sigmap= -1._wp
                   Aux= (cos(Phi-S1*sigma*sigmap*Phip))**2
                   VintPhipA(i)= VintPhipA(i)+Aux*Iqpm
                   VintPhipB(i)= VintPhipB(i)-Aux
                end do
             end do
             call Simpson(VPhip,VintPhipA,nph,Aux)
             CoefLLLsA= CoefLLLsA + (2._wp/3._wp)/VRNeNs(m)*(Zlq)**3*AuxSE*Geff*Aux
             call Simpson(VPhip,VintPhipB,nph,Aux)
             CoefLLLsB= CoefLLLsB + (2._wp/3._wp)/VRNeNs(m)*(Zlq)**3*AuxSE*Geff*Aux
             ! Induced scattering due to the ions:
             VintPhipB= 0._wp
             do i= 1,nph
                Phip= VPhip(i)
                Sphip= sin(Phip)
                Cphip= cos(Phip)
                do sigmap= -1,1,2
                   do S1= -1,1,2
                      Aux1= (Q2+Qstar**2-2._wp*S1*Q*Qstar &
                           * cos(Phi-S1*sigma*sigmap*Phip))
                      if (Aux1<1.e-6_wp) then
                         AuxNum= sqrt(1.e-6_wp)
                         Aalpha= RMiMe/VRTiTs(m)*(9._wp/4._wp)*(Qstar/Zlq)**2/1.e-6_wp 
                      else
                         AuxNum= sqrt(Aux1)
                         Aalpha= RMiMe/VRTiTs(m)*(9._wp/4._wp)*(Qstar/Zlq)**2/Aux1 
                      end if
                      Aux= (cos(Phi-S1*sigma*sigmap*Phip))**2*AuxNum
                      do SS= -1,1,2
                         Qs= Qstar+SS/sqrt(2._wp*Aalpha)
                         Qxp= Qs*Sphip
                         Qzp= Qs*Cphip
                         if( (Qxp>=VQx(nqx) .or. Qxp<=VQx(1)) .or. &
                              (Qzp>=VQz(nqz) .or. Qzp<=VQz(1)) ) then
                            Qp2= Qxp**2+Qzp**2
                            Qp= sqrt(Qp2)
                            if(Qp<=5.E-3_wp) Qp=5.E-3_wp
                            call Locate(VQQ,nqcd,Qp,ires)
                            call Aitp1d2(nqcd,VQQ,BremL1D,Qp,Aux,ires)
                            BremL0= Aux
                            call Aitp1d2(nqcd,VQQ,BremS1D,Qp,Aux,ires)
                            BremS0= Aux
                            call Aitp1d2(nqcd,VQQ,GcollL1D,Qp,Aux,ires)
                            GcollL0= Aux
                            call Aitp1d2(nqcd,VQQ,GcollS1D,Qp,Aux,ires)
                            GcollS0= Aux
                            call Iwave_Init(Qx,Qz,BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
                            Iqp= ILinit
                         else
                            S2= sigma*sigmap
                            ires= IQxLLL2(iqx,kqz,i,S2,S1,SS)
                            kres= IQzLLL2(iqx,kqz,i,S2,S1,SS)
                            if(sigmap==1) then
                               call Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxp,Qzp,Iqp,ires,kres)
                            else
                               call Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxp,Qzp,Iqp,ires,kres)
                            end if
                         end if
                         VintPhipB(i)= VintPhipB(i)+Aux*Iqp * SS*Qs/Qstar
                      end do !SS
                   end do !S1
                end do !sigmap
             end do
             call Simpson(VPhip,VintPhipB,nph,Aux)
             CoefLLLsB= CoefLLLsB + (4._wp/3._wp)/VRNeNs(m)/sqrt(VRTiTs(m)*RMiMe) &
                  * exp(-0.5_wp)/sqrt(2._wp) * (Zlq)**2 * Aux
          else
             CoefLLLsA= 0._wp
             CoefLLLsB= 0._wp
          end if

          if(LscatLT=="Yes") then
             ! Contribution due to spontaneous and induced scattering involving L and T 
             ! waves:
             EcalT(:,:,-1)= ITm(:,:)
             EcalT(:,:,+1)= ITp(:,:)
             ! Electron and ion contributions:
             Qstar= sqrt(3._wp/2._wp*Ve2C2)*Q
             VintPhipA= 0._wp
             VintPhipB= 0._wp
             do i= 1,nph
                Phip= VPhip(i)
                Qxp= Qstar*sin(Phip)
                Qzp= Qstar*cos(Phip)
                if( (Qxp>=VQx(nqx) .or. Qxp<=VQx(1)) .or. &
                     (Qzp>=VQz(nqz) .or. Qzp<=VQz(1)) ) then
                   Qp2= Qxp**2+Qzp**2
                   Qp= sqrt(Qp2)
                   if(Qp<=5.e-3_wp) Qp=5.e-3_wp
                   call Locate(VQQ,nqcd,Qp,ires)
                   call Aitp1d2(nqcd,VQQ,BremL1D,Qp,Aux,ires)
                   BremL0= Aux
                   call Aitp1d2(nqcd,VQQ,BremS1D,Qp,Aux,ires)
                   BremS0= Aux
                   call Aitp1d2(nqcd,VQQ,GcollL1D,Qp,Aux,ires)
                   GcollL0= Aux
                   call Aitp1d2(nqcd,VQQ,GcollS1D,Qp,Aux,ires)
                   GcollS0= Aux
                   call Iwave_Init(Qxp,Qzp,BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
                   Iqpp= ITinit
                   Iqpm= ITinit
                else
                   ires= IQxLLT1(iqx,kqz,i)
                   kres= IQzLLT1(iqx,kqz,i)
                   call Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxp,Qzp,Iqpp,ires,kres)
                   call Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxp,Qzp,Iqpm,ires,kres)
                end if
                do S1= -1,1,2
                   sigmap= 1._wp
                   Aux= (sin(Phi-S1*sigma*sigmap*Phip))**2
                   VintPhipA(i)= VintPhipA(i)+Aux*Iqpp/2._wp
                   VintPhipB(i)= VintPhipB(i)-Aux
                   sigmap= -1._wp
                   Aux= (sin(Phi-S1*sigma*sigmap*Phip))**2
                   VintPhipA(i)= VintPhipA(i)+Aux*Iqpm/2._wp
                   VintPhipB(i)= VintPhipB(i)-Aux
                end do
             end do
             call Simpson(VPhip,VintPhipA,nph,Aux)
             CoefLLTsA= 2._wp*Ve2C2/VRNeNs(m) * (Zlq)**3 * AuxSE*Geff*Aux &
                  + 2._wp*Ve2C2/VRNeNs(m) * (Zlq)**3 * AuxSE*AuxScatElSpo*Geff*Aux
             call Simpson(VPhip,VintPhipB,nph,Aux)
             CoefLLTsB= 2._wp*Ve2C2/VRNeNs(m) * (Zlq)**3 * AuxSE*Geff*Aux &
                  + 2._wp*Ve2C2/VRNeNs(m) * (Zlq)**3 * AuxSE*AuxScatElSpo*Geff*Aux
             ! Induced scattering due to the ions:
             VintPhipB= 0._wp
             do i= 1,nph
                Phip= VPhip(i)
                Sphip= sin(Phip)
                Cphip= cos(Phip)
                do sigmap= -1,1,2
                   do S1= -1,1,2
                      Aux1= (Q2+Qstar**2-2._wp*S1*Q*Qstar &
                           * cos(Phi-S1*sigma*sigmap*Phip))
                      if (Aux1<1.e-6_wp) then
                         AuxNum= sqrt(1.e-6_wp)
                         Aalpha= RMiMe/VRTiTs(m)*(1._wp/Ve2C2**2) &
                              * (Qstar/Zlq)**2/1.e-6_wp
                      else
                         AuxNum= sqrt(Aux1)
                         Aalpha= RMiMe/VRTiTs(m)*(1._wp/Ve2C2**2) &
                              * (Qstar/Zlq)**2/Aux1 
                      end if
                      Aux= (sin(Phi-S1*sigma*sigmap*Phip))**2*sqrt(AuxNum)
                      do SS= -1,1,2
                         Qs= Qstar+SS/sqrt(2._wp*Aalpha)
                         Qxp= Qs*Sphip
                         Qzp= Qs*Cphip
                         if( (Qxp>=VQx(nqx) .or. Qxp<=VQx(1)) .or. &
                              (Qzp>=VQz(nqz) .or. Qzp<=VQz(1)) ) then
                            Qp2= Qxp**2+Qzp**2
                            Qp= sqrt(Qp2)
                            if(Qp<=5.e-3_wp) Qp=5.e-3_wp
                            call Locate(VQQ,nqcd,Qp,ires)
                            call Aitp1d2(nqcd,VQQ,BremL1D,Qp,Aux,ires)
                            BremL0= Aux
                            call Aitp1d2(nqcd,VQQ,BremS1D,Qp,Aux,ires)
                            BremS0= Aux
                            call Aitp1d2(nqcd,VQQ,GcollL1D,Qp,Aux,ires)
                            GcollL0= Aux
                            call Aitp1d2(nqcd,VQQ,GcollS1D,Qp,Aux,ires)
                            GcollS0= Aux
                            call Iwave_Init(Qxp,Qzp,BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
                            Iqp= ITinit
                         else
                            S2= sigma*sigmap
                            ires= IQxLLT2(iqx,kqz,i,S2,S1,SS)
                            kres= IQzLLT2(iqx,kqz,i,S2,S1,SS)
                            if(sigmap==1) then
                               call Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxp,Qzp,Iqp,ires,kres)
                            else
                               call Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxp,Qzp,Iqp,ires,kres)
                            end if
                         end if
                         VintPhipB(i)= VintPhipB(i)+Aux*Iqp/2._wp * SS*Qs/Qstar
                      end do !SS
                   end do !S1
                end do !sigmap
             end do
             call Simpson(VPhip,VintPhipB,nph,Aux)
             CoefLLTsB= CoefLLTsB + 2._wp*Ve2C2/VRNeNs(m)/sqrt(VRTiTs(m)*RMiMe) &
                  * exp(-0.5_wp)/sqrt(2._wp) * (Zlq)**2 * Aux
          else
             CoefLLTsA= 0._wp
             CoefLLTsB= 0._wp
          end if

          if(GcollEvol== "Yes")then
             GcollLB=GcollLp(iqx,kqz)!-GcollLm(iqx,kqz)
          else
             GcollLB=0._wp
          end if

          !Brem aqui

          if(BremssEvol== "Yes")then
             BremLA= BremLp(iqx,kqz)
          else
             BremLA= 0._wp
          end if

          CoefA(iqx,kqz)= CoefEA+CoefLLSdA+CoefLLTdA+CoefLSTdA+CoefLTTdA &
               +CoefLLLsA+CoefLLTsA+BremLA
          CoefB(iqx,kqz)= CoefEB+CoefLLSdB+CoefLLTdB+CoefLSTdB+CoefLTTdB &
               +CoefLLLsB+CoefLLTsB+2._wp*GcollLB
       end do
    end do

    ! CALL Output_Coef("Lwave",CoefA,CoefB)

    return
  end subroutine Coef_Lwave

  subroutine Coef_Swave(sigma,Dfdux,Dfduz,CoefA,CoefB)
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Ux,Qx,Qz,Qx2,Qz2,Q,Q2
    real(wp) :: Zlq,Zsq,Zsoq,Uresp,Uresm
    real(wp) :: Aux,AuxSE,D1,AuxCoef
    real(wp) :: CoefEA,CoefSLLdA,CoefSLTdA
    real(wp) :: CoefEB,CoefSLLdB,CoefSLTdB
    real(wp) :: BremSA,GcollSB
    real(wp) :: Zlqp,Zlqdif,Ztqdif
    real(wp) :: Muq
    real(wp) :: Qxp,Qzp
    real(wp) :: Qxp2,Qzp2
    real(wp) :: Qxdif,Qzdif
    real(wp) :: Iqp,Iqdif
    real(wp) :: Fesum,Femax,Fef,Feb
    ! real(wp) :: Fekappa
    real(wp) :: Uz
    real(wp) :: Beta
    real(wp), dimension(nux) :: VintuxA,VauxuxA
    real(wp), dimension(nux) :: VintuxB,VauxuxB
    real(wp), dimension(nuz) :: VauxuzA,VintuzA
    real(wp), dimension(nuz) :: VauxuzB,VintuzB
    real(wp), dimension(nux,nuz) :: Dfdux,Dfduz
    real(wp), dimension(nqz) :: VauxqzA,VauxqzB
    real(wp), dimension(nqx) :: Vauxqxp,Vauxqxdif
    real(wp), dimension(nqx,nqz) :: CoefA,CoefB
    real(wp), dimension(nqx,nqz,-1:1) :: EcalL
    real(wp), dimension(nqx,nqz,-1:1) :: EcalT
    integer :: l,sigma,sigmap,sigmapp
    integer :: iqx,kqz,kqzp
    integer :: iresp,iresm
    integer :: ires,iresdif,iresdif2
    integer :: m
    integer :: S1,S2,S3
    integer :: Asp,Aspp

    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    if(SpontEmis=="Yes") then
       AuxSE= 1.E0_wp
    else
       AuxSE= 0.E0_wp
    end if

    BremSA= 0._wp
    GcollSB= 0._wp
    CoefA= 0._wp
    CoefB= 0._wp
    do iqx= 1,nqx
       Qx= VQx(iqx)
       do kqz= 1,nqz
          Qz= VQz(kqz)
          Qx2= Qx**2
          Qz2= Qz**2
          Q2= Qx2+Qz2
          Q= sqrt(Q2)
          Muq= Q**3*AA/2._wp
          Zlq= ZL(Qx,Qz)
          Zsq= ZS(Qx,Qz)
          Zsoq= AA*sqrt(VRTeTs(m))/sqrt(1._wp+Q2/2._wp*VRTeTs(m)/VRNeNs(m))

          if(Semis=="Yes") then
             ! Contribution due to spontaneous and induced emission:
             if (Qz>Qx) then
                do l= 1,nux
                   Ux= VUx(l)
                   Uresp= UzrpSql(iqx,kqz,l,sigma)
                   Uresm= UzrmSql(iqx,kqz,l,sigma)
                   iresp= IUzrpSql(iqx,kqz,l,sigma)
                   iresm= IUzrmSql(iqx,kqz,l,sigma)
                   if (iresp==0 .or. iresp==nuz) then
                      Uz= Uresp
                      call Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
                      VintuxA(l)= AuxSE * VRNeNs(m)*Geff*Fesum
                      ! VintuxB(l)= -(sigma*Zlq) &
                      ! * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) & 
                      ! - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                      ! * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m)/(Kappae-1._wp))**(-1._wp)*Fef &
                      ! - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                      ! * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m)/(Kappae-1._wp))**(-1._wp)*Feb &
                      ! + Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                      ! + Qz*2._wp*(Uz-Uf)*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                      ! * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m)/(Kappae-1._wp))**(-1._wp)*Fef &
                      ! + Qz*2._wp*(Uz-Ub)*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                      ! * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m)/(Kappae-1._wp))**(-1._wp)*Feb )
                      VintuxB(l)= (sigma*Zlq) &
                           * ( Qx*2._wp*Ux*Femax/VRTeTs(m)+Qx*2._wp*Ux*Fef/VRTfTs(m) &
                           + Qx*2._wp*Ux*Feb/VRTbTs(m) &
                           - Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m)-Qz*2._wp*(Uz-Uf)*Fef/VRTfTs(m) &
                           - Qz*2._wp*(Uz-Ub)*Feb/VRTbTs(m) )
                   else
                      VauxuzA(:)= AuxSE * VRNeNs(m)*Geff*Fe(l,:)
                      VauxuzB(:)= (sigma*Zlq)*(-Qx*Dfdux(l,:) &
                           +Qz*Dfduz(l,:))
                      call Aitp1d2(nuz,VUz,VauxuzA,Uresp,Aux,iresp)
                      VintuxA(l)= Aux
                      call Aitp1d2(nuz,VUz,VauxuzB,Uresp,Aux,iresp)
                      VintuxB(l)= Aux
                   end if
                   if (iresm==0 .or. iresm==nuz) then
                      Uz= Uresm
                      call Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
                      VintuxA(l)= VintuxA(l) + AuxSE*VRNeNs(m)*Geff*Fesum
                      ! VintuxB(l)= -VintuxB(l)+ (sigma*Zlq) &
                      !  * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) & 
                      !  - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                      !  * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m)/(Kappae-1._wp))**(-1._wp)*Fef &
                      !  - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                      !  * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m)/(Kappae-1._wp))**(-1._wp)*Feb &
                      !  - Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                      !  - Qz*2._wp*(Uz-Uf)*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                      !  * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m)/(Kappae-1._wp))**(-1._wp)*Fef &
                      !  - Qz*2._wp*(Uz-Ub)*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                      !  * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m)/(Kappae-1._wp))**(-1._wp)*Feb )
                      VintuxB(l)= VintuxB(l) + (sigma*Zlq) &
                           * ( -Qx*2._wp*Ux*Femax/VRTeTs(m)-Qx*2._wp*Ux*Fef/VRTfTs(m) &
                           - Qx*2._wp*Ux*Feb/VRTbTs(m) &
                           - Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m)-Qz*2._wp*(Uz-Uf)*Fef/VRTfTs(m) &
                           - Qz*2._wp*(Uz-Ub)*Feb/VRTbTs(m) )
                   else
                      VauxuzA(:)= AuxSE * VRNeNs(m)*Geff*Fe(l,:)
                      VauxuzB(:)= (sigma*Zlq)*(Qx*Dfdux(l,:) &
                           +Qz*Dfduz(l,:))
                      call Aitp1d2(nuz,VUz,VauxuzA,Uresm,Aux,iresm)
                      VintuxA(l)= VintuxA(l) + Aux
                      call Aitp1d2(nuz,VUz,VauxuzB,Uresm,Aux,iresm)
                      VintuxB(l)= VintuxB(l) + Aux
                   end if
                end do
                call Simpson(VUx,VintuxA,nux,CoefEA)
                call Simpson(VUx,VintuxB,nux,CoefEB)
                CoefEA= CoefEA/(abs(Qz))
                CoefEB= CoefEB/(abs(Qz))
             else
                do l= 1,nuz
                   Uz= VUz(l)
                   Uresp= UxrpSql(iqx,kqz,l,sigma)
                   Uresm= UxrmSql(iqx,kqz,l,sigma)
                   iresp= IUxrpSql(iqx,kqz,l,sigma)
                   iresm= IUxrmSql(iqx,kqz,l,sigma)
                   if (iresm==0 .or. iresm==nux) then
                      if (iresm==nux) then
                         Ux= Uresm
                         call Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
                         VintuzA(l)= AuxSE * VRNeNs(m)*Geff*Fesum
                         ! VintuzB(l)= -(sigma*Zlq) &
                         !  * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) & 
                         !  - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                         !  * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m)/(Kappae-1._wp))**(-1._wp)*Fef &
                         !  - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                         !  * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m)/(Kappae-1._wp))**(-1._wp)*Feb &
                         !  + Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                         !  + Qz*2._wp*(Uz-Uf)*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                         !  * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m)/(Kappae-1._wp))**(-1._wp)*Fef &
                         !  + Qz*2._wp*(Uz-Ub)*(Kappae)/(Kappae-1._wp5_wp)/VRTbTs(m) &
                         !  * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m)/(Kappae-1._wp))**(-1._wp)*Feb )
                         VintuzB(l)= (sigma*Zlq) &
                              * ( Qx*2._wp*Ux*Femax/VRTeTs(m)+Qx*2._wp*Ux*Fef/VRTfTs(m) &
                              + Qx*2._wp*Ux*Feb/VRTbTs(m) &
                              - Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m)-Qz*2._wp*(Uz-Uf)*Fef/VRTfTs(m) &
                              - Qz*2._wp*(Uz-Ub)*Feb/VRTbTs(m) )
                      else
                         VintuzA(l)= 0.0_wp
                         VintuzB(l)= 0.0_wp
                      end if
                   else
                      VauxuxA(:)= AuxSE * VRNeNs(m)*Geff*Fe(:,l)
                      VauxuxB(:)= (sigma*Zlq)*(-Qx*Dfdux(:,l) &
                           +Qz*Dfduz(:,l))
                      call Aitp1d2(nux,VUx,VauxuxA,Uresm,Aux,iresm)
                      VintuzA(l)= Aux
                      call Aitp1d2(nux,VUx,VauxuxB,Uresm,Aux,iresm)
                      VintuzB(l)= Aux
                   end if
                   if (iresp==0 .or. iresp==nux) then
                      if (iresp==nux) then
                         Ux= Uresp
                         call Fe_Init(Ux,Uz,Fesum,Femax,Fef,Feb)
                         VintuzA(l)= VintuzA(l) + AuxSE * VRNeNs(m)*Geff*Fesum
                         ! VintuzB(l)= VintuzB(l)+ (sigma*Zlq) &
                         !  * ( -Qx*2._wp*Ux*Femax/VRTeTs(m) & 
                         !  - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                         !  * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m)/(Kappae-1._wp))**(-1._wp)*Fef &
                         !  - Qx*2._wp*Ux*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                         !  * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m)/(Kappae-1._wp))**(-1._wp)*Feb &
                         !  - Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m) &
                         !  - Qz*2._wp*(Uz-Uf)*(Kappae)/(Kappae-1._wp)/VRTfTs(m) &
                         !  * (1._wp+(Ux**2+(Uz-Uf)**2)/VRTfTs(m)/(Kappae-1._wp))**(-1._wp)*Fef &
                         !  - Qz*2._wp*(Uz-Ub)*(Kappae)/(Kappae-1._wp)/VRTbTs(m) &
                         !  * (1._wp+(Ux**2+(Uz-Ub)**2)/VRTbTs(m)/(Kappae-1._wp))**(-1._wp)*Feb )
                         VintuzB(l)= VintuzB(l) + (sigma*Zlq) &
                              * ( -Qx*2._wp*Ux*Femax/VRTeTs(m)-Qx*2._wp*Ux*Fef/VRTfTs(m) &
                              - Qx*2._wp*Ux*Feb/VRTbTs(m) &
                              - Qz*2._wp*(Uz-U0)*Femax/VRTeTs(m)-Qz*2._wp*(Uz-Uf)*Fef/VRTfTs(m) &
                              - Qz*2._wp*(Uz-Ub)*Feb/VRTbTs(m) )
                      else
                         VintuzA(l)= VintuzA(l)+ 0.0_wp
                         VintuzB(l)= VintuzB(l)+ 0.0_wp
                      end if
                   else
                      VauxuxA(:)= AuxSE * VRNeNs(m)*Geff*Fe(:,l)
                      VauxuxB(:)= (sigma*Zlq)*(Qx*Dfdux(:,l) &
                           +Qz*Dfduz(:,l))
                      call Aitp1d2(nux,VUx,VauxuxA,Uresp,Aux,iresp)
                      VintuzA(l)= VintuzA(l) + Aux
                      call Aitp1d2(nux,VUx,VauxuxB,Uresp,Aux,iresp)
                      VintuzB(l)= VintuzB(l) + Aux
                   end if
                end do
                call Simpson(VUz,VintuzA,nuz,CoefEA)
                call Simpson(VUz,VintuzB,nuz,CoefEB)
                CoefEA= CoefEA/(abs(Qx))
                CoefEB= CoefEB/(abs(Qx))
             end if
             CoefEA= Pi * Q*AA/2._wp &
                  * (VRTeTs(m))**2/sqrt(VRTeTs(m)*VRNeNs(m)) * CoefEA
             CoefEB= Pi * Q*AA/2._wp &
                  * (VRTeTs(m))**2/sqrt(VRTeTs(m)*VRNeNs(m)) * CoefEB
             Aux= sqrt(Pi)*AA/2._wp * VRTeTs(m)/sqrt(VRNeNs(m)) &
                  * sqrt(RMiMe*RTeTi)*exp(-RMiMe/VRTiTs(m)*(Zsoq)**2)
             CoefEA= CoefEA + Aux*(AuxSE * VRNeNs(m)*Geff)
             CoefEB= CoefEB + Aux*(- 2._wp/VRTiTs(m)*Zlq*Zsq)
          else
             CoefEA= 0._wp
             CoefEB= 0._wp
          end if

          if(SdecayLL=="Yes") then
             ! Contribution due to spontaneous and induced decay:
             EcalL(:,:,-1)= ILm(:,:)
             EcalL(:,:,+1)= ILp(:,:)
             !   IF (Qz>Qx) THEN
             !   IF (Qz<0._wp) THEN
             !   ELSE
             VauxqzA= 0._wp   ! Initializes the integrand of the Qzp integrals
             VauxqzB= 0._wp
             do kqzp= 1,nqz
                Qzp= VQz(kqzp)
                Qzp2= Qzp**2
                do sigmap= -1,1,2
                   Asp= sigmap/sigma
                   Vauxqxp(:)= EcalL(:,kqzp,sigmap)
                   do sigmapp= -1,1,2
                      Aspp= sigmapp/sigma
                      do S1= -1,1,2
                         do S2= -1,1,2
                            Qxp= QxpSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            ires= IQxpSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            iresdif= IQxpdifSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            Qxp2= Qxp**2
                            Qxdif= Qx-S1*Qxp
                            Qzdif= Qz-S2*Qzp
                            S3= sign(1._wp,Qzdif)
                            Zlqp= ZL(Qxp,Qzp)
                            Zlqdif= ZL(Qxdif,abs(Qzdif))
                            D1= -S2*sigmap*1.5_wp*Beta*Qxp &
                                 /sqrt(1._wp+1.5_wp*Beta*(Qxp2+Qzp2)) &
                                 + S1*S3*sigmapp*1.5_wp*Beta*Qxdif &
                                 /sqrt(1._wp+1.5_wp*Beta*(Qxdif**2+Qzdif**2))
                            AuxCoef= (Qxp*Qxdif+Qzp*Qzdif)**2 &
                                 /(Qxp2+Qzp2+EpsMin)/(Qxdif**2+Qzdif**2+EpsMin)/(abs(D1)+EpsMin)
                            iresdif2= kqz-S2*kqzp
                            Vauxqxdif(:)= EcalL(:,abs(iresdif2+S2*sign(1,iresdif2)),sigmapp)
                            call Aux_Coef_Decay_z(nqx,nqz,ires,iresdif,&
                                 Qxp,Qzp,Qxdif,Qzdif,"L","L",VQx,VQz,Vauxqxp,Vauxqxdif,&
                                 Iqp,Iqdif,nqcd,VQQ,BremL1D,BremS1D,GcollL1D,GcollS1D)
                            VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*Zlq*Iqp*Iqdif
                            VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
                                 +S3*sigmapp*Zlqdif*Iqp)
                         end do
                      end do
                   end do
                end do
             end do
             call Simpson(VQz,VauxqzA,nqz,Aux)
             CoefSLLdA= Aux
             call Simpson(VQz,VauxqzB,nqz,Aux)
             CoefSLLdB= Aux
             !   END IF
             CoefSLLdA= (AA/2._wp)*(sigma*Zlq)*Q &
                  /VRNeNs(m)/sqrt(VRNeNs(m)*VRTeTs(m)) * CoefSLLdA
             CoefSLLdB= (AA/2._wp)*(sigma*Zlq)*Q &
                  /VRNeNs(m)/sqrt(VRNeNs(m)*VRTeTs(m)) * CoefSLLdB
          else
             CoefSLLdA= 0._wp
             CoefSLLdB= 0._wp
          end if

          if(SdecayLT=="Yes") then
             ! Contribution due to spontaneous and induced decay involving L and T waves:
             EcalL(:,:,-1)= ILm(:,:)
             EcalL(:,:,+1)= ILp(:,:)
             EcalT(:,:,-1)= ITm(:,:)
             EcalT(:,:,+1)= ITp(:,:)
             !   IF (Qz>Qx) THEN
             !   IF (Qz<0._wp) THEN
             !   ELSE
             VauxqzA= 0._wp   ! Initializes the integrand of the Qzp integrals
             VauxqzB= 0._wp
             do kqzp= 1,nqz
                Qzp= VQz(kqzp)
                Qzp2= Qzp**2
                do sigmap= -1,1,2
                   Asp= sigmap/sigma
                   Vauxqxp(:)= EcalL(:,kqzp,sigmap)
                   do sigmapp= -1,1,2
                      Aspp= sigmapp/sigma
                      do S1= -1,1,2
                         do S2= -1,1,2
                            Qzdif= Qz-S2*Qzp
                            Qxp= QxpSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            ires= IQxpSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            iresdif= IQxpdifSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            Qxp2= Qxp**2
                            Qxdif= Qx-S1*Qxp
                            S3= sign(1._wp,Qzdif)
                            Zlqp= ZL(Qxp,Qzp)
                            Ztqdif= ZT(Qxdif,abs(Qzdif))
                            D1= -S2*sigmap*(1.5_wp*Qxp/sqrt(1._wp+1.5_wp*(Qxp2+Qzp2)) &
                                 +S1*S3*sigmapp*Qxdif/Ve2C2/sqrt(1._wp+(Qxdif**2+Qzdif**2)/Ve2C2))
                            AuxCoef= (Q2*(Qxp2+Qzp2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
                                 /(Qxp2+Qzp2+EpsMin)/(Qxdif**2+Qzdif**2+EpsMin)/(abs(D1)+EpsMin)
                            iresdif2= kqz-S2*kqzp
                            Vauxqxdif(:)= EcalT(:,abs(iresdif2+S2*sign(1,iresdif2)),sigmapp)
                            call Aux_Coef_Decay_z(nqx,nqz,ires,iresdif,&
                                 Qxp,Qzp,Qxdif,Qzdif,"L","T",VQx,VQz,Vauxqxp,Vauxqxdif,&
                                 Iqp,Iqdif,nqcd,VQQ,BremL1D,BremS1D,GcollL1D,GcollS1D)
                            VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*Zlq*Iqp*Iqdif
                            VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
                                 +S3*sigmapp*2._wp*Ztqdif*Iqp)
                         end do
                      end do
                   end do
                end do
             end do
             call Simpson(VQz,VauxqzA,nqz,Aux)
             CoefSLTdA= Aux
             call Simpson(VQz,VauxqzB,nqz,Aux)
             CoefSLTdB= Aux
             !   END IF
             CoefSLTdA= (AA/2._wp)/VRNeNs(m)/sqrt(VRNeNs(m)*VRTeTs(m)) &
                  *Q*(sigma*Zlq) * CoefSLTdA
             CoefSLTdB= (AA/2._wp)/VRNeNs(m)/sqrt(VRNeNs(m)*VRTeTs(m)) &
                  *Q*(sigma*Zlq) * CoefSLTdB
          else
             CoefSLTdA= 0._wp
             CoefSLTdB= 0._wp
          end if

          if(GcollEvol== "Yes")then
             GcollSB=GcollSp(iqx,kqz)
          else
             GcollSB=0._wp
          end if

          if(BremssEvol== "Yes")then
             BremSA= BremSp(iqx,kqz)
          else
             BremSA= 0._wp
          end if

          CoefA(iqx,kqz)= CoefEA+CoefSLLdA+CoefSLTdA+BremSA
          CoefB(iqx,kqz)= CoefEB+CoefSLLdB+CoefSLTdB+2._wp*GcollSB
       end do
    end do


    ! CALL Output_Coef("Swave",CoefA,CoefB)

    return
  end subroutine Coef_Swave

  subroutine Coef_Twave(sigma,CoefA,CoefB)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Qx,Qz
    real(wp) :: Qx2,Qz2,Q2,Q
    real(wp) :: Zlq,Ztq,Zlqp,Zlqdif,Ztqp
    real(wp) :: Aux,AuxSE,AuxCoef,Aux1
    !REAL(WP) :: Aux0,Auxe,Auxi
    real(wp) :: AuxScatElSpo
    !REAL(WP) :: AuxScatElInd
    real(wp) :: CoefTLLdA,CoefTLSdA,CoefTTLdA,CoefTLTsA
    real(wp) :: CoefTLLdB,CoefTLSdB,CoefTTLdB,CoefTLTsB
    real(wp) :: Qxp,Qzp,Qxp2,Qzp2
    real(wp) :: Qxdif,Qzdif
    real(wp) :: Qp,Qp2
    real(wp) :: Iqp,Iqdif,Iqpp,Iqpm
    real(wp) :: D1
    real(wp) :: ILinit,ISinit,ITinit
    real(wp) :: Phip,Phi,Sphip,Cphip
    real(wp) :: Aalpha,AuxNum
    real(wp) :: Qstar,Qs
    real(wp) :: BremL0,BremS0,GcollL0,GcollS0
    !REAL(WP), DIMENSION(nqx) :: VintqxA,VintqxB
    !REAL(WP), DIMENSION(nqz) :: VintqzA,VintqzB
    !REAL(WP), DIMENSION(nqz) :: VauxqzA,VauxqzB
    !REAL(WP), DIMENSION(nqx) :: VauxqxA,VauxqxB
    real(wp), dimension(nqx,nqz) :: CoefA,CoefB
    real(wp), dimension(nqx2,nqz2,-1:1) :: EcalL2
    real(wp), dimension(nqx2,nqz2,-1:1) :: EcalS2
    real(wp), dimension(nqx2,nqz2,-1:1) :: EcalT2
    real(wp), dimension(nqz2) :: VauxqzA,VauxqzB
    real(wp), dimension(nqx2) :: Vauxqxp,Vauxqxdif
    !REAL(WP), DIMENSION(nqz) :: Vauxqzp,Vauxqzdif
    real(wp), dimension(nqx2,nqz2) :: CoefAaux,CoefBaux
    real(wp), dimension(nph) :: VintPhipA,VintPhipB
    integer :: i,sigma
    integer :: m
    integer :: sigmap,sigmapp
    integer :: iqx,kqz,kqzp
    integer :: ires,iresdif,iresdif2,kres
    integer :: iresq,kresq
    integer :: S1,S2,S3,SS
    integer :: Asp,Aspp

    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    if(SpontEmis=="Yes") then
       AuxSE= 1.e0_wp
    else
       AuxSE= 0.e0_wp
    end if
    if(ScatElSpo=="Yes") then
       AuxScatElSpo= 1.e0_wp
    else
       AuxScatElSpo= 0.e0_wp
    end if
    !IF(ScatElInd=="Yes") THEN
    ! AuxScatElInd= 1.e0_wp
    !ELSE
    ! AuxScatElInd= 0.e0_wp
    !END IF

    do iqx= 1,nqx2
       Qx= VQx2(iqx)
       iresq= IQresx2(iqx)
       do kqz= 1,nqz2
          Qz= VQz2(kqz)
          kresq= IQresz2(kqz)
          call Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qx,Qz,Aux,iresq,kresq)
          EcalL2(iqx,kqz,-1)= Aux
          call Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qx,Qz,Aux,iresq,kresq)
          EcalL2(iqx,kqz,+1)= Aux
          call Aitp2d2(nqx,nqz,VQx,VQz,ISm,Qx,Qz,Aux,iresq,kresq)
          EcalS2(iqx,kqz,-1)= Aux
          call Aitp2d2(nqx,nqz,VQx,VQz,ISp,Qx,Qz,Aux,iresq,kresq)
          EcalS2(iqx,kqz,+1)= Aux
          call Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qx,Qz,Aux,iresq,kresq)
          EcalT2(iqx,kqz,-1)= Aux
          call Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qx,Qz,Aux,iresq,kresq)
          EcalT2(iqx,kqz,+1)= Aux
       end do
    end do

    ! Evaluation of the decay and scattering terms using the (nqx2,nqz2) grid, for 
    ! interpolation over the (nqx,nqz) grid:
    CoefAaux= 0._wp
    CoefBaux= 0._wp
    do iqx= 1,nqx2
       Qx= VQx2(iqx)
       do kqz= 1,nqz2
          Qz= VQz2(kqz)
          Qx2= Qx**2
          Qz2= Qz**2
          Q2= Qx2+Qz2
          Q= sqrt(Q2)
          Zlq= ZL(Qx,Qz)
          Ztq= ZT(Qx,Qz)
          Phi= acos(Qz/Q)

          if(TdecayLL=="Yes") then
             ! Contribution due to spontaneous and induced decay involving L waves:
             !   IF (Qz>Qx) THEN
             !   IF (Qz<0._wp) THEN
             !   ELSE
             VauxqzA= 0._wp   ! Initializes the integrand of the Qzp integrals
             VauxqzB= 0._wp
             do kqzp= 1,nqz2
                Qzp= VQz2(kqzp)
                Qzp2= Qzp**2
                do sigmap= -1,1,2
                   Asp= sigmap/sigma
                   Vauxqxp(:)= EcalL2(:,kqzp,sigmap)
                   do sigmapp= -1,1,2
                      Aspp= sigmapp/sigma
                      do S1= -1,1,2
                         do S2= -1,1,2
                            Qxp= QxpTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            ires= IQxpTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            iresdif= IQxpdifTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            Qxp2= Qxp**2
                            Qxdif= Qx-S1*Qxp
                            Qzdif= Qz-S2*Qzp
                            S3= sign(1._wp,Qzdif)
                            Zlqp= ZL(Qxp,Qzp)
                            Zlqdif= ZL(Qxdif,abs(Qzdif))
                            D1= -S2*sigmap*1.5_wp*Qxp/sqrt(1._wp+1.5_wp*(Qxp2+Qzp2)) &
                                 +S1*S3*sigmapp*1.5_wp*Qxdif/sqrt(1._wp+1.5_wp*(Qxdif**2+Qzdif**2))
                            AuxCoef= (Q2*(Qxp2+Qzp2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
                                 /(Qxp2+Qzp2+EpsMin)/(Qxdif**2+Qzdif**2+EpsMin) &
                                 *((Qxp2+Qzp2)/S2/sigmap/Zlqp &
                                 -(Qxdif**2+Qzdif**2)/S3/sigmapp/Zlqdif)**2/(abs(D1)+EpsMin)
                            iresdif2= kqz-S2*kqzp
                            Vauxqxdif(:)= EcalL2(:,abs(iresdif2+S2*sign(1,iresdif2)),sigmapp)
                            call Aux_Coef_Decay_z(nqx2,nqz2,ires,iresdif,&
                                 Qxp,Qzp,Qxdif,Qzdif,"L","L",VQx2,VQz2,Vauxqxp,Vauxqxdif,&
                                 Iqp,Iqdif,nqcd,VQQ,BremL1D,BremS1D,GcollL1D,GcollS1D)
                            VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*2._wp*Ztq*Iqp*Iqdif
                            VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
                                 +S3*sigmapp*Zlqdif*Iqp)
                         end do
                      end do
                   end do
                end do
             end do
             call Simpson(VQz2,VauxqzA,nqz2,Aux)
             CoefTLLdA= Aux
             call Simpson(VQz2,VauxqzB,nqz2,Aux)
             CoefTLLdB= Aux
             !   END IF
             CoefTLLdA= (1._wp/32._wp) * (sigma*Ztq/Q2) * CoefTLLdA
             CoefTLLdB= (1._wp/32._wp) * (sigma*Ztq/Q2) * CoefTLLdB
          else
             CoefTLLdA= 0._wp
             CoefTLLdB= 0._wp
          end if

          if(TdecayLS=="Yes") then
             ! Contribution due to spontaneous and induced decay involving L and S waves:
             !   IF (Qz>Qx) THEN
             !   IF (Qz<0._wp) THEN
             !   ELSE
             VauxqzA= 0._wp   ! Initializes the integrand of the Qzp integrals
             VauxqzB= 0._wp
             do kqzp= 1,nqz2
                Qzp= VQz2(kqzp)
                Qzp2= Qzp**2
                do sigmap= -1,1,2
                   Asp= sigmap/sigma
                   Vauxqxp(:)= EcalL2(:,kqzp,sigmap)
                   do sigmapp= -1,1,2
                      Aspp= sigmapp/sigma
                      do S1= -1,1,2
                         do S2= -1,1,2
                            Qxp= QxpTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            ires= IQxpTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            iresdif= IQxpdifTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            Qxp2= Qxp**2
                            Qxdif= Qx-S1*Qxp
                            Qzdif= Qz-S2*Qzp
                            S3= sign(1._wp,Qzdif)
                            Zlqp= ZL(Qxp,Qzp)
                            Zlqdif= ZL(Qxdif,abs(Qzdif))
                            D1= -S2*sigmap*(1.5_wp*Qxp/sqrt(1._wp+1.5_wp*(Qxp2+Qzp2))) &
                                 + S1*S3*sigmapp*AA*Qxdif/sqrt(Qxdif**2+Qzdif**2+EpsMin) &
                                 /(1._wp+(Qxdif**2+Qzdif**2)/2._wp)**(1.5_wp)
                            AuxCoef= sqrt((Qxdif)**2+(Qzdif)**2) &
                                 *(Q2*(Qxp2+Qzp2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
                                 /(Qxp2+Qzp2+EpsMin)/(abs(D1)+EpsMin)
                            iresdif2= kqz-S2*kqzp
                            Vauxqxdif(:)= EcalS2(:,abs(iresdif2+S2*sign(1,iresdif2)),sigmapp)
                            call Aux_Coef_Decay_z(nqx2,nqz2,ires,iresdif,&
                                 Qxp,Qzp,Qxdif,Qzdif,"L","S",VQx2,VQz2,Vauxqxp,Vauxqxdif,&
                                 Iqp,Iqdif,nqcd,VQQ,BremL1D,BremS1D,GcollL1D,GcollS1D)
                            VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*2._wp*Ztq*Iqp*Iqdif
                            VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
                                 +S3*sigmapp*Zlqdif*Iqp)
                         end do
                      end do
                   end do
                end do
             end do
             call Simpson(VQz2,VauxqzA,nqz2,Aux)
             CoefTLSdA= Aux
             call Simpson(VQz2,VauxqzB,nqz2,Aux)
             CoefTLSdB= Aux
             !   END IF
             CoefTLSdA= (AA/2._wp) * (sigma*Ztq/Q2) * CoefTLSdA
             CoefTLSdB= (AA/2._wp) * (sigma*Ztq/Q2) * CoefTLSdB
          else
             CoefTLSdA= 0._wp
             CoefTLSdB= 0._wp
          end if

          if(TdecayTL=="Yes") then
             ! Contribution due to spontaneous and induced decay involving L and S waves:
             !   IF (Qz>Qx) THEN
             !   IF (Qz<0._wp) THEN
             !   ELSE
             VauxqzA= 0._wp   ! Initializes the integrand of the Qzp integrals
             VauxqzB= 0._wp
             do kqzp= 1,nqz2
                Qzp= VQz2(kqzp)
                Qzp2= Qzp**2
                do sigmap= -1,1,2
                   Asp= sigmap/sigma
                   Vauxqxp(:)= EcalT2(:,kqzp,sigmap)
                   do sigmapp= -1,1,2
                      Aspp= sigmapp/sigma
                      do S1= -1,1,2
                         do S2= -1,1,2
                            Qxp= QxpTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            ires= IQxpTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            iresdif= IQxpdifTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
                            Qxp2= Qxp**2
                            Qxdif= Qx-S1*Qxp
                            Qzdif= Qz-S2*Qzp
                            S3= sign(1._wp,Qzdif)
                            Ztqp= ZT(Qxp,Qzp)
                            Zlqdif= ZL(Qxdif,abs(Qzdif))
                            D1= -S2*sigmap*(Qxp/Ve2C2/sqrt(1._wp+(Qxp2+Qzp2)/Ve2C2)) &
                                 + S1*S3*sigmapp*1.5_wp*Qxdif/sqrt(1._wp+1.5_wp*(Qxdif**2+Qzdif**2))
                            AuxCoef= ((Qxdif)**2+(Qzdif)**2)/Ztqp**2 &
                                 *(Q2+(S1*Qx*Qxp+S2*Qz*Qzp)**2/(Qxp2+Qzp2+EpsMin))/(abs(D1)+EpsMin)
                            iresdif2= kqz-S2*kqzp
                            Vauxqxdif(:)= EcalL2(:,abs(iresdif2+S2*sign(1,iresdif2)),sigmapp)
                            call Aux_Coef_Decay_z(nqx2,nqz2,ires,iresdif,&
                                 Qxp,Qzp,Qxdif,Qzdif,"T","L",VQx2,VQz2,Vauxqxp,Vauxqxdif,&
                                 Iqp,Iqdif,nqcd,VQQ,BremL1D,BremS1D,GcollL1D,GcollS1D)
                            VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*2._wp*Ztq*Iqp*Iqdif
                            VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*2._wp*Ztqp*Iqdif &
                                 +S3*sigmapp*Zlqdif*Iqp)
                         end do
                      end do
                   end do
                end do
             end do
             call Simpson(VQz2,VauxqzA,nqz2,Aux)
             CoefTTLdA= Aux
             call Simpson(VQz2,VauxqzB,nqz2,Aux)
             CoefTTLdB= Aux
             !   END IF
             CoefTTLdA= (sigma*Ztq/Ztq**2/8._wp/Q2) * CoefTTLdA
             CoefTTLdB= (sigma*Ztq/Ztq**2/8._wp/Q2) * CoefTTLdB
          else
             CoefTTLdA= 0._wp
             CoefTTLdB= 0._wp
          end if

          if(TscatLT=="Yes") then
             ! Contribution due to spontaneous and induced scattering involving L and T
             ! waves:
             ! Electron and Ion contributions:
             Qstar= sqrt(2._wp/3._wp/Ve2C2)*Q
             VintPhipA= 0._wp
             VintPhipB= 0._wp
             do i= 1,nph
                Phip= VPhip(i)
                Qxp= Qstar*sin(Phip)
                Qzp= Qstar*cos(Phip)
                if( (Qxp>=VQx(nqx) .or. Qxp<=VQx(1)) .or. &
                     (Qzp>=VQz(nqz) .or. Qzp<=VQz(1)) ) then
                   Qp2= Qxp**2+Qzp**2
                   Qp= sqrt(Qp2)
                   if(Qp<=5.e-3_wp) Qp=5.e-3_wp
                   call Locate(VQQ,nqcd,Qp,ires)
                   call Aitp1d2(nqcd,VQQ,BremL1D,Qp,Aux,ires)
                   BremL0= Aux
                   call Aitp1d2(nqcd,VQQ,BremS1D,Qp,Aux,ires)
                   BremS0= Aux
                   call Aitp1d2(nqcd,VQQ,GcollL1D,Qp,Aux,ires)
                   GcollL0= Aux
                   call Aitp1d2(nqcd,VQQ,GcollS1D,Qp,Aux,ires)
                   GcollS0= Aux
                   call Iwave_Init(Qxp,Qzp,BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
                   Iqpp= ILinit
                   Iqpm= ILinit
                else
                   ires= IQxTLT1(iqx,kqz,i)
                   kres= IQzTLT1(iqx,kqz,i)
                   call Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxp,Qzp,Iqpp,ires,kres)
                   call Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxp,Qzp,Iqpm,ires,kres)
                end if
                do S1= -1,1,2
                   sigmap= 1._wp
                   Aux= (sin(Phi-S1*sigma*sigmap*Phip))**2
                   VintPhipA(i)= VintPhipA(i)+Aux*Iqpp
                   VintPhipB(i)= VintPhipB(i)-Aux/2._wp
                   sigmap= -1._wp
                   Aux= (sin(Phi-S1*sigma*sigmap*Phip))**2
                   VintPhipA(i)= VintPhipA(i)+Aux*Iqpm
                   VintPhipB(i)= VintPhipB(i)-Aux/2._wp
                end do
             end do
             call Simpson(VPhip,VintPhipA,nph,Aux)
             CoefTLTsA= (2._wp/3._wp)/VRNeNs(m) * (Ztq)**3 * AuxSE*Geff*Aux &
                  + (2._wp/3._wp)/VRNeNs(m) * (Ztq)**3 * AuxSE*AuxScatElSpo*Geff*Aux 
             call Simpson(VPhip,VintPhipB,nph,Aux)
             CoefTLTsB= (2._wp/3._wp)/VRNeNs(m) * (Ztq)**3 * AuxSE*Geff*Aux &
                  + (2._wp/3._wp)/VRNeNs(m) * (Ztq)**3 * AuxSE*AuxScatElSpo*Geff*Aux 
             ! Induced scattering due to the ions:
             VintPhipB= 0._wp
             do i= 1,nph
                Phip= VPhip(i)
                Sphip= sin(Phip)
                Cphip= cos(Phip)
                do sigmap= -1,1,2
                   do S1= -1,1,2
                      Aux1= (Q2+Qstar**2-2._wp*S1*Q*Qstar*cos(Phi-S1*sigma*sigmap*Phip))
                      if (Aux1<1.e-6_wp) then
                         AuxNum= sqrt(1.e-6_wp)
                         Aalpha= RMiMe/VRTiTs(m)*(9._wp/4._wp)*(Qstar/Ztq)**2/1.e-6_wp
                      else
                         AuxNum= sqrt(Aux1)
                         Aalpha= RMiMe/VRTiTs(m)*(9._wp/4._wp)*(Qstar/Ztq)**2/Aux1 
                      end if
                      Aux= (sin(Phi-S1*sigma*sigmap*Phip))**2*sqrt(AuxNum)
                      do SS= -1,1,2
                         Qs= Qstar+SS/sqrt(2._wp*Aalpha)
                         Qxp= Qs*Sphip
                         Qzp= Qs*Cphip
                         if( (Qxp>=VQx(nqx) .or. Qxp<=VQx(1)) .or. &
                              (Qzp>=VQz(nqz) .or. Qzp<=VQz(1)) ) then
                            Qp2= Qxp**2+Qzp**2
                            Qp= sqrt(Qp2)
                            if(Qp<=5.e-3_wp) Qp=5.e-3_wp
                            call Locate(VQQ,nqcd,Qp,ires)
                            call Aitp1d2(nqcd,VQQ,BremL1D,Qp,Aux,ires)
                            BremL0= Aux
                            call Aitp1d2(nqcd,VQQ,BremS1D,Qp,Aux,ires)
                            BremS0= Aux
                            call Aitp1d2(nqcd,VQQ,GcollL1D,Qp,Aux,ires)
                            GcollL0= Aux
                            call Aitp1d2(nqcd,VQQ,GcollS1D,Qp,Aux,ires)
                            GcollS0= Aux
                            call Iwave_Init(Qxp,Qzp,BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
                            Iqp= ILinit
                         else
                            S2= sigma*sigmap
                            ires= IQxTLT2(iqx,kqz,i,S2,S1,SS)
                            kres= IQzTLT2(iqx,kqz,i,S2,S1,SS)
                            if(sigmap==1) then
                               call Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxp,Qzp,Iqp,ires,kres)
                            else
                               call Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxp,Qzp,Iqp,ires,kres)
                            end if
                         end if
                         VintPhipB(i)= VintPhipB(i)+Aux*Iqp * SS*Qs/Qstar
                      end do !SS
                   end do !S1
                end do !sigmap
             end do
             call Simpson(VPhip,VintPhipB,nph,Aux)
             CoefTLTsB= CoefTLTsB + (4._wp/3._wp)/VRNeNs(m)/sqrt(VRTiTs(m)*RMiMe) &
                  * exp(-0.5_wp)/sqrt(2._wp) * (Ztq)**2 * Aux/2._wp
          else
             CoefTLTsA= 0._wp
             CoefTLTsB= 0._wp
          end if

          CoefAaux(iqx,kqz)= CoefTLLdA+CoefTLSdA+CoefTTLdA+CoefTLTsA
          CoefBaux(iqx,kqz)= CoefTLLdB+CoefTLSdB+CoefTTLdB+CoefTLTsB
       end do
    end do

    CoefA= 0._wp
    CoefB= 0._wp
    do iqx= 1,nqx
       Qx= VQx(iqx)
       iresq= IQresx(iqx)
       do kqz= 1,nqz
          Qz= VQz(kqz)
          kresq= IQresz(kqz)

          ! Interpolation of the contribution due to decay terms (TLLd, TLSd, TTLd),
          ! and the contribution of spontaneous and induced scattering involving L and
          ! T waves (TLTs):
          ! Routine "Aitp2d2b" corrects the interpolation at the extreme points.
          call Aitp2d2b(nqx2,nqz2,VQx2,VQz2,CoefAaux,Qx,Qz,Aux,iresq,kresq)
          CoefA(iqx,kqz)= Aux
          call Aitp2d2b(nqx2,nqz2,VQx2,VQz2,CoefBaux,Qx,Qz,Aux,iresq,kresq)
          CoefB(iqx,kqz)= Aux
       end do
    end do

    ! CALL Output_Coef("Twave",CoefA,CoefB)

    return
  end subroutine Coef_Twave

  subroutine Aux_Coef_Decay_z(nqx,nqz,ires,iresdif,&
       Qxp,Qzp,Qxdif,Qzdif,Wavep,Wavedif,VQx,VQz,Vauxqxp,Vauxqxdif,&
       Iqp,Iqdif,nqcd,VQQ,BremL1D,BremS1D,GcollL1D,GcollS1D)
    !USE Common_Params
    !USE Common_Arrays
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Math_Constants
    use Phys_Constants
    implicit none
    integer :: nqx,nqz,ires,iresdif
    integer :: nqcd
    real(wp) :: Qxp,Qzp,Qxdif,Qzdif,Iqp,Iqdif,ILinit,ISinit,ITinit
    real(wp) :: Qp,Qp2
    real(wp) :: BremL0,BremS0,GcollL0,GcollS0
    real(wp) :: Aux
    real(wp), dimension(nqcd) :: BremL1D,BremS1D,GcollL1D,GcollS1D
    real(wp), dimension(nqcd) :: VQQ
    real(wp), dimension(nqx) :: VQx
    real(wp), dimension(nqz) :: VQz
    real(wp), dimension(nqx) :: Vauxqxp,Vauxqxdif
    character(LEN=1) :: Wavep,Wavedif

    Iqp= 0._wp
    Iqdif= 0._wp
    if (Qxp>0._wp) then
       if (ires>0 .and. ires<nqx) then
          call Aitp1d2(nqx,VQx,Vauxqxp,Qxp,Iqp,ires)
       else
          if (ires>=nqx) then
             Qp2= Qxp**2+Qzp**2
             Qp= sqrt(Qp2)
             if(Qp<=5.e-3_wp) Qp=5.e-3_wp
             call Locate(VQQ,nqcd,Qp,ires)
             call Aitp1d2(nqcd,VQQ,BremL1D,Qp,Aux,ires)
             BremL0= Aux
             call Aitp1d2(nqcd,VQQ,BremS1D,Qp,Aux,ires)
             BremS0= Aux
             call Aitp1d2(nqcd,VQQ,GcollL1D,Qp,Aux,ires)
             GcollL0= Aux
             call Aitp1d2(nqcd,VQQ,GcollS1D,Qp,Aux,ires)
             GcollS0= Aux
             call Iwave_Init(Qxp,Qzp,BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
             select case(Wavep)
             case("L")
                Iqp= ILinit
             case("S")
                Iqp= ISinit
             case("T")
                Iqp= ITinit
             case DEFAULT
             end select
          else
             Iqp= Vauxqxp(1)
          end if
       end if
       if(abs(Qzdif)<VQz(nqz)) then
          if (iresdif>0 .and. iresdif<nqx) then
             call Aitp1d2(nqx,VQx,Vauxqxdif,abs(Qxdif),Iqdif,iresdif)
          else
             if (iresdif>=nqx) then
                Qp2= Qxdif**2+Qzdif**2
                Qp= sqrt(Qp2)
                if(Qp<=5.e-3_wp) Qp=5.e-3_wp
                call Locate(VQQ,nqcd,Qp,ires)
                call Aitp1d2(nqcd,VQQ,BremL1D,Qp,Aux,ires)
                BremL0= Aux
                call Aitp1d2(nqcd,VQQ,BremS1D,Qp,Aux,ires)
                BremS0= Aux
                call Aitp1d2(nqcd,VQQ,GcollL1D,Qp,Aux,ires)
                GcollL0= Aux
                call Aitp1d2(nqcd,VQQ,GcollS1D,Qp,Aux,ires)
                GcollS0= Aux
                call Iwave_Init(abs(Qxdif),abs(Qzdif),BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
                select case(Wavedif)
                case("L")
                   Iqdif= ILinit
                case("S")
                   Iqdif= ISinit
                case("T")
                   Iqdif= ITinit
                case DEFAULT
                end select
             else
                Iqdif= Vauxqxdif(1)
             end if
          end if
       else
          Qp2= Qxdif**2+Qzdif**2
          Qp= sqrt(Qp2)
          if(Qp<=5.e-3_wp) Qp=5.e-3_wp
          call Locate(VQQ,nqcd,Qp,ires)
          call Aitp1d2(nqcd,VQQ,BremL1D,Qp,Aux,ires)
          BremL0= Aux
          call Aitp1d2(nqcd,VQQ,BremS1D,Qp,Aux,ires)
          BremS0= Aux
          call Aitp1d2(nqcd,VQQ,GcollL1D,Qp,Aux,ires)
          GcollL0= Aux
          call Aitp1d2(nqcd,VQQ,GcollS1D,Qp,Aux,ires)
          GcollS0= Aux
          call Iwave_Init(abs(Qxdif),abs(Qzdif),BremL0,BremS0,GcollL0,GcollS0,ILinit,ISinit,ITinit)
          select case(Wavedif)
          case("L")
             Iqdif= ILinit
          case("S")
             Iqdif= ISinit
          case("T")
             Iqdif= ITinit
          case DEFAULT
          end select
       end if
    else
       Iqp= 0._wp
       Iqdif= 0._wp
    end if

    return
  end subroutine Aux_Coef_Decay_z

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

  real(wp) function ZS(Qx,Qz)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    implicit none
    real(wp), intent(in) :: Qx,Qz
    real(wp) :: Q,Q2
    integer :: m

    m= 1   ! Auxiliary quantity, for the spatial profile.
    Q2= Qx**2+Qz**2
    Q= sqrt(Q2)
    ZS= Q*AA*sqrt(VRTeTs(m))/sqrt(1._wp+Q2/2._wp*VRTeTs(m)/VRNeNs(m))
    if (Qz < 0._wp) then
       ZS= -ZS
    else
    end if
    return
  end function ZS

  real(wp) function ZT(Qx,Qz)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    implicit none
    real(wp), intent(in) :: Qx,Qz
    real(wp) :: Q2

    Q2= Qx**2+Qz**2
    ZT= sqrt(1._wp+Q2/Ve2C2)
    if (Qz < 0._wp) then
       ZT= -ZT
    else
    end if
    return
  end function ZT

  subroutine ZBRAC2(FUNC,X1,X2,SUCCES) 
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    real(wp), intent(INOUT) :: X1,X2
    real(wp), parameter :: FACTOR= 1.6_wp
    real(wp) :: F1,F2 
    integer, parameter :: NTRY= 10
    integer :: J
    logical, intent(OUT) :: SUCCES 
    !      REAL(wp) :: X1,X2,FUNC,FACTOR
    !      INTEGER :: NTRY,J
    !      LOGICAL :: SUCCES 
    !      EXTERNAL FUNC
    !      PARAMETER (FACTOR=1.6_wp,NTRY=10)
    !        Given a function FUNC and an initial guessed range X1 to X2, the
    !        routine expands the range geometrically until a root is bracketed by
    !        the returned values X1 and X2 (in which case succes returns as .true.)
    !        or until the range becomes unacceptably large (in which case succes
    !        returns as .false.).
    !        This is a version of ZBRAC. It expands only the upper limit of the 
    !        range.
    interface 
       real(wp) function FUNC(X)
        use,intrinsic:: iso_fortran_env, only: wp=>real64
        real(wp), intent(in) :: X
       end function FUNC
    end interface
    if(X1.eq.X2)then
       open(98,FILE='Warning_Zbrac2.wt')
       write(98,*) 'You have to guess an initial range in zbrac'
       close(98)
       stop
    else
    end if
    F1=FUNC(X1) 
    F2=FUNC(X2) 
    SUCCES=.true. 
    do J=1,NTRY 
       if(F1*F2.lt.0._wp)return 
       X2=X2+FACTOR*(X2-X1) 
       F2=FUNC(X2) 
    end do
    SUCCES=.false. 
    return 
  end subroutine ZBRAC2

  real(wp) function RTBIS(FUNC,X1,X2,XACC) 
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    integer JMAX,J
    real(wp) X1,X2,XACC
    real(wp) F,FMID,DX,XMID
    !      REAL RTBIS,X1,X2,XACC,FUNC
    !      EXTERNAL FUNC
    parameter (JMAX=40)     ! Maximum allowed number of bisections.
    !         Using bisection, find the root of a function FUNC known to lie
    !         between x1 and x2. The root, returned as RTBIS, will be refined
    !         until its accuracy is +-XACC.
    interface 
       real(wp) function FUNC(X)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: X
       end function FUNC
    end interface
    FMID=FUNC(X2) 
    F=FUNC(X1) 
    if(F*FMID.ge.0._wp)then
       open(98,FILE='Warning_Rtbis.wt')
       write(98,*)'Root must be bracketed in rtbis.'
       close(98)
       stop
    else
    end if
    if(F.lt.0._wp)then         ! Orient the search so that f>0 lies at x+dx
       RTBIS=X1 
       DX=X2-X1 
    else 
       RTBIS=X2 
       DX=X1-X2 
    endif
    do J=1,JMAX           ! Bisection loop
       DX=DX*.5_wp
       XMID=RTBIS+DX 
       FMID=FUNC(XMID) 
       if(FMID.lt.0._wp)RTBIS=XMID 
       if(abs(DX).lt.XACC .or. FMID.eq.0._wp) return 
    end do
    open(98,FILE='Warning_Rtbis.wt')
    write(98,*)'Too many bisections in rtbis' 
    close(98)
    stop
  end function RTBIS

  real(wp) function Funcx_RcdLLS(Qxp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp), intent(in) :: Qxp
    real(wp) :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
    integer :: Asp,Aspp,S1,S2,S3 
    integer :: m

    Qx= Aux1_Rcd(1)
    Qz= Aux1_Rcd(2)
    Qzp= Aux1_Rcd(3)
    Asp= Aux2_Rcd(1)
    Aspp= Aux2_Rcd(2)
    S1= Aux2_Rcd(3)
    S2= Aux2_Rcd(4)
    Qxdif= Qx-S1*Qxp
    Qzdif= Qz-S2*Qzp
    S3= sign(1._wp,Qzdif)
    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    Funcx_RcdLLS= ZL(Qx,abs(Qz))-Asp*S2*ZL(Qxp,abs(Qzp)) &
         -Aspp*S3*ZS(Qxdif,abs(Qzdif))
    return
  end function Funcx_RcdLLS

  real(wp) function Funcx_RcdLLT(Qxp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp), intent(in) :: Qxp
    real(wp) :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
    integer :: Asp,Aspp,S1,S2,S3 
    integer :: m

    Qx= Aux1_Rcd(1)
    Qz= Aux1_Rcd(2)
    Qzp= Aux1_Rcd(3)
    Asp= Aux2_Rcd(1)
    Aspp= Aux2_Rcd(2)
    S1= Aux2_Rcd(3)
    S2= Aux2_Rcd(4)
    Qxdif= Qx-S1*Qxp
    Qzdif= Qz-S2*Qzp
    S3= sign(1._wp,Qzdif)
    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    Funcx_RcdLLT= ZL(Qx,abs(Qz))-Asp*S2*ZL(Qxp,abs(Qzp)) &
         -Aspp*S3*ZT(Qxdif,abs(Qzdif))
    return
  end function Funcx_RcdLLT

  real(wp) function Funcx_RcdLST(Qxp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp), intent(in) :: Qxp
    real(wp) :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
    integer :: Asp,Aspp,S1,S2,S3 
    integer :: m

    Qx= Aux1_Rcd(1)
    Qz= Aux1_Rcd(2)
    Qzp= Aux1_Rcd(3)
    Asp= Aux2_Rcd(1)
    Aspp= Aux2_Rcd(2)
    S1= Aux2_Rcd(3)
    S2= Aux2_Rcd(4)
    Qxdif= Qx-S1*Qxp
    Qzdif= Qz-S2*Qzp
    S3= sign(1._wp,Qzdif)
    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    Funcx_RcdLST= ZL(Qx,abs(Qz))-Asp*S2*ZS(Qxp,abs(Qzp)) &
         -Aspp*S3*ZT(Qxdif,abs(Qzdif))
    return
  end function Funcx_RcdLST

  real(wp) function Funcx_RcdLTT(Qxp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp), intent(in) :: Qxp
    real(wp) :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
    integer :: Asp,Aspp,S1,S2,S3 
    integer :: m

    Qx= Aux1_Rcd(1)
    Qz= Aux1_Rcd(2)
    Qzp= Aux1_Rcd(3)
    Asp= Aux2_Rcd(1)
    Aspp= Aux2_Rcd(2)
    S1= Aux2_Rcd(3)
    S2= Aux2_Rcd(4)
    Qxdif= Qx-S1*Qxp
    Qzdif= Qz-S2*Qzp
    S3= sign(1._wp,Qzdif)
    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    Funcx_RcdLTT= ZL(Qx,abs(Qz))-Asp*S2*ZT(Qxp,abs(Qzp)) &
         -Aspp*S3*ZT(Qxdif,abs(Qzdif))
    return
  end function Funcx_RcdLTT

  real(wp) function Funcx_RcdSLL(Qxp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp), intent(in) :: Qxp
    real(wp) :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
    integer :: Asp,Aspp,S1,S2,S3 
    integer :: m

    Qx= Aux1_Rcd(1)
    Qz= Aux1_Rcd(2)
    Qzp= Aux1_Rcd(3)
    Asp= Aux2_Rcd(1)
    Aspp= Aux2_Rcd(2)
    S1= Aux2_Rcd(3)
    S2= Aux2_Rcd(4)
    Qxdif= Qx-S1*Qxp
    Qzdif= Qz-S2*Qzp
    S3= sign(1._wp,Qzdif)
    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    Funcx_RcdSLL= ZS(Qx,abs(Qz))-Asp*S2*ZL(Qxp,abs(Qzp)) &
         -Aspp*S3*ZL(Qxdif,abs(Qzdif))
    return
  end function Funcx_RcdSLL

  real(wp) function Funcx_RcdSLT(Qxp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp), intent(in) :: Qxp
    real(wp) :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
    integer :: Asp,Aspp,S1,S2,S3 
    integer :: m

    Qx= Aux1_Rcd(1)
    Qz= Aux1_Rcd(2)
    Qzp= Aux1_Rcd(3)
    Asp= Aux2_Rcd(1)
    Aspp= Aux2_Rcd(2)
    S1= Aux2_Rcd(3)
    S2= Aux2_Rcd(4)
    Qxdif= Qx-S1*Qxp
    Qzdif= Qz-S2*Qzp
    S3= sign(1._wp,Qzdif)
    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    Funcx_RcdSLT= ZS(Qx,abs(Qz))-Asp*S2*ZL(Qxp,abs(Qzp)) &
         -Aspp*S3*ZT(Qxdif,abs(Qzdif))
    return
  end function Funcx_RcdSLT

  real(wp) function Funcx_RcdTLL(Qxp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp), intent(in) :: Qxp
    real(wp) :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
    integer :: Asp,Aspp,S1,S2,S3 
    integer :: m

    Qx= Aux1_Rcd(1)
    Qz= Aux1_Rcd(2)
    Qzp= Aux1_Rcd(3)
    Asp= Aux2_Rcd(1)
    Aspp= Aux2_Rcd(2)
    S1= Aux2_Rcd(3)
    S2= Aux2_Rcd(4)
    Qxdif= Qx-S1*Qxp
    Qzdif= Qz-S2*Qzp
    S3= sign(1._wp,Qzdif)
    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    Funcx_RcdTLL= ZT(Qx,abs(Qz))-Asp*S2*ZL(Qxp,abs(Qzp)) &
         -Aspp*S3*ZL(Qxdif,abs(Qzdif))
    return
  end function Funcx_RcdTLL

  real(wp) function Funcx_RcdTLS(Qxp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp), intent(in) :: Qxp
    real(wp) :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
    integer :: Asp,Aspp,S1,S2,S3 
    integer :: m

    Qx= Aux1_Rcd(1)
    Qz= Aux1_Rcd(2)
    Qzp= Aux1_Rcd(3)
    Asp= Aux2_Rcd(1)
    Aspp= Aux2_Rcd(2)
    S1= Aux2_Rcd(3)
    S2= Aux2_Rcd(4)
    Qxdif= Qx-S1*Qxp
    Qzdif= Qz-S2*Qzp
    S3= sign(1._wp,Qzdif)
    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    Funcx_RcdTLS= ZT(Qx,abs(Qz))-Asp*S2*ZL(Qxp,abs(Qzp)) &
         -Aspp*S3*ZS(Qxdif,abs(Qzdif))
    return
  end function Funcx_RcdTLS

  real(wp) function Funcx_RcdTTL(Qxp)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    real(wp), intent(in) :: Qxp
    real(wp) :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
    integer :: Asp,Aspp,S1,S2,S3 
    integer :: m

    Qx= Aux1_Rcd(1)
    Qz= Aux1_Rcd(2)
    Qzp= Aux1_Rcd(3)
    Asp= Aux2_Rcd(1)
    Aspp= Aux2_Rcd(2)
    S1= Aux2_Rcd(3)
    S2= Aux2_Rcd(4)
    Qxdif= Qx-S1*Qxp
    Qzdif= Qz-S2*Qzp
    S3= sign(1._wp,Qzdif)
    m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
    Beta= VRTeTs(m)/VRNeNs(m)

    Funcx_RcdTTL= ZT(Qx,abs(Qz))-Asp*S2*ZT(Qxp,abs(Qzp)) &
         -Aspp*S3*ZL(Qxdif,abs(Qzdif))
    return
  end function Funcx_RcdTTL

  subroutine Fnorm(Anorm)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    integer :: l,m
    real(wp) :: Anorm
    real(wp), dimension(nux) :: Vintux
    real(wp), dimension(nuz) :: Vintuz

    do l= 1,nux
       forall(m=1:nuz)
          Vintuz(m)= Fe(l,m)
       end forall
       call Simpson(VUz,Vintuz,nuz,Anorm)
       Vintux(l)= Anorm
    end do
    call Simpson(VUx,Vintux,nux,Anorm)
    Anorm= 2._wp*Anorm  ! Symmetry in Ux.
    if (RenormFe=="Yes") then 
       Fe= Fe/Anorm
    else
    end if
    return
  end subroutine Fnorm

  subroutine Energy(Epart,Ewave,EppEw)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    integer :: i,k
    real(wp) :: Ewave,Epart,EppEw
    real(wp) :: Qx2,Uz2,Aux
    real(wp), dimension(nux) :: Vintux
    real(wp), dimension(nuz) :: Vintuz
    real(wp), dimension(nqx) :: Vintqx
    real(wp), dimension(nqz) :: Vintqz

    do k= 1,nuz
       Uz2= (VUz(k))**2
       forall(i=1:nux)
          Vintux(i)= (Uz2+(VUx(i))**2)*Fe(i,k)
       end forall
       call Simpson(VUx,Vintux,nux,Aux)
       Vintuz(k)= Aux
    end do
    call Simpson(VUz,Vintuz,nuz,Aux)
    Epart= Aux/2._wp

    do i= 1,nqx
       Qx2= (VQx(i))**2
       forall(k=1:nqz)
          Vintqz(k)= ILp(i,k)+ILm(i,k) &
               + (sqrt(Qx2+(VQz(k))**2))**3*AA/2._wp*(ISp(i,k)+ISm(i,k)) &
               + ITp(i,k)+ITm(i,k)
       end forall
       call Simpson(VQz,Vintqz,nqz,Aux)
       Vintqx(i)= Aux
    end do
    call Simpson(VQx,Vintqx,nqx,Aux)
    Ewave= Aux

    EppEw= Epart+Ewave

    return
  end subroutine Energy

  subroutine Energy2(EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    integer :: i,k
    real(wp) :: EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
    real(wp) :: Qx2,Aux
    real(wp) :: Q,Dq,Qxaux,Qzaux
    real(wp), dimension(nqx) :: Vintqx
    real(wp), dimension(nqz) :: Vintqz
    real(wp), dimension(nph) :: VQ,Vintq,VZtq,VintZq

    do i= 1,nqx
       forall(k=1:nqz)
          Vintqz(k)= ILp(i,k)+ILm(i,k)
       end forall
       call Simpson(VQz,Vintqz,nqz,Aux)
       Vintqx(i)= Aux
    end do
    call Simpson(VQx,Vintqx,nqx,Aux)
    EwaveL= Aux

    do i= 1,nqx
       Qx2= (VQx(i))**2
       forall(k=1:nqz)
          Vintqz(k)= (sqrt(Qx2+(VQz(k))**2))**3*AA/2._wp*(ISp(i,k)+ISm(i,k))
       end forall
       call Simpson(VQz,Vintqz,nqz,Aux)
       Vintqx(i)= Aux
    end do
    call Simpson(VQx,Vintqx,nqx,Aux)
    EwaveS= Aux

    do i= 1,nqx
       forall(k=1:nqz)
          Vintqz(k)= ITp(i,k)+ITm(i,k)
       end forall
       call Simpson(VQz,Vintqz,nqz,Aux)
       Vintqx(i)= Aux
    end do
    call Simpson(VQx,Vintqx,nqx,Aux)
    EwaveT= Aux

    Dq= (Qzf-Qmin)/(nph-1)
    do i= 1,nph
       Q= Qmin+(i-1)*Dq
       VQ(i)= Q
       VZtq(i)= sqrt(1._wp+Q**2/Ve2C2)
    end do
    do i= 1,nph
       Q= VQ(i)
       do k= 1,nph
          Qxaux= Q*sin(VPhip(k))
          Qzaux= Q*cos(VPhip(k))
          call Aitp2d(nqx,nqz,VQx,VQz,ITp,Qxaux,Qzaux,Aux)
          Vintq(k)= Aux
          call Aitp2d(nqx,nqz,VQx,VQz,ITm,Qxaux,Qzaux,Aux)
          Vintq(k)= Vintq(k)+Aux
       end do
       call Simpson(VPhip,Vintq,nph,Aux)
       VintZq(i)= Aux*VQ(i)
    end do
    EwaveTF= 0._wp
    EwaveTH= 0._wp
    EwaveT3= 0._wp
    do i= 1,nph-1
       if (VZtq(i)<1.5_wp) then
          EwaveTF= EwaveTF+Dq*(VintZq(i)+VintZq(i+1))/2._wp   
       else
          if (VZtq(i)<2.5_wp) then
             EwaveTH= EwaveTH+Dq*(VintZq(i)+VintZq(i+1))/2._wp   
          else
             if (VZtq(i)<3.5_wp) then
                EwaveT3= EwaveT3+Dq*(VintZq(i)+VintZq(i+1))/2._wp   
             else
             end if
          end if
       end if
    end do

    return
  end subroutine Energy2

  subroutine Output_Coef(CoefChoice,CoefA,CoefB)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    integer :: i,k
    integer, parameter :: Step=1
    real(wp) :: Qx,Qz
    real(wp), dimension(nqx,nqz) :: CoefA,CoefB
    character(LEN=5) :: CoefChoice

    select case(CoefChoice)
    case("Lwave")
       open(1,FILE='Lwave.wt')
    case("Swave")
       open(1,FILE='Swave.wt')
    case("Twave")
       open(1,FILE='Twave.wt')
    case DEFAULT
    end select

    do i= 1,nqx
       Qx= VQx(nqx+1-i)
       write(1,*)' '
       do k= 1,nqz
          Qz= VQz(nqz+1-k)
          write(1,*) -Qx,-Qz,CoefA(nqx+1-i,nqz+1-k),CoefB(nqx+1-i,nqz+1-k)
       end do
       do k= 1,nqz
          Qz= VQz(k)
          write(1,*) -Qx,Qz,CoefA(nqx+1-i,k),CoefB(nqx+1-i,k)
       end do
    end do
    do i= 1,nqx
       Qx= VQx(i)
       write(1,*)' '
       do k= 1,nqz
          Qz= VQz(nqz+1-k)
          write(1,*) Qx,-Qz,CoefA(i,nqz+1-k),CoefB(i,nqz+1-k)
       end do
       do k= 1,nqz
          Qz= VQz(k)
          write(1,*) Qx,Qz,CoefA(i,k),CoefB(i,k)
       end do
    end do
    close(1)

    return
  end subroutine Output_Coef

  subroutine Output(WriteChoice)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    integer :: i,k
    integer, parameter :: Step=1, np= 128 !odd
    real(wp) :: Aux
    real(wp) :: Qx,Qz,Q,Muq
    real(wp) :: U,Phi,DPhi,Du
    real(wp) :: UxAux,UzAux
    real(wp), dimension(nqx) :: Vintqx
    real(wp), dimension(nqz) :: Iqzm,Iqzp,Iqzm2,Iqzp2
    real(wp), dimension(nux) :: Vintux
    real(wp), dimension(nuz) :: Feuz
    real(wp), dimension(np) :: VU,VPhi
    real(wp), dimension(np) :: Feu1,Feu2
    real(wp), dimension(np) :: Feu0,Feu00
    real(wp), dimension(np) :: Vintu0,Vintu00,Vintu1,Vintu2
    character(LEN=3) :: WriteChoice

    select case(WriteChoice)

    case("Feu")
       open(1,FILE="Feu.wt")
       open(2,FILE="Feu0.wt")
       Du= (Ulim-EpsMin)/(np-1)
       do i=1,np
          U= EpsMin+(i-1)*Du
          VU(i)= U
       end do
       DPhi= (Pi-EpsMin)/(np-1)
       do i=1,np
          Phi= EpsMin+(i-1)*DPhi
          VPhi(i)= Phi
       end do
       do i= 1,np
          U= VU(i)
          do k= 1,np
             Phi= VPhi(k)
             UxAux= U*sin(Phi)
             UzAux= U*cos(Phi)
             call Aitp2d(nux,nuz,VUx,VUz,Fe0,Uxaux,Uzaux,Aux)
             Vintu0(k)= Aux
             Vintu00(k)= U*Aux
             call Aitp2d(nux,nuz,VUx,VUz,Fe,Uxaux,Uzaux,Aux)
             Vintu1(k)= Aux
             Vintu2(k)= U*Aux
          end do
          call Simpson(VPhi,Vintu0,np,Aux)
          Feu0(i)= Aux
          call Simpson(VPhi,Vintu00,np,Aux)
          Feu00(i)= Aux
          call Simpson(VPhi,Vintu1,np,Aux)
          Feu1(i)= Aux
          call Simpson(VPhi,Vintu2,np,Aux)
          Feu2(i)= Aux
       end do

       do k= 1,np!,Step
          write(1,*) VU(k),Feu1(k),Feu2(k)
          write(2,*) VU(k),Feu0(k),Feu00(k)
       end do
       close(1)
       close(2)
       
    case("Fi ")
       open(1,FILE='Fi.wt')
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fi(nux+1-i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fi(nux+1-i,k)
             end if
             write(1,*) -VUx(nux+1-i),VUz(k),Aux
          end do
       end do
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fi(i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fi(i,k)
             end if
             write(1,*) VUx(i),VUz(k),Aux
          end do
       end do
       close(1)
       open(1,FILE= "Fi1D.wt")
       do k= 1,nuz
          if(abs(Fi(1,k))<=EpsMin) then
             Aux= EpsMin
          else
             Aux= Fi(1,k)
          end if
          write(1,*) VUz(k),Aux
       end do
       close(1)     

    case("Fe0")
       open(1,FILE='Fe0.wt')
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe0(nux+1-i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe0(nux+1-i,k)
             end if
             write(1,*) -VUx(nux+1-i),VUz(k),Aux
          end do
       end do
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe0(i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe0(i,k)
             end if
             write(1,*) VUx(i),VUz(k),Aux
          end do
       end do
       close(1)
       open(1,FILE= "Fe01D.wt")
       do k= 1,nuz
          if(abs(Fe0(1,k))<=EpsMin) then
             Aux= EpsMin
          else
             Aux= Fe0(1,k)
          end if
          write(1,*) VUz(k),Aux
       end do
       close(1)
       open(1,FILE='Fe1D0.wt')
       do k= 1,nuz
          do i= 1,nux
             Vintux(i)= Fe(i,k)
          end do
          call Simpson(VUx,Vintux,nux,Aux)
          Feuz(k)= Aux
       end do
       do k= 1,nuz,Step
          write(1,*) VUz(k),Feuz(k)
       end do
       close(1)

    case("Fe ")
       open(1,FILE='Fe.wt')
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe(nux+1-i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe(nux+1-i,k)
             end if
             write(1,*) -VUx(nux+1-i),VUz(k),Aux
          end do
       end do
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe(i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe(i,k)
             end if
             write(1,*) VUx(i),VUz(k),Aux
          end do
       end do
       close(1)
       open(1,FILE= "Fe1D.wt")
       do k= 1,nuz
          if(abs(Fe(1,k))<=EpsMin) then
             Aux= EpsMin
          else
             Aux= Fe(1,k)
          end if
          write(1,*) VUz(k),Aux
       end do
       close(1)

    case("Ai ")
       open(1,FILE='Ai.wt')
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             write(1,*) -VUx(nux+1-i),VUz(k),Ax(nux+1-i,k),Az(nux+1-i,k)
          end do
       end do
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             write(1,*) VUx(i),VUz(k),Ax(i,k),Az(i,k)
          end do
       end do
       close(1)

    case("Dij")
       open(1,FILE='Dij.wt')
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             write(1,*) -VUx(nux+1-i),VUz(k),Dxx(nux+1-i,k),Dxz(nux+1-i,k),&
                  Dzz(nux+1-i,k)
          end do
       end do
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             write(1,*) VUx(i),VUz(k),Dxx(i,k),Dxz(i,k),Dzz(i,k)
          end do
       end do
       close(1)

    case("IL ")
       open(1,FILE='IL.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,ILm(nqx+1-i,nqz+1-k),ZL(-Qx,-Qz)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,ILp(nqx+1-i,k),ZL(-Qx,Qz)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,ILm(i,nqz+1-k),ZL(Qx,-Qz)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,ILp(i,k),ZL(Qx,Qz)
          end do
       end do
       close(1)
       open(1,FILE= "IL1D.wt")
       do k= 1,nqz
          write(1,*) VQz(k),ILp(1,k)
       end do
       close(1)

    case("IS ")
       open(1,FILE='IS.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) -Qx,-Qz,Muq*ISm(nqx+1-i,nqz+1-k),ZS(-Qx,-Qz)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) -Qx,Qz,Muq*ISp(nqx+1-i,k),ZS(-Qx,Qz)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Qx,-Qz,Muq*ISm(i,nqz+1-k),ZS(Qx,-Qz)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Qx,Qz,Muq*ISp(i,k),ZS(Qx,Qz)
          end do
       end do
       close(1)
       open(1,FILE= "IS1D.wt")
       do k= 1,nqz
          Qz= VQz(k)
          Q= sqrt(VQx(1)**2+Qz**2)
          Muq= Q**3*AA/2._wp
          write(1,*) VQz(k),Muq*ISp(1,k)
       end do
       close(1)

    case("IT ")
       open(1,FILE='IT.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,ITm(nqx+1-i,nqz+1-k),ZT(-Qx,-Qz)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,ITp(nqx+1-i,k),ZT(-Qx,Qz)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,ITm(i,nqz+1-k),ZT(Qx,-Qz)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,ITp(i,k),ZT(Qx,Qz)
          end do
       end do
       close(1)

    case("IL0")
       open(1,FILE='IL0.wt')
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),ILm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),VQz(k),ILp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) VQx(i),-VQz(nqz+1-k),ILm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) VQx(i),VQz(k),ILp(i,k)
          end do
       end do
       close(1)
       open(1,FILE= "IL01D.wt")
       do k= 1,nqz
          write(1,*) VQz(k),ILp(1,k)
       end do
       close(1)
       open(1,FILE='IL1D0.wt')
       do k= 1,nqz
          do i= 1,nqx
             Vintqx(i)= ILm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm(k)= Aux
          do i= 1,nqx
             Vintqx(i)= ILp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp(k)= Aux
       end do
       do k= 1,nqz,Step
          write(1,*) -VQz(nqz+1-k),Iqzm(nqz+1-k)
       end do
       do k= 1,nqz,Step
          write(1,*) VQz(k),Iqzp(k)
       end do
       close(1)

    case("IS0")
       open(1,FILE='IS0.wt')
       do i= 1,nqx,Step
          Qx= -VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Qx,-Qz,Muq*ISm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Qx,Qz,Muq*ISp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Qx,-Qz,Muq*ISm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Qx,Qz,Muq*ISp(i,k)
          end do
       end do
       close(1)
       open(1,FILE= "IS01D.wt")
       do k= 1,nqz
          Qz= VQz(k)
          Q= sqrt(VQx(1)**2+Qz**2)
          Muq= Q**3*AA/2._wp
          write(1,*) VQz(k),Muq*ISp(1,k)
       end do
       close(1)
       open(1,FILE='IS1D0.wt')
       do k= 1,nqz
          Qz= VQz(k)
          do i= 1,nqx
             Vintqx(i)= ISm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm(k)= Aux
          do i= 1,nqx
             Qx= VQx(i)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             Vintqx(i)= Muq*ISm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm2(k)= Aux
          do i= 1,nqx
             Vintqx(i)= ISp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp(k)= Aux
          do i= 1,nqx
             Qx= VQx(i)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             Vintqx(i)= Muq*ISp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp2(k)= Aux
       end do
       do k= 1,nqz,Step
          write(1,*) -VQz(nqz+1-k),Iqzm2(nqz+1-k),Iqzm(nqz+1-k)
       end do
       do k= 1,nqz,Step
          write(1,*) VQz(k),Iqzp2(k),Iqzp(k)
       end do
       close(1)

    case("IT0")
       open(1,FILE='IT0.wt')
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),ITm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) -VQx(nqx+1-i),VQz(k),ITp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) VQx(i),-VQz(nqz+1-k),ITm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) VQx(i),VQz(k),ITp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='IT1D0.wt')
       do k= 1,nqz
          do i= 1,nqx
             Vintqx(i)= ITm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm(k)= Aux
          do i= 1,nqx
             Vintqx(i)= ITp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp(k)= Aux
       end do
       do k= 1,nqz,Step
          write(1,*) -VQz(nqz+1-k),Iqzm(nqz+1-k)
       end do
       do k= 1,nqz,Step
          write(1,*) VQz(k),Iqzp(k)
       end do
       close(1)

       ! CASE("Fe1")
       !  OPEN(1,FILE='Fe1D.wt')
       !  DO k= 1,nuz
       !   DO i= 1,nux
       !    Vintux(i)= Fe(i,k)
       !   END DO
       !   CALL Simpson(VUx,Vintux,nux,Aux)
       !   Feuz(k)= Aux
       !  END DO
       !  DO k= 1,nuz,Step
       !   WRITE(1,*) VUz(k),Feuz(k)
       !  END DO
       !  CLOSE(1)

       ! CASE("IL1")
       !  OPEN(1,FILE='IL1D.wt')
       !  DO k= 1,nqz
       !   DO i= 1,nqx
       !    Vintqx(i)= ILm(i,k)
       !   END DO
       !   CALL Simpson(VQx,Vintqx,nqx,Aux)
       !   Iqzm(k)= Aux
       !   DO i= 1,nqx
       !    Vintqx(i)= ILp(i,k)
       !   END DO
       !   CALL Simpson(VQx,Vintqx,nqx,Aux)
       !   Iqzp(k)= Aux
       !  END DO
       !  DO k= 1,nqz,Step
       !   WRITE(1,*) -VQz(nqz+1-k),Iqzm(nqz+1-k)
       !  END DO
       !  DO k= 1,nqz,Step
       !   WRITE(1,*) VQz(k),Iqzp(k)
       !  END DO
       !  CLOSE(1)

       ! CASE("IS1")
       !  OPEN(1,FILE='IS1D.wt')
       !  DO k= 1,nqz
       !   Qz= VQz(k)
       !   DO i= 1,nqx
       !    Vintqx(i)= ISm(i,k)
       !   END DO
       !   CALL Simpson(VQx,Vintqx,nqx,Aux)
       !   Iqzm(k)= Aux
       !   DO i= 1,nqx
       !    Qx= VQx(i)
       !    Q= SQRT(Qx**2+Qz**2)
       !    Muq= Q**3*AA/2._wp
       !    Vintqx(i)= Muq*ISm(i,k)
       !   END DO
       !   CALL Simpson(VQx,Vintqx,nqx,Aux)
       !   Iqzm2(k)= Aux
       !   DO i= 1,nqx
       !    Vintqx(i)= ISp(i,k)
       !   END DO
       !   CALL Simpson(VQx,Vintqx,nqx,Aux)
       !   Iqzp(k)= Aux
       !   DO i= 1,nqx
       !    Qx= VQx(i)
       !    Q= SQRT(Qx**2+Qz**2)
       !    Muq= Q**3*AA/2._wp
       !    Vintqx(i)= Muq*ISp(i,k)
       !   END DO
       !   CALL Simpson(VQx,Vintqx,nqx,Aux)
       !   Iqzp2(k)= Aux
       !  END DO
       !  DO k= 1,nqz,Step
       !   WRITE(1,*) -VQz(nqz+1-k),Iqzm2(nqz+1-k),Iqzm(nqz+1-k)
       !  END DO
       !  DO k= 1,nqz,Step
       !   WRITE(1,*) VQz(k),Iqzp2(k),Iqzp(k)
       !  END DO
       !  CLOSE(1)

    case("IT1")
       open(1,FILE='IT1D.wt')
       do k= 1,nqz
          do i= 1,nqx
             Vintqx(i)= ITm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm(k)= Aux
          do i= 1,nqx
             Vintqx(i)= ITp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp(k)= Aux
       end do
       do k= 1,nqz,Step
          write(1,*) -VQz(nqz+1-k),Iqzm(nqz+1-k)
       end do
       do k= 1,nqz,Step
          write(1,*) VQz(k),Iqzp(k)
       end do
       close(1)

    case("GcL")
       open(1,FILE='GcollL.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,GcollLm(nqx+1-i,nqz+1-k)! ,&
                  ! BremLm(nqx+1-i,nqz+1-k)/GcollLm(nqx+1-i,nqz+1-k),&
                  ! GcollLm(nqx+1-i,nqz+1-k)/GqlLm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,GcollLp(nqx+1-i,k)! ,&
                  ! -BremLp(nqx+1-i,k)/GcollLp(nqx+1-i,k),&
                  ! GcollLp(nqx+1-i,k)/GqlLp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,GcollLm(i,nqz+1-k)! ,&
                  ! BremLm(i,nqz+1-k)/GcollLm(i,nqz+1-k),&
                  ! GcollLm(i,nqz+1-k)/GqlLm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,GcollLp(i,k)! , &
                  ! -BremLp(i,k)/GcollLp(i,k), &
                  ! GcollLp(i,k)/GqlLp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='GcollL1D.wt')
       do i= 1,nqcd
          write(1,*) VQQ(i),GcollL1D(i),-BremL1D(i)/GcollL1D(i)
       end do
       close(1)
       ! open(1,FILE='TGcollL1D.wt')
       ! do i= 1,nqcd
       !    write(1,*) VQQ(i),GcollL1D(i),GcollLp1D(i)-GcollLm1D(i)
       ! end do
       ! close(1)

    case("GqL")
       open(1,FILE='GqlL.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,GqlLm(nqx+1-i,nqz+1-k),-GqlLm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,GqlLp(nqx+1-i,k),-GqlLp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,GqlLm(i,nqz+1-k),-GqlLm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,GqlLp(i,k),-GqlLp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='GqlL1D.wt')
       do i= 1,nqcd
          write(1,*) VQQ(i),GqlL1D(i)
       end do
       close(1)

    case("GcS")
       open(1,FILE='GcollS.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,GcollSm(nqx+1-i,nqz+1-k),-GcollSm(nqx+1-i,nqz+1-k),&
                  GcollSm(nqx+1-i,nqz+1-k)/GqlSm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,GcollSp(nqx+1-i,k),-GcollSp(nqx+1-i,k),&
                  GcollSp(nqx+1-i,k)/GqlSp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,GcollSm(i,nqz+1-k),-GcollSm(i,nqz+1-k),&
                  GcollSm(i,nqz+1-k)/GqlSm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,GcollSp(i,k),-GcollSp(i,k),&
                  GcollSp(i,k)/GqlSp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='GcollS1D.wt')
       do i= 1,nqcd
          write(1,*) VQQ(i),GcollS1D(i)
       end do
       close(1)

    case("GqS")
       open(1,FILE='GqlS.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,GqlSm(nqx+1-i,nqz+1-k),-GqlSm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,GqlSp(nqx+1-i,k),-GqlSp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,GqlSm(i,nqz+1-k),-GqlSm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,GqlSp(i,k),-GqlSp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='GqlS1D.wt')
       do i= 1,nqcd
          write(1,*) VQQ(i),GqlS1D(i)
       end do
       close(1)

       !400 FORMAT(1x,5(e13.5e3,1x))

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

    case("BrS")
       open(1,FILE='BremS.wt')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) -Qx,-Qz,BremSm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) -Qx,Qz,BremSp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) Qx,-Qz,BremSm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) Qx,Qz,BremSp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='BremS1D.wt')
       do i= 1,nqcd
          write(1,*) VQQ(i),BremS1D(i)
       end do
       close(1)

    case DEFAULT

    end select

    return

  end subroutine Output

  subroutine Output2(WriteChoice)
    ! Subroutine generating output in the format suitable for MatLab (Nov/2013)
    ! Nov. 28, 2013
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    implicit none
    integer :: i,k
    integer, parameter :: Step=1
    real(wp) :: Aux
    real(wp) :: Qx,Qz,Q,Muq
    real(wp) :: Dq,Qxaux,Qzaux
    real(wp), dimension(nqx) :: Vintqx
    real(wp), dimension(nqz) :: Iqzm,Iqzp,Iqzm2,Iqzp2
    real(wp), dimension(nux) :: Vintux
    real(wp), dimension(nuz) :: Feuz
    real(wp), dimension(nph) :: VQ,Vintq
    character(LEN=3) :: WriteChoice

    select case(WriteChoice)

    case("Ui ")
       open(1,FILE='Ux')
       do i= 1,nux,Step
          write(1,*) -VUx(nux+1-i)
       end do
       do i= 1,nux,Step
          write(1,*) VUx(i)
       end do
       close(1)
       open(1,FILE='Uz')
       do i= 1,nuz,Step
          write(1,*) VUz(i)
       end do
       close(1)

    case("Qi ")
       open(1,FILE='Qx')
       do i= 1,nqx,Step
          write(1,*) -VQx(nqx+1-i)
       end do
       do i= 1,nqx,Step
          write(1,*) VQx(i)
       end do
       close(1)
       open(1,FILE='Qz')
       do i= 1,nqz,Step
          write(1,*) -VQz(nqz+1-i)
       end do
       do i= 1,nqz,Step
          write(1,*) VQz(i)
       end do
       close(1)

    case("Q  ")
       open(1,FILE='Q')
       Dq= (Qzf-Qmin)/(nph-1)
       do i= 1,nph
          Q= Qmin+(i-1)*Dq
          write(1,*) Q
       end do
       close(1)
       open(1,FILE='ZTQ')
       do i= 1,nph
          Q= Qmin+(i-1)*Dq
          write(1,*) sqrt(1._wp+Q**2/Ve2C2)
       end do
       close(1)

    case("Fe0")
       open(1,FILE='Fe02D')
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe0(nux+1-i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe0(nux+1-i,k)
             end if
             write(1,*) Aux
          end do
       end do
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe0(i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe0(i,k)
             end if
             write(1,*) Aux
          end do
       end do
       close(1)
       open(1,FILE='Fe01D')
       do k= 1,nuz
          do i= 1,nux
             Vintux(i)= Fe(i,k)
          end do
          call Simpson(VUx,Vintux,nux,Aux)
          Feuz(k)= Aux
       end do
       do k= 1,nuz,Step
          write(1,*) Feuz(k)
       end do
       close(1)

    case("Fe ")
       open(1,FILE='Fe2D')
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe(nux+1-i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe(nux+1-i,k)
             end if
             write(1,*) Aux
          end do
       end do
       do i= 1,nux,Step
          write(1,*)' '
          do k= 1,nuz,Step
             if(abs(Fe(i,k))<=EpsMin) then
                Aux= EpsMin
             else
                Aux= Fe(i,k)
             end if
             write(1,*) Aux
          end do
       end do
       close(1)

    case("IL ")
       open(1,FILE='IL2D')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) ILm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) ILp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) ILm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) ILp(i,k)
          end do
       end do
       close(1)

    case("IS ")
       open(1,FILE='IS2D')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Muq*ISm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Muq*ISp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Muq*ISm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Muq*ISp(i,k)
          end do
       end do
       close(1)

    case("IT ")
       open(1,FILE='IT2D')
       do i= 1,nqx,Step
          Qx= VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) ITm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) ITp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             write(1,*) ITm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             write(1,*) ITp(i,k)
          end do
       end do
       close(1)

    case("IL0")
       open(1,FILE='IL02D')
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) ILm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) ILp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) ILm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) ILp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='IL01D')
       do k= 1,nqz
          do i= 1,nqx
             Vintqx(i)= ILm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm(k)= Aux
          do i= 1,nqx
             Vintqx(i)= ILp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp(k)= Aux
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzm(nqz+1-k)
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzp(k)
       end do
       close(1)

    case("IS0")
       open(1,FILE='IS02D')
       do i= 1,nqx,Step
          Qx= -VQx(nqx+1-i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Muq*ISm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Muq*ISp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          Qx= VQx(i)
          write(1,*)' '
          do k= 1,nqz,Step
             Qz= VQz(nqz+1-k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Muq*ISm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             Qz= VQz(k)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             write(1,*) Muq*ISp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='IS01D')
       do k= 1,nqz
          Qz= VQz(k)
          do i= 1,nqx
             Qx= VQx(i)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             Vintqx(i)= Muq*ISm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm2(k)= Aux
          do i= 1,nqx
             Qx= VQx(i)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             Vintqx(i)= Muq*ISp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp2(k)= Aux
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzm2(nqz+1-k)
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzp2(k)
       end do
       close(1)

    case("IT0")
       open(1,FILE='IT02D')
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) ITm(nqx+1-i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) ITp(nqx+1-i,k)
          end do
       end do
       do i= 1,nqx,Step
          write(1,*)' '
          do k= 1,nqz,Step
             write(1,*) ITm(i,nqz+1-k)
          end do
          do k= 1,nqz,Step
             write(1,*) ITp(i,k)
          end do
       end do
       close(1)
       open(1,FILE='IT01D')
       do k= 1,nqz
          do i= 1,nqx
             Vintqx(i)= ITm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm(k)= Aux
          do i= 1,nqx
             Vintqx(i)= ITp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp(k)= Aux
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzm(nqz+1-k)
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzp(k)
       end do
       close(1)

       ! CASE("Fe1")
       !  OPEN(1,FILE='Fe1D')
       !  DO k= 1,nuz
       !   DO i= 1,nux
       !    Vintux(i)= Fe(i,k)
       !   END DO
       !   CALL Simpson(VUx,Vintux,nux,Aux)
       !   Feuz(k)= Aux
       !  END DO
       !  DO k= 1,nuz,Step
       !   WRITE(1,*) Feuz(k)
       !  END DO
       !  CLOSE(1)

    case("IL1")
       open(1,FILE='IL1D')
       do k= 1,nqz
          do i= 1,nqx
             Vintqx(i)= ILm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm(k)= Aux
          do i= 1,nqx
             Vintqx(i)= ILp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp(k)= Aux
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzm(nqz+1-k)
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzp(k)
       end do
       close(1)

    case("IS1")
       open(1,FILE='IS1D')
       do k= 1,nqz
          Qz= VQz(k)
          do i= 1,nqx
             Qx= VQx(i)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             Vintqx(i)= Muq*ISm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm2(k)= Aux
          do i= 1,nqx
             Qx= VQx(i)
             Q= sqrt(Qx**2+Qz**2)
             Muq= Q**3*AA/2._wp
             Vintqx(i)= Muq*ISp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp2(k)= Aux
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzm2(nqz+1-k)
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzp2(k)
       end do
       close(1)

    case("IT1")
       open(1,FILE='IT1D')
       do k= 1,nqz
          do i= 1,nqx
             Vintqx(i)= ITm(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzm(k)= Aux
          do i= 1,nqx
             Vintqx(i)= ITp(i,k)
          end do
          call Simpson(VQx,Vintqx,nqx,Aux)
          Iqzp(k)= Aux
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzm(nqz+1-k)
       end do
       do k= 1,nqz,Step
          write(1,*) Iqzp(k)
       end do
       close(1)

    case("ITQ")
       open(1,FILE='ITQ1')
       open(2,FILE='ITQ2')
       Dq= (Qzf-Qmin)/(nph-1)
       do i= 1,nph
          Q= Qmin+(i-1)*Dq
          VQ(i)= Q
       end do
       do i= 1,nph
          Q= VQ(i)
          do k= 1,nph
             Qxaux= Q*sin(VPhip(k))
             Qzaux= Q*cos(VPhip(k))
             call Aitp2d(nqx,nqz,VQx,VQz,ITp,Qxaux,Qzaux,Aux)
             Vintq(k)= Aux
             call Aitp2d(nqx,nqz,VQx,VQz,ITm,Qxaux,Qzaux,Aux)
             Vintq(k)= Vintq(k)+Aux
          end do
          call Simpson(VPhip,Vintq,nph,Aux)
          write(1,*) Aux
          write(2,*) Aux*Q
       end do
       close(1)
       close(2)

    case DEFAULT

    end select

    return

  end subroutine Output2

  !---------------------------------------------------------------------
  ! Mathematical Routines:
  !---------------------------------------------------------------------

  subroutine Split(Dux,Duz,DTau)

    ! See Notes D-10.

    use,intrinsic:: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: Auxx,Auxz,Auxxx,Auxzz
    real(wp) :: DTau,Dux,Duz,U,Ux,Uz
    real(wp), dimension(nux) :: Alpha1,Beta1,Gamma1,Psi1
    real(wp), dimension(nuz) :: Alpha2,Beta2,Gamma2,Psi2
    real(wp), dimension(nux) :: Faux1,Faux1b
    real(wp), dimension(nuz) :: Faux2,Faux2b
    real(wp), dimension(nux) :: Baux1,Gaux1
    real(wp), dimension(nuz) :: Baux2,Gaux2
    real(wp), dimension(nux,nuz) :: Dfdux,Dfduz
    real(wp), dimension(nux,nuz) :: Dfduxdux,Dfduxduz
    real(wp), dimension(nux,nuz) :: Dfduzdux,Dfduzduz
    real(wp), dimension(nux,nuz) :: DAxdux
    real(wp), dimension(nux,nuz) :: DAzduz
    real(wp), dimension(nux,nuz) :: DDxxdux,DDxzdux
    real(wp), dimension(nux,nuz) :: DDzxduz,DDzzduz
    real(wp), dimension(nux,nuz) :: CAx,CAz
    real(wp), dimension(nux,nuz) :: CDxx,CDxz,CDzx,CDzz
    real(wp), dimension(nux) :: Auxi
    real(wp), dimension(nuz) :: Auxk
    real(wp), dimension(nux,nuz) :: Fex,Fez
    integer :: i,k

    ! CALL Coef_A
    call Coef_D
    Auxx= DTau/4._wp/Dux
    Auxz= DTau/4._wp/Duz
    Auxxx= DTau/2._wp/Dux/Dux
    Auxzz= DTau/2._wp/Duz/Duz

    CAx= Ax+ColAx
    CAz= Az+ColAz
    CDxx= Dxx+ColDxx
    CDxz= Dxz+ColDxz
    CDzx= Dzx+ColDxz
    CDzz= Dzz+ColDzz

    call Derivx5p2d(nux,nuz,VUx,CDxx,DDxxdux)
    call Derivx5p2d(nux,nuz,VUx,CDxz,DDxzdux)
    call Derivy5p2d(nux,nuz,VUz,CDzx,DDzxduz)
    call Derivy5p2d(nux,nuz,VUz,CDzz,DDzzduz)
    call Derivx5p2d(nux,nuz,VUx,CAx,DAxdux)
    call Derivy5p2d(nux,nuz,VUz,CAz,DAzduz)

    select case(DerivLn)
    case("Yes")
       call Derivxy5pln2d(nux,nuz,VUx,VUz,Fe,Dfdux,Dfduz)
        Dfdux(1,:)= 0._wp
    case("No ")
       call Derivxy5p2d(nux,nuz,VUx,VUz,Fe,Dfdux,Dfduz)
        Dfdux(1,:)= 0._wp
    end select
    call Derivxy5p2d(nux,nuz,VUx,VUz,Dfdux,Dfduxdux,Dfduxduz)
    call Derivxy5p2d(nux,nuz,VUx,VUz,Dfduz,Dfduzdux,Dfduzduz)
     Dfduxdux(1,:)= 0._wp
     Dfduxduz(1,:)= 0._wp
     Dfduzdux(1,:)= 0._wp

    ! Operator L-ux:
    do k= 1,nuz
       forall(i=2:nux-1)
          Auxi(i)= (tanh(sqrt((VUx(i))**2+(VUz(k))**2)))**2
          Alpha1(i)= ( ( CAx(i,k)+DDxxdux(i,k) ) * Auxx - CDxx(i,k)*Auxxx ) *Auxi(i)
          Beta1(i)= 1._wp + (- DAxdux(i,k)*DTau/2._wp + 2._wp*CDxx(i,k)*Auxxx ) *Auxi(i)
          Gamma1(i)= (- ( CAx(i,k)+DDxxdux(i,k) ) * Auxx - CDxx(i,k)*Auxxx ) *Auxi(i)
          Psi1(i)= Fe(i,k) + ( DAxdux(i,k)*Fe(i,k)*(DTau/2._wp) &
               + ( CDxz(i,k)*Dfduxduz(i,k) + DDxzdux(i,k)*Dfduz(i,k) ) * (DTau) &
               + ( CAx(i,k)+DDxxdux(i,k) ) * Auxx * (Fe(i+1,k)-Fe(i-1,k)) &
               + CDxx(i,k)*Auxxx * (Fe(i+1,k)-2._wp*Fe(i,k)+Fe(i-1,k)) ) *Auxi(i)
          Faux1b(i)= Fe(i,k)
       end forall

       i= 1
       Auxi(i)= (tanh(sqrt((VUx(i))**2+(VUz(k))**2)))**2
       Alpha1(i)= 0._wp
       Beta1(i)= 1._wp  + ( - DAxdux(i,k)*DTau/2._wp ) * Auxi(i) !&
            ! + ( CAx(i,k)+DDxxdux(i,k) ) * (DTau/2._wp/Dux) ) * Auxi(i)
       Gamma1(i)= 0._wp ! - ( CAx(i,k)+DDxxdux(i,k) ) * (DTau/2._wp/Dux) * Auxi(i)
       Psi1(i)= Fe(i,k)  + ( DAxdux(i,k)*Fe(i,k)*(DTau/2._wp) ) * Auxi(i) !&
            ! + ( CDxz(i,k)*Dfduxduz(i,k) + DDxzdux(i,k)*Dfduz(i,k) ) &
            ! * (DTau)) * Auxi(i)       
            ! + ( CAx(i,k)+DDxxdux(i,k) ) * (DTau/2._wp/Dux) * (Fe(i+1,k)-Fe(i,k)) &
            ! + CDxx(i,k)*(DTau/Dux/Dux) * (Fe(i+2,k)-2._wp*Fe(i+1,k)+Fe(i,k))
       Faux1b(i)= Fe(i,k)

       ! i= 2
       ! Psi1(i)= Psi1(i)-Alpha1(i)*Fe(i-1,k)
       ! Alpha1(i)= 0._wp
       ! Faux1b(i)= Fe(i,k)

       i= nux
       Auxi(i)= (tanh(sqrt((VUx(i))**2+(VUz(k))**2)))**2
       Alpha1(i)= ( CAx(i,k)+DDxxdux(i,k) ) * (DTau/2._wp/Dux) *Auxi(i)
       Beta1(i)= 1._wp + (- DAxdux(i,k)*DTau/2._wp &
            - ( CAx(i,k)+DDxxdux(i,k) ) * (DTau/2._wp/Dux) ) *Auxi(i)
       Gamma1(i)= 0._wp
       Psi1(i)= Fe(i,k) + ( DAxdux(i,k)*Fe(i,k)*(DTau/2._wp) &
            + ( CDxz(i,k)*Dfduxduz(i,k) + DDxzdux(i,k)*Dfduz(i,k) ) * (DTau) &
            + ( CAx(i,k)+DDxxdux(i,k) ) * (DTau/2._wp/Dux) * (Fe(i,k)-Fe(i-1,k)) &
            + CDxx(i,k)*(DTau/Dux/Dux) * (Fe(i,k)-2._wp*Fe(i-1,k)+Fe(i-2,k)) ) *Auxi(i)
       Faux1b(i)= Fe(i,k)

       call Tridag(1,nux,Alpha1,Beta1,Gamma1,Psi1,Faux1,Baux1,Gaux1)

       ! ! Temporary boundary condition:
        Faux1(1)= Fe(1,k)
        Faux1b(1)= Faux1(2)

       ! Correcting instabilities:
       call Cor_Ampli_2(Faux1b,Faux1,nux)

       ! ! Boundary condition:
       Uz=VUz(k)
       Ux=VUx(1)
       U=sqrt(Ux**2+Uz**2)
       if(abs(U).lt.Ucrit) then
          Faux1(1)= Faux1(1)
       else
          Faux1(1)= Faux1(2)
       end if

       forall(i=1:nux)
          Fex(i,k)= Faux1(i)
       end forall
    end do

    ! Operator L-uz:
    do i= 1,nux
       forall(k=2:nuz-1)
          Auxk(k)= (tanh(sqrt((VUx(i))**2+(VUz(k))**2)))**2
          Alpha2(k)= ( ( CAz(i,k)+DDzzduz(i,k) ) * Auxz - CDzz(i,k)*Auxzz ) *Auxk(k)
          Beta2(k)= 1._wp + (- DAzduz(i,k)*DTau/2._wp*Auxk(k) + 2._wp*CDzz(i,k)*Auxzz) *Auxk(k)
          Gamma2(k)= (- ( CAz(i,k)+DDzzduz(i,k) ) * Auxz - CDzz(i,k)*Auxzz ) *Auxk(k)
          Psi2(k)= Fe(i,k) + ( DAzduz(i,k)*Fe(i,k)*(DTau/2._wp) & 
               + ( CDzx(i,k)*Dfduzdux(i,k) + DDzxduz(i,k)*Dfdux(i,k) ) * (DTau) &
               + ( CAz(i,k)+DDzzduz(i,k) ) * Auxz * (Fe(i,k+1)-Fe(i,k-1)) &
               + CDzz(i,k)*Auxzz * (Fe(i,k+1)-2._wp*Fe(i,k)+Fe(i,k-1)) )*Auxk(k)
          Faux2b(k)= Fe(i,k)
       end forall
       k= 1
       Auxk(k)= (tanh(sqrt((VUx(i))**2+(VUz(k))**2)))**2
       Alpha2(k)= 0._wp
       Beta2(k)= 1._wp + ( - DAzduz(i,k)*DTau/2._wp &
            + ( CAz(i,k)+DDzzduz(i,k) ) * (DTau/2._wp/Duz) ) *Auxk(k)
       Gamma2(k)= - ( CAz(i,k)+DDzzduz(i,k) ) * (DTau/2._wp/Duz)*Auxk(k)
       Psi2(k)= Fe(i,k) + ( DAzduz(i,k)*Fe(i,k)*(DTau/2._wp) &
            + ( CDzx(i,k)*Dfduzdux(i,k) + DDzxduz(i,k)*Dfdux(i,k) ) * (DTau) &
            + ( CAz(i,k)+DDzzduz(i,k) ) * (DTau/2._wp/Duz) * (Fe(i,k+1)-Fe(i,k)) &
            + CDzz(i,k)*(DTau/Duz/Duz) * (Fe(i,k+2)-2._wp*Fe(i,k+1)+Fe(i,k)) ) *Auxk(k)
       Faux2b(k)= Fe(i,k)
       k= nuz
       Auxk(k)= (tanh(sqrt((VUx(i))**2+(VUz(k))**2)))**2
       Alpha2(k)= ( CAz(i,k)+DDzzduz(i,k) ) * (DTau/2._wp/Duz)*Auxk(k)
       Beta2(k)= 1._wp + (- DAzduz(i,k)*DTau/2._wp &
            - ( CAz(i,k)+DDzzduz(i,k) ) * (DTau/2._wp/Duz) ) *Auxk(k)
       Gamma2(k)= 0._wp
       Psi2(k)= Fe(i,k) + ( DAzduz(i,k)*Fe(i,k)*(DTau/2._wp) &
            + ( CDzx(i,k)*Dfduzdux(i,k) + DDzxduz(i,k)*Dfdux(i,k) ) * (DTau) &
            + ( CAz(i,k)+DDzzduz(i,k) ) * (DTau/2._wp/Duz) * (Fe(i,k)-Fe(i,k-1)) &
            + CDzz(i,k)*(DTau/Duz/Duz) * (Fe(i,k)-2._wp*Fe(i,k-1)+Fe(i,k-2)) ) *Auxk(k)
       Faux2b(k)= Fe(i,k)

       call Tridag(1,nuz,Alpha2,Beta2,Gamma2,Psi2,Faux2,Baux2,Gaux2)

       ! Correcting instabilities:
       call Cor_Ampli_2(Faux2b,Faux2,nuz)

       forall(k=1:nuz)
          Fez(i,k)= Faux2(k)
       end forall
    end do

    Fe= Fex+Fez-Fe
 
    return
  end subroutine Split

  subroutine Evol_Iwave(DTau,Tau,WaveType)
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    ! See Notes B-11.
    use Common_Params
    use Common_Arrays
    use Math_Constants
    use Phys_Constants
    implicit none
    real(wp) :: DTau,Tau
    integer :: i,k,sigma
    integer :: Nvar,ind
    real(wp), dimension(nux,nuz) :: Dfdux,Dfduz
    real(wp), dimension(nqx,nqz) :: Iwave,Iwavenew
    real(wp), dimension(nqx,nqz) :: CoefA,CoefB
    !REAL(WP), DIMENSION(nqx,nqz) :: BremssL,BremssS
    real(wp), dimension(nqx*nqz) :: Ft,Dfdt,Ftout
    character(LEN=6) :: WaveType

    select case(DerivLn)
    case("Yes")
       call Derivxy5pln2d(nux,nuz,VUx,VUz,Fe,Dfdux,Dfduz)
    case("No ")
       call Derivxy5p2d(nux,nuz,VUx,VUz,Fe,Dfdux,Dfduz)
    end select
    Dfdux(1,:)= 0._wp

    select case(WaveType)
    case("Lwavep")
       sigma= 1
       Iwave= ILp
       call Coef_Lwave(sigma,Dfdux,Dfduz,CoefA,CoefB)
    case("Lwavem")
       sigma= -1
       Iwave= ILm
       call Coef_Lwave(sigma,Dfdux,Dfduz,CoefA,CoefB)
    case("Swavep")
       sigma= 1
       Iwave= ISp
       call Coef_Swave(sigma,Dfdux,Dfduz,CoefA,CoefB)
    case("Swavem")
       sigma= -1
       Iwave= ISm
       call Coef_Swave(sigma,Dfdux,Dfduz,CoefA,CoefB)
    case("Twavep")
       sigma= 1
       Iwave= ITp
       call Coef_Twave(sigma,CoefA,CoefB)
    case("Twavem")
       sigma= -1
       Iwave= ITm
       call Coef_Twave(sigma,CoefA,CoefB)
    end select
    CoefARK= CoefA
    CoefBRK= CoefB

    !H= DTau
    Nvar= nqx*nqz
    ind= 0
    do i= 1,nqx
       do k= 1,nqz
          ind= ind+1
          Ft(ind)= Iwave(i,k)
       end do
    end do
    !IF( (Tau+H-Tau2)*(Tau+H-Tau1) > 0._wp ) H=Tau2-Tau
    call DERIVS(Tau+DTau,Ft,Dfdt,Nvar)
    call RK4(Ft,Dfdt,Nvar,Tau,DTau,Ftout)
    call Cor_Ampli_2(Ft,Ftout,Nvar)
    Ft= Ftout
    ind= 0
    do i= 1,nqx
       do k= 1,nqz
          ind= ind+1
          Iwavenew(i,k)= Ft(ind)
       end do
    end do

    select case(WaveType)
    case("Lwavep")
       ILp= Iwavenew
    case("Lwavem")
       ILm= Iwavenew
    case("Swavep")
       ISp= Iwavenew
    case("Swavem")
       ISm= Iwavenew
    case("Twavep")
       ITp= Iwavenew
    case("Twavem")
       ITm= Iwavenew
    end select

    return
  end subroutine Evol_Iwave

  !
  ! Subroutine for solving a system of linear simultaneous
  ! equations having a tridiagonal coefficient matrix.
  ! The equations are numbered from If through Lim, and their
  ! sub-diagonal, diagonal, and super-diagonal coefficients 
  ! are stored in the arrays A, B, and C. The computed solution
  ! vector V(If)...V(Lim) is stored in the array V.
  !* Carnahan, Luther and Wilkes, Applied Numerical Methods,
  !  John Wiley, 1969, pag. 446.
  !
  !* Mod. Oct/90: Beta and Gamma appear in the variable list (LFZ).
  !* Mod. Aug/06: Adapted to Fortran 95.
  !* Mod. Aug/10: Fortran 95, FORALL instead of DO.
  !
  subroutine Tridag(if,k,A,B,C,D,V,Beta,Gamma)
    ! 	PARAMETER(K=100)
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    implicit none
    integer, intent(IN) :: if,k 
    integer :: Lim,Ifp1,Last,i,j
    real(wp), intent(in) :: A(k),B(k),C(k),D(k)
    real(wp), intent(out) :: V(k)
    real(wp) :: Beta(k),gamma(k)
    !
    Lim= k
    !... Compute intermediate arrays Beta and Gamma ...
    Beta(if)=B(if)
    gamma(if)=D(if)/Beta(if)
    Ifp1=if+1
    do i=Ifp1,Lim
       Beta(i)=B(i)-A(i)*C(i-1)/Beta(i-1)
       gamma(i)=(D(i)-A(i)*gamma(i-1))/Beta(i)
    end do

    !... Compute final solution vector V .....
    V(Lim)=gamma(Lim)
    Last=Lim-if
    do j=1,Last
       i=Lim-j
       V(i)=gamma(i)-C(i)*V(i+1)/Beta(i)
    end do
    return
  end subroutine Tridag

  subroutine DERIVS(Tau,Ft,Dfdt,N)
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use Common_Params
    use Common_Arrays
    implicit none
    integer, parameter :: NMAX=6561    ! NMAX= 81*81
    integer, intent(in) :: N
    integer :: i,k,ind
    real(wp) :: Tau
    real(wp), dimension(N) :: Ft,DFdt
    ind= 0
    do i= 1,nqx
       do k= 1,nqz
          ind= ind+1
          Dfdt(ind)= CoefARK(i,k)+CoefBRK(i,k)*Ft(ind)
       end do
    end do
    Tau= 0._wp    ! Just to avoid compilation warning
    return
  end subroutine DERIVS

  subroutine RK4(Y,DYDX,N,X,H,YOUT)
    ! Given values for N variables Y and their derivatives DYDX
    ! known at X, use the fourth-order Runge-Kutta method to advan-
    ! ce the solution over an interval H and return the incremen-
    ! ted variables as YOUT, which need not be a distinct array from
    ! Y. The user supplies the subroutine DERIVS(X,Y,DYDX,N) which
    ! returns derivatives DYDX at X.
    ! See NUMERICAL RECIPES.
    ! Version Fortran 95, Feb. 2007.
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    implicit none
    integer, parameter :: NMAX=6561    ! NMAX= 81*81
    integer, intent(in) :: N
    integer :: I
    real(wp), dimension(N) :: Y,DYDX,YOUT
    real(wp), dimension(NMAX) :: YT,DYT,DYM
    real(wp) :: X,H,HH,H6,XH

    HH=H*0.5_wp
    H6=H/6._wp
    XH=X+HH
    ! CALL  DERIVS(X,Y,DYDX)
    do I=1,N
       YT(I)=Y(I)+HH*DYDX(I)
    end do
    call DERIVS(XH,YT,DYT,N)
    do I=1,N
       YT(I)=Y(I)+HH*DYT(I)
    end do
    call DERIVS(XH,YT,DYM,N)
    do I=1,N
       YT(I)=Y(I)+H*DYM(I)
       DYM(I)=DYT(I)+DYM(I)
    end do
    call DERIVS(X+H,YT,DYT,N)
    do I=1,N
       YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2._wp*DYM(I))
    end do
    return
  end subroutine RK4

  subroutine Cor_Ampli(Vf,n)
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use Math_Constants
    implicit none
    integer :: n,i
    real(wp), dimension (n) :: Vf(n)

    ! Corrects negative amplitudes when they are not acceptable.
    ! For a given real vector, replaces negative values by zero.
    ! Version Fortran 95, Aug 2006.

    do i= 1,n
       if ( Vf(i) < 0._wp ) then
          Vf(i)= 0._wp
       else
       end if
    end do

    return
  end subroutine Cor_Ampli

  subroutine Cor_Ampli_2(Vf,Vfnew,n)
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    use Math_Constants
    implicit none
    integer :: n,i
    real(wp), dimension (n) :: Vf(n),Vfnew(n)
    ! Corrects negative amplitudes when they are not acceptable.
    ! For a given real vector "Vfnew", replaces negative values by the previous 
    ! value ("Vf").

    do i= 1,n
       ! IF ( Vfnew(i) < EpsMin ) THEN
       !  Vfnew(i)= Vf(i) 
       ! ELSE
       ! END IF
       if ( Vfnew(i) < 0._wp ) then
          !  Vfnew(i)= 0._wp 
          Vfnew(i)= Vf(i) 
       else
       end if
    end do

    return
  end subroutine Cor_Ampli_2

  subroutine Aitp1d2(nx,Vx,Fx,Xp,Fp,i)
    ! Uses linear interpolation in order to obtain the value of a function
    ! F, at the point Xp. The function F is given as a set of points Fx(x).
    ! Uses subroutine Locate (Numerical Recipes, P. 96)
    ! Version Fortran 95, Aug, 2006.
    ! Version of Aitp1d, without the call to Locate, which is called outside.
    use,intrinsic :: iso_fortran_env, only: wp=>real64
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
    use,intrinsic :: iso_fortran_env, only: wp=>real64
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

  subroutine Aitp2d(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp)
    ! Uses linear interpolation in order to obtain the value of a function
    ! F, at the point (Xp,Yp). 
    ! The function F is given as a set of points Fxy(x,y).
    ! Uses subroutine Locate (Numerical Recipes, P. 96)
    ! Version Fortran 95, Sep., 2006.
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    implicit none
    real(wp), dimension(nx) :: Vx
    real(wp), dimension(ny) :: Vy
    real(wp), dimension(nx,ny) :: Fxy
    real(wp) :: Xp,Yp,Fp
    real(wp) :: F1,F2,F3,F4,T,U
    integer :: nx,ny,i,j

    call Locate(Vx,nx,Xp,i)
    call Locate(Vy,ny,Yp,j)

    if ( (i==0 .or. i==nx) .or. (j==0 .or. j==ny) ) then
       FP= 0._wp
    else
       F1= Fxy(i,j)
       F2= Fxy(i+1,j)
       F3= Fxy(i+1,j+1)
       F4= Fxy(i,j+1)
       T= ( Xp-Vx(i) )/( Vx(i+1)-Vx(i) )
       U= ( Yp-Vy(j) )/( Vy(j+1)-Vy(j) )
       FP= (1._wp-T)*(1._wp-U)*F1 + T*(1._wp-U)*F2 + T*U*F3 + (1._wp-T)*U*F4
       !write(11,2001)i,Vx(i),Vx(i+1),Xp,T
       !write(12,2001)j,Vy(j),Vy(j+1),Yp,U
       !2001 FORMAT(1x,i4,1x,4(f8.3,1x))
    end if
    return
  end subroutine Aitp2d

  subroutine Aitp2d2(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp,i,j)
    ! Uses linear interpolation in order to obtain the value of a function
    ! F, at the point (Xp,Yp). 
    ! The function F is given as a set of points Fxy(x,y).
    ! Uses subroutine Locate (Numerical Recipes, P. 96)
    ! Version Fortran 95, Sep., 2006.
    ! Version of Aitp2d, without the calls to Locate, which is called outside.
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    implicit none
    real(wp), dimension(nx) :: Vx
    real(wp), dimension(ny) :: Vy
    real(wp), dimension(nx,ny) :: Fxy
    real(wp) :: Xp,Yp,Fp
    real(wp) :: F1,F2,F3,F4,T,U
    integer :: nx,ny,i,j

    !CALL Locate(Vx,nx,Xp,i)
    !CALL Locate(Vy,ny,Yp,j)

    if ( (i==0 .or. i==nx) .or. (j==0 .or. j==ny) ) then
       FP= 0.
    else
       F1= Fxy(i,j)
       F2= Fxy(i+1,j)
       F3= Fxy(i+1,j+1)
       F4= Fxy(i,j+1)
       T= ( Xp-Vx(i) )/( Vx(i+1)-Vx(i) )
       U= ( Yp-Vy(j) )/( Vy(j+1)-Vy(j) )
       FP= (1._wp-T)*(1._wp-U)*F1 + T*(1._wp-U)*F2 + T*U*F3 + (1._wp-T)*U*F4
    end if
    return
  end subroutine Aitp2d2

  subroutine Aitp2d2b(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp,i,j)
    ! Uses linear interpolation in order to obtain the value of a function
    ! F, at the point (Xp,Yp). 
    ! The function F is given as a set of points Fxy(x,y).
    ! Uses subroutine Locate (Numerical Recipes, P. 96)
    ! Version Fortran 95, Sep., 2006.
    ! Version of Aitp2d, without the calls to Locate, which is called outside.
    ! Corrects the interpolation at the extreme points, making the value at the
    ! points outside of the grid equal to the nearest point inside the grid.
    ! It is only for interpolation using similar grids, to avoid errors at the edge.
    use,intrinsic :: iso_fortran_env, only: wp=>real64
    implicit none
    real(wp), dimension(nx) :: Vx
    real(wp), dimension(ny) :: Vy
    real(wp), dimension(nx,ny) :: Fxy
    real(wp) :: Xp,Yp,Fp
    real(wp) :: F1,F2,F3,F4,T,U
    integer :: nx,ny,i,j

    !CALL Locate(Vx,nx,Xp,i)
    !CALL Locate(Vy,ny,Yp,j)

    if (i==0) then
       if (j==0) then
          FP= Fxy(1,1)
       else
          if (j==ny) then
             FP= Fxy(1,ny)
          else
             FP= Fxy(1,j)+ (Yp-Vy(j))/(Vy(j+1)-Vy(j))*(Fxy(1,j+1)-Fxy(1,j)) 
          end if
       end if
    else
       if (i==nx) then
          if (j==0) then
             FP= Fxy(nx,1)
          else
             if (j==ny) then
                FP= Fxy(nx,ny)
             else
                FP= Fxy(nx,j)+ (Yp-Vy(j))/(Vy(j+1)-Vy(j))*(Fxy(nx,j+1)-Fxy(nx,j)) 
             end if
          end if
       else
          if (j==0) then
             FP= Fxy(i,1)+ (Xp-Vx(i))/(Vx(i+1)-Vx(i))*(Fxy(i+1,1)-Fxy(i,1)) 
          else
             if (j==ny) then
                FP= Fxy(i,ny)+ (Xp-Vx(i))/(Vx(i+1)-Vx(i))*(Fxy(i+1,ny)-Fxy(i,ny)) 
             else
                F1= Fxy(i,j)
                F2= Fxy(i+1,j)
                F3= Fxy(i+1,j+1)
                F4= Fxy(i,j+1)
                T= ( Xp-Vx(i) )/( Vx(i+1)-Vx(i) )
                U= ( Yp-Vy(j) )/( Vy(j+1)-Vy(j) )
                FP= (1._wp-T)*(1._wp-U)*F1 + T*(1._wp-U)*F2 + T*U*F3 + (1._wp-T)*U*F4
             end if
          end if
       end if
    end if
    return
  end subroutine Aitp2d2b

  subroutine Simpson(Vx,F,N,Res)
    ! Version Fortran 95, Aug 2006.
    ! The Jun 2018 version is suitable for
    ! both odd and even values of N.
    ! If N is odd, the subroutine uses the
    ! regular Simpson formula.
    ! If N is even, the subroutine uses the
    ! regular Simpson rule for i=2:N-4 and
    ! the 3/8 Simpson formula for i=N-4,N.
    !! use precision_m
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    integer :: N,Nm1,Nm2,Nm3,Nm4,Nm5,i
    real(wp), dimension(N) :: F,Vx
    real(wp) :: Res,Res1,H

    if ( N < 5 ) then
       open(1,FILE='Warning_Simpson.wt')
       if(mod(N,2) /= 0) then
          write(1,*) ' N must be larger or equal 5!'
       else
          write(1,*) ' N must be larger or equal 6!'
       end if
       close(1)
       stop
    else
    end if

    H= ( Vx(N) - Vx(1) )/(N-1)

    Nm1= N-1
    Nm2= N-2
    Nm3= N-3
    Nm4= N-4
    Nm5= N-5
    Res= 0._wp
    Res1= 0._wp
    if(mod(N,2) .ne. 0) then
       ! print*, "odd"
       do i= 2,Nm1,2
          Res= Res + 4._wp*F(i)
       end do
       do i= 3,Nm2,2
          Res= Res + 2._wp*F(i)
       end do
       Res= Res + F(1) + F(N)
       Res= Res*H/3._wp
    else
       ! print*, "even"
       do i= 2,Nm4,2
          Res= Res + 4._wp*F(i)
       end do
       do i= 3,Nm5,2
          Res= Res + 2._wp*F(i)
       end do
       Res= Res + F(1) + F(Nm3)
       Res1= F(Nm3)+ 3._wp*(F(Nm2)+F(Nm1))+F(N)
       Res= Res*H/3._wp + Res1*3._wp*H/8._wp
    end if
    return

  end subroutine Simpson

  ! subroutine Simpson(Vx,F,N,Res)
  !   ! Version Fortran 95, Aug 2006.
  !   implicit none
  !   integer :: N,Nm1,Nm2,i
  !   real(wp), dimension(N) :: F,Vx
  !   real(wp) :: Res,H

  !   ! N must be odd.
  !   ! The intervals must be equal.

  !   H= ( Vx(N) - Vx(1) )/(N-1)
  !   Res=0._wp
  !   if ( N < 5 ) then
  !      open(1,FILE='Warning_Simpson.wt')
  !      write(1,*) ' N must be larger or equal 5, and odd! '
  !      close(1)
  !      stop
  !   else
  !   end if
  !   Nm1= N-1
  !   Nm2= N-2
  !   do i= 2,Nm1,2
  !      Res= Res + 4.*F(i)
  !   end do
  !   do i= 3,Nm2,2
  !      Res= Res + 2._wp*F(i)
  !   end do
  !   Res= Res + F(1) + F(N)
  !   Res= Res*H/3.
  !   return

  ! end subroutine Simpson

 !  subroutine Qsimp(func,a,b,s)
!     implicit none
!     ! Version 22/Aug/2015
!     real(wp) :: a,b,s
!     real(wp) :: EPS= 1.E-4
!     integer :: j
!     integer :: jmax= 15
!     real(wp) :: os,ost,st
!     !      EXTERNAL func
!     interface 
!        real(wp) function FUNC(X)
!          real(wp), intent(in) :: X
!        end function FUNC
!     end interface
!     ost=-1.E30
!     os= -1.E30
!     do j=1,JMAX
!        call trapzd(func,a,b,st,j)
!        s=(4.e0_wp*st-ost)/3.e0_wp
!        if (j.gt.5) then
!           if (abs(s-os).lt.EPS*abs(os).or. &
!                (s.eq.0..and.os.eq.0.)) return
!        endif
!        os=s
!        ost=st
!     enddo
!     open(98,file="WARNING_qsimp")
!     write(98,*) " Too many steps in qsimp! "
!     close(98)
!     !      STOP
!   end subroutine Qsimp


!   subroutine Trapzd(func,a,b,s,n)
!     implicit none
!     ! Version 22/Aug/2015
!     integer :: n
!     integer :: it,j
!     real(wp) :: a,b,s
!     real(wp) :: del,suma,tnm,x
!     !      EXTERNAL func
!     interface 
!        real(wp) function FUNC(X)
!          real(wp), intent(in) :: X
!        end function FUNC
!     end interface
!     if (n.eq.1) then
!        s=0.5E0*(b-a)*(func(a)+func(b))
!     else
!        it=2**(n-2)
!        tnm=it
!        del=(b-a)/tnm
!        x=a+0.5E0*del
!        suma=0.e0_wp
!        do j=1,it
!           suma=suma+func(x)
!           x=x+del
!        enddo
!        s=0.5E0*(s+(b-a)*suma/tnm)
!     endif
!     return
!   end subroutine Trapzd

!   subroutine Qsimpb(func,a,b,s)
!     implicit none
!     ! Version 22/Aug/2015
!     ! To be used if "a' is finite, positive, and "b" goes to infinity, or if
!     ! "b" is finite, negative, and "a" goes to -infinity.
!     real(wp) :: a,b,s
!     real(wp) :: EPS= 1.E-4
!     integer :: j
!     integer :: jmax= 15
!     real(wp) :: os,ost,st
!     !      EXTERNAL func
!     interface 
!        real(wp) function FUNC(X)
!          real(wp), intent(in) :: X
!        end function FUNC
!     end interface
!     ost=-1.E30
!     os= -1.E30
!     do j=1,JMAX
!        call midinf(func,a,b,st,j)
!        s=(4.e0_wp*st-ost)/3.e0_wp
!        if (j.gt.5) then
!           if (abs(s-os).lt.EPS*abs(os).or. &
!                (s.eq.0.e0_wp.and.os.eq.0.e0_wp)) return
!        endif
!        os=s
!        ost=st
!     enddo
!     open(98,file="WARNING_qsimpb")
!     write(98,*) " Too many steps in qsimpb! "
!     close(98)
!     !      STOP
!   end subroutine Qsimpb

!   subroutine MIDINF(FUNK,AA,BB,S,N) 
!     !	Numerical Recipes, Cambridge, 1989, p. 118:
!     !	This routine is an exact replacement for MIDPNT, i.e. return as S 
!     !	the Nth stage of refinement of the integral of FUNK from AA to BB,
!     !	except that the function is evaluated at evenly spaced points in 1/x
!     !	rather than in x. This allows the upper limit BB to be as large and
!     !	positive as the computer allows, or the lower limit AA to be as large
!     !	and negative, but not both. AA and BB must have the same sign.
!     ! Version 22/Aug/2015
!     implicit none
!     integer :: n
!     integer :: it,j
!     real(wp) :: aa,bb,s,a,b
!     real(wp) :: del,sum,tnm,x,ddel,func
!     !      EXTERNAL func
!     interface 
!        real(wp) function FUNK(X)
!          real(wp), intent(in) :: X
!        end function FUNK
!     end interface
!     FUNC(X)=FUNK(1.e0_wp/X)/X**2   ! This is a statement function which effects
!     ! the change of variable. 
!     B=1.e0_wp/AA 		! These two statements change the limits of
!     ! integration accordingly.
!     A=1.e0_wp/BB 
!     if (N.eq.1) then 		! From this point on, the routine is exactly
!        ! identical to MIDPNT.
!        S=(B-A)*FUNC(0.5E0*(A+B)) 
!        IT=1 
!     else 
!        TNM=IT 
!        DEL=(B-A)/(3.e0_wp*TNM) 
!        DDEL=DEL+DEL 
!        X=A+0.5E0*DEL 
!        SUM=0.e0_wp 
!        do 11 J=1,IT 
!           SUM=SUM+FUNC(X) 
!           X=X+DDEL 
!           SUM=SUM+FUNC(X) 
!           X=X+DEL 
! 11        continue 
!           S=(S+(B-A)*SUM/TNM)/3.e0_wp 
!           IT=3*IT 
!        endif 
!        return 
!      end subroutine MIDINF 

     subroutine Derivxy5p2d(nx,ny,Vx,Vy,Fxy,Dfdx,Dfdy)
       ! Version Jan. 2010: 
       ! Evaluates the x-derivative and the y-derivative of two-dimensional
       ! function Fxy(x,y), given by an array on (nx,ny) elements

       ! Version Fortran 95, Sept. 2006.
       ! Uses 5-point derivative for the internal points.
       ! Version Fortran 95, Aug. 2010 (uses FORALL instead of DO).
       ! See Abramowitz & Stegun, 1970, Eq. 25.3.6, for 'p=0'
       ! Based on "deriv5p.f", College Park, 16/Feb/2001 (L. F. Ziebell)
       ! Needs equally spaced points!
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       implicit none
       integer :: i,j,nx,ny
       integer :: nxm1,nym1,nxm2,nym2
       real(wp), dimension(nx) :: Vx
       real(wp), dimension(ny) :: Vy
       real(wp), dimension(nx,ny) :: Fxy(nx,ny),Dfdx(nx,ny),Dfdy(nx,ny)
       real(wp) :: Dx,Dy

       Dx= Vx(2)-Vx(1)
       Dy= Vy(2)-Vy(1)

       nxm1= nx-1
       nxm2= nx-2
       forall(j=1:ny)
          Dfdx(1,j)= (Fxy(2,j)-Fxy(1,j))/Dx
          Dfdx(2,j)= (Fxy(3,j)-Fxy(1,j))/(2._wp*Dx)
          Dfdx(nx,j)= (Fxy(nx,j)-Fxy(nxm1,j))/Dx
          Dfdx(nxm1,j)= (Fxy(nx,j)-Fxy(nxm2,j))/(2._wp*Dx)
          forall(i=3:nxm2)
             Dfdx(i,j)= ( (Fxy(i+1,j)-Fxy(i-1,j))*2._wp/3._wp &
                  - (Fxy(i+2,j)-Fxy(i-2,j))/12._wp) / Dx
          end forall
       end forall

       nym1= ny-1
       nym2= ny-2
       forall(i=1:nx)
          Dfdy(i,1)= (Fxy(i,2)-Fxy(i,1))/Dy
          Dfdy(i,2)= (Fxy(i,3)-Fxy(i,1))/(2._wp*Dy)
          Dfdy(i,ny)= (Fxy(i,ny)-Fxy(i,nym1))/Dy
          Dfdy(i,nym1)= (Fxy(i,ny)-Fxy(i,nym2))/(2._wp*Dy)
          forall(j=3:nym2)
             Dfdy(i,j)= ( (Fxy(i,j+1)-Fxy(i,j-1))*2._wp/3._wp &
                  - (Fxy(i,j+2)-Fxy(i,j-2))/12._wp) / Dy
          end forall
       end forall

       return
     end subroutine Derivxy5p2d
     !
     !
     subroutine Derivxy5pln2d(nx,ny,Vx,Vy,Fxy,Dfdx,Dfdy)
       ! Version Jan. 2010: Introduces logarithm approach to improve derivatives
       ! of functions with small values.

       ! Version Fortran 95, Sept. 2006.
       ! Uses 5-point derivative for the internal points.
       ! See Abramowitz & Stegun, 1970, Eq. 25.3.6, for 'p=0'
       ! Based on "deriv5p.f", College Park, 16/Feb/2001 (L. F. Ziebell)
       ! Needs equally spaced points!

       ! Evaluates the x-derivative and the y derivative of two-dimensional
       ! function Fxy(x,y), given by an array on (nx,ny) elements
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       implicit none
       integer :: i,j,nx,ny
       integer :: nxm1,nym1,nxm2,nym2
       real(wp), dimension(nx) :: Vx
       real(wp), dimension(ny) :: Vy
       real(wp), dimension(nx,ny) :: Fxy(nx,ny),Dfdx(nx,ny),Dfdy(nx,ny)
       real(wp), dimension(nx,ny) :: Fold(nx,ny),Faux(nx,ny)
       integer, dimension(2) :: Imin
       real(wp) :: Dx,Dy
       real(wp) :: Fxymin

       Fold= Fxy
       Imin= minloc(Fxy)
       Fxymin= Fxy(Imin(1),Imin(2))
       Fxy= Fxy+abs(Fxymin)
       Fxy= Fxy+1.0E-30_wp
       Faux= Fxy
       Fxy= log(Fxy)

       Dx= Vx(2)-Vx(1)
       Dy= Vy(2)-Vy(1)

       nxm1= nx-1
       nxm2= nx-2
       forall(j=1:ny)
          Dfdx(1,j)= (Fxy(2,j)-Fxy(1,j))/Dx * Faux(1,j)
          Dfdx(2,j)= (Fxy(3,j)-Fxy(1,j))/(2._wp*Dx) * Faux(2,j)
          Dfdx(nx,j)= (Fxy(nx,j)-Fxy(nxm1,j))/Dx * Faux(nx,j)
          Dfdx(nxm1,j)= (Fxy(nx,j)-Fxy(nxm2,j))/(2._wp*Dx) *Faux(nxm1,j)
          forall(i=3:nxm2)
             Dfdx(i,j)= ( (Fxy(i+1,j)-Fxy(i-1,j))*2._wp/3._wp &
                  - (Fxy(i+2,j)-Fxy(i-2,j))/12._wp) / Dx * Faux(i,j)
          end forall
       end forall

       nym1= ny-1
       nym2= ny-2
       forall(i=1:nx)
          Dfdy(i,1)= (Fxy(i,2)-Fxy(i,1))/Dy * Faux(i,1)
          Dfdy(i,2)= (Fxy(i,3)-Fxy(i,1))/(2._wp*Dy) * Faux(i,2)
          Dfdy(i,ny)= (Fxy(i,ny)-Fxy(i,nym1))/Dy * Faux(i,ny)
          Dfdy(i,nym1)= (Fxy(i,ny)-Fxy(i,nym2))/(2._wp*Dy) * Faux(i,nym1)
          forall(j=3:nym2)
             Dfdy(i,j)= ( (Fxy(i,j+1)-Fxy(i,j-1))*2._wp/3._wp &
                  - (Fxy(i,j+2)-Fxy(i,j-2))/12._wp) / Dy * Faux(i,j)
          end forall
       end forall

       Fxy= Fold
       return
     end subroutine Derivxy5pln2d
     !

     subroutine Derivx5p2d(nx,ny,Vx,Fxy,Dfdx)
       ! Version Jan. 2010: 
       ! Evaluates the x-derivative of two-dimensional
       ! function Fxy(x,y), given by an array on (nx,ny) elements

       ! Based on "Derivxy5p2d", version Fortran 95, Sept. 2006.
       ! Uses 5-point derivative for the internal points.
       ! Version Fortran 95, Aug. 2010 (uses FORALL instead of DO).
       ! See Abramowitz & Stegun, 1970, Eq. 25.3.6, for 'p=0'
       ! Based on "deriv5p.f", College Park, 16/Feb/2001 (L. F. Ziebell)
       ! Needs equally spaced points!
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       implicit none
       integer :: i,j,nx,ny
       integer :: nxm1,nxm2
       real(wp), dimension(nx) :: Vx
       real(wp), dimension(nx,ny) :: Fxy(nx,ny),Dfdx(nx,ny)
       real(wp) :: Dx

       Dx= Vx(2)-Vx(1)

       nxm1= nx-1
       nxm2= nx-2
       forall(j=1:ny)
          Dfdx(1,j)= (Fxy(2,j)-Fxy(1,j))/Dx
          Dfdx(2,j)= (Fxy(3,j)-Fxy(1,j))/(2._wp*Dx)
          Dfdx(nx,j)= (Fxy(nx,j)-Fxy(nxm1,j))/Dx
          Dfdx(nxm1,j)= (Fxy(nx,j)-Fxy(nxm2,j))/(2._wp*Dx)
          forall(i=3:nxm2)
             Dfdx(i,j)= ( (Fxy(i+1,j)-Fxy(i-1,j))*2._wp/3._wp &
                  - (Fxy(i+2,j)-Fxy(i-2,j))/12._wp) / Dx
          end forall
       end forall
       return
     end subroutine Derivx5p2d

     !
     subroutine Derivy5p2d(nx,ny,Vy,Fxy,Dfdy)
       ! Version Jan. 2010: 
       ! Evaluates the y-derivative of two-dimensional
       ! function Fxy(x,y), given by an array on (nx,ny) elements

       ! Based on "Derivxy5p2d", version Fortran 95, Sept. 2006.
       ! Uses 5-point derivative for the internal points.
       ! Version Fortran 95, Aug. 2010 (uses FORALL instead of DO).
       ! See Abramowitz & Stegun, 1970, Eq. 25.3.6, for 'p=0'
       ! Based on "deriv5p.f", College Park, 16/Feb/2001 (L. F. Ziebell)
       ! Needs equally spaced points!
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       implicit none
       integer :: i,j,nx,ny
       integer :: nym1,nym2
       real(wp), dimension(ny) :: Vy
       real(wp), dimension(nx,ny) :: Fxy(nx,ny),Dfdy(nx,ny)
       real(wp) :: Dy

       Dy= Vy(2)-Vy(1)

       nym1= ny-1
       nym2= ny-2
       forall(i=1:nx)
          Dfdy(i,1)= (Fxy(i,2)-Fxy(i,1))/Dy
          Dfdy(i,2)= (Fxy(i,3)-Fxy(i,1))/(2._wp*Dy)
          Dfdy(i,ny)= (Fxy(i,ny)-Fxy(i,nym1))/Dy
          Dfdy(i,nym1)= (Fxy(i,ny)-Fxy(i,nym2))/(2._wp*Dy)
          forall(j=3:nym2)
             Dfdy(i,j)= ( (Fxy(i,j+1)-Fxy(i,j-1))*2._wp/3._wp &
                  - (Fxy(i,j+2)-Fxy(i,j-2))/12._wp) / Dy
          end forall
       end forall

       return
     end subroutine Derivy5p2d
     !
     !
     real(wp) function PERF(X)
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       ! Error function of argument X, based on PERF
       ! (from R. L. Meyer, Universite de Nancy, Franca)
       ! Version Feb. 18, 2009.
       implicit none
       !REAL(WP) :: PERFC
       real(wp) :: X,X2,SNUM,SDEN
       real(wp), dimension(5) :: P,Q
       data P /3.209377589138469472562E+03_wp,3.774852376853020208137E+02_wp,&
            1.138641541510501556495E+02_wp,3.161123743870565596947E+00_wp,&
            1.857777061846031526730E-01_wp/
       data Q /2.844236833439170622273E+03_wp,1.282616526077372275645E+03_wp,&
            2.440246379344441733056E+02_wp,2.360129095234412093499E+01_wp,1._wp/
       if (X <= 0.5_wp) then
          X2=X*X
          SNUM=X*(P(1)+X2*(P(2)+X2*(P(3)+X2*(P(4)+X2*P(5)))))
          SDEN=Q(1)+X2*(Q(2)+X2*(Q(3)+X2*(Q(4)+X2*Q(5))))
          PERF=SNUM/SDEN
       else
          PERF=1._wp-PERFC(X)
       end if
       return
     end function PERF

     real(wp) function PERFC(X)
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       ! Complementary error function of argument X, based on PERF 
       ! (from R. L. Meyer, Universite de Nancy, Franca)
       ! Version Feb. 18, 2009.
       implicit none
       !REAL(WP) :: PERF
       real(wp) :: X,X2,Y2,SNUM,SDEN
       real(wp), dimension(9) :: P,Q
       real(wp), dimension(6) :: R,S
       data P /1.23033935479799725272E+03_wp,2.05107837782607146532E+03_wp,&
            1.71204761263407058314E+03_wp,8.81952221241769090411E+02_wp,&
            2.98635138197400131132E+02_wp,6.61191906371416294775E+01_wp,&
            8.88314979438837594118E+00_wp,5.64188496988670089180E-01_wp,&
            2.15311535474403846343E-08_wp/
       data Q /1.23033935480374942043E+03_wp,3.43936767414372163696E+03_wp,&
            4.36261909014324715820E+03_wp,3.29079923573345962678E+03_wp,&
            1.62138957456669018874E+03_wp,5.37181101862009857509E+02_wp,&
            1.17693950891312499305E+02_wp,1.57449261107098347253E+01_wp,&
            1._wp/
       data R /-6.58749161529837803157E-04_wp,-1.60837851487422766278E-02_wp,&
            -1.25781726111229246204E-01_wp,-3.60344899949804439429E-01_wp,&
            -3.05326634961232344035E-01_wp,-1.63153871373020978498E-02_wp/
       data S /2.33520497626869185443E-03_wp,6.05183413124413191178E-02_wp,&
            5.27905102951428412248E-01_wp,1.87295284992346047209E+00_wp,&
            2.56852019228982242072E+00_wp,1._wp/
       if (X > 90._wp) then
          PERFC=0._wp
       else
          if (X <= 0.5_wp) then
             PERFC=1._wp-PERF(X)
          else
             X2=X*X
             if (X <= 4._wp) then
                SNUM=exp(-X2)*(P(1)+X*(P(2)+X*(P(3)+X*(P(4)+X*(P(5)+X*(P(6)&
                     +X*(P(7)+X*(P(8)+X*P(9)))))))))
                SDEN=Q(1)+X*(Q(2)+X*(Q(3)+X*(Q(4)+X*(Q(5)+X*(Q(6)+X*(Q(7)&
                     +X*(Q(8)+X*Q(9))))))))
                PERFC=SNUM/SDEN
             else
                Y2=1._wp/X2
                SNUM= R(1)+Y2*(R(2)+Y2*(R(3)+Y2*(R(4)+Y2*(R(5)+Y2*R(6)))))
                SDEN=S(1)+Y2*(S(2)+Y2*(S(3)+Y2*(S(4)+Y2*(S(5)+Y2*S(6)))))
                PERFC=exp(-X2)*(0.56418958354775628694_wp+SNUM/SDEN*Y2)/X
             end if
          end if
       end if
       return
     end function PERFC

     real(wp) function PERFCE(X)
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       ! Evaluates EXP(X**2)*ERFC(X
       ! Based on DERFCE(X)
       ! (from R. L. Meyer, Universite de Nancy, Franca)
       ! Version Oct. 09, 2015
       implicit none
       real(wp) :: X,X2,SNUM,SDEN,Y2
       real(wp), dimension(9) :: P,Q
       real(wp), dimension(6) :: R,S
       data P/1.23033935479799725272D+03,2.05107837782607146532D+03,&
            1.71204761263407058314D+03,8.81952221241769090411D+02,&
            2.98635138197400131132D+02,6.61191906371416294775D+01,&
            8.88314979438837594118D+00,5.64188496988670089180D-01,&
            2.15311535474403846343D-08/
       data Q/1.23033935480374942043D+03,3.43936767414372163696D+03,&
            4.36261909014324715820D+03,3.29079923573345962678D+03,&
            1.62138957456669018874D+03,5.37181101862009857509D+02,&
            1.17693950891312499305D+02,1.57449261107098347253D+01,&
            1./
       data R/-6.58749161529837803157D-04,-1.60837851487422766278D-02,&
            -1.25781726111229246204D-01,-3.60344899949804439429D-01,&
            -3.05326634961232344035D-01,-1.63153871373020978498D-02/
       data S/2.33520497626869185443D-03,6.05183413124413191178D-02,&
            5.27905102951428412248D-01,1.87295284992346047209D+00,&
            2.56852019228982242072D+00,1./
       if (X <= 0.5_wp) then
          PERFCE=(1._wp-PERF(X))*exp(X*X)
       else
          X2=X*X
          if (X <=4._wp) then
             SNUM= (P(1)+X*(P(2)+X*(P(3)+X*(P(4)+X*(P(5)+X*(P(6)+X*(P(7) &
                  +X*(P(8)+X*P(9)))))))))
             SDEN= Q(1)+X*(Q(2)+X*(Q(3)+X*(Q(4)+X*(Q(5)+X*(Q(6)+X*(Q(7)+X*(Q(8)+ &
                  X*Q(9))))))))
             PERFCE=SNUM/SDEN
          else
             Y2=1._wp/X2
             SNUM= R(1)+Y2*(R(2)+Y2*(R(3)+Y2*(R(4)+Y2*(R(5)+Y2*R(6)))))
             SDEN=S(1)+Y2*(S(2)+Y2*(S(3)+Y2*(S(4)+Y2*(S(5)+Y2*S(6)))))
             PERFCE= (0.56418958354775628694_wp+SNUM/SDEN*Y2)/X
          end if
       end if
       return
     end function PERFCE

     subroutine Coef_Coll
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       use Common_Params
       use Common_Arrays
       use Math_Constants
       use Phys_Constants
       implicit none
       real(wp) :: Ux,Uz,U,U2,U3,U5
       real(wp) :: Cee,Cei,RMeMe,RVeVi,RVeVe,RViVe
       real(wp) :: PsiArg,PhiArg,PhipArg,Arg,Arg2
       real(wp) :: ColAxe,ColAze,ColAxi,ColAzi
       real(wp) :: ColDxxe,ColDxze,ColDzze,ColDxxi,ColDxzi,ColDzzi
       !REAL(WP) :: PERF
       real(wp) :: Aux,ZZ
       integer :: l,m

       ZZ= 1  ! Version which works only for ions with unit charge.
       Cee= (2._wp*Pi)*Geff * log(3._wp/(8._wp*Pi)/sqrt(2._wp)/Geff)
       Cei= Cee*ZZ**2
       RMeMe= 1._wp
       RVeVi= sqrt(RTeTi*RMiMe)
       RVeVe= 1._wp
       RViVe= 1._wp/RVeVi

       select case(CollTerm)

       case("Yes")
          select case(CollTermForm)

          case("Complete")
             do m= 1,nuz
                Uz= VUz(m)
                do l= 1,nux
                   Ux= VUx(l)
                   U2= Ux**2+Uz**2
                   U= sqrt(U2)
                   U3= U2*U
                   U5= U3*U2
                   U3= U3+EpsMin
                   U5= U5+EpsMin
                   Aux= (tanh(U))**2

                   Arg= U*RVeVe
                   Arg2= Arg**2
                   PhiArg= PERF(Arg)
                   PhipArg= 2._wp*exp(-Arg2)/sqrt(Pi)
                   PsiArg= PhiArg-Arg*PhipArg
                   ColAxe= Cee*(2._wp/RMeMe*PsiArg*Ux/U3)
                   ColAze= Cee*(2._wp/RMeMe*PsiArg*Uz/U3)
                   ColDxxe= Cee*((U2*PhiArg-RVeVe**2*PsiArg/2._wp)*Uz*Uz/U5 &
                        + RVeVe**2*PsiArg*Ux*Ux/U5)
                   ColDxze= - Cee*((U2*PhiArg-RVeVe**2*PsiArg/2._wp)*Ux*Uz/U5 &
                        - RVeVe**2*PsiArg*Ux*Uz/U5)
                   ColDzze= Cee*((U2*PhiArg-RVeVe**2*PsiArg/2._wp)*Ux*Ux/U5 &
                        + RVeVe**2*PsiArg*Uz*Uz/U5)

                   Arg= U*RVeVi
                   Arg2= Arg**2
                   PhiArg= PERF(Arg)
                   PhipArg= 2._wp*exp(-Arg2)/sqrt(Pi)
                   PsiArg= PhiArg-Arg*PhipArg
                   ColAxi= Cei*(2._wp/RMiMe*PsiArg*Ux/U3)
                   ColAzi= Cei*(2._wp/RMiMe*PsiArg*Uz/U3)
                   ColDxxi= Cei*((U2*PhiArg-RViVe**2*PsiArg/2._wp)*Uz*Uz/U5 &
                        + RViVe**2*PsiArg*Ux*Ux/U5)
                   ColDxzi= - Cei*((U2*PhiArg-RViVe**2*PsiArg/2._wp)*Ux*Uz/U5 &
                        - RViVe**2*PsiArg*Ux*Uz/U5)
                   ColDzzi= Cei*((U2*PhiArg-RViVe**2*PsiArg/2._wp)*Ux*Ux/U5 &
                        + RViVe**2*PsiArg*Uz*Uz/U5)

                   ColAx(l,m)= Aux*(ColAxe+ColAxi)
                   ColAz(l,m)= Aux*(ColAze+ColAzi)
                   ColDxx(l,m)= Aux*(ColDxxe+ColDxxi)
                   ColDxz(l,m)= Aux*(ColDxze+ColDxzi)
                   ColDzz(l,m)= Aux*(ColDzze+ColDzzi)
                end do
             end do

          case("Expanded")
             do m= 1,nuz
                Uz= VUz(m)
                do l= 1,nux
                   Ux= VUx(l)
                   U2= Ux**2+Uz**2
                   U= sqrt(U2)
                   U3= U2*U
                   U5= U3*U2
                   U3= U3+EpsMin
                   U5= U5+EpsMin
                   Aux= (tanh(U))**2

                   Arg= U*RVeVe
                   Arg2= Arg**2
                   PhiArg= 1._wp
                   PsiArg= 1._wp
                   ColAxe= Cee*(2._wp/RMeMe*PsiArg*Ux/U3)
                   ColAze= Cee*(2._wp/RMeMe*PsiArg*Uz/U3)
                   ColDxxe= Cee*((U2*PhiArg-RVeVe**2*PsiArg/2._wp)*Uz*Uz/U5 &
                        + RVeVe**2*PsiArg*Ux*Ux/U5)
                   ColDxze= - Cee*((U2*PhiArg-RVeVe**2*PsiArg/2._wp)*Ux*Uz/U5 &
                        - RVeVe**2*PsiArg*Ux*Uz/U5)
                   ColDzze= Cee*((U2*PhiArg-RVeVe**2*PsiArg/2._wp)*Ux*Ux/U5 &
                        + RVeVe**2*PsiArg*Uz*Uz/U5)

                   Arg= U*RVeVi
                   Arg2= Arg**2
                   PhiArg= 1._wp
                   PsiArg= 1._wp
                   ColAxi= Cei*(2._wp/RMiMe*PsiArg*Ux/U3)
                   ColAzi= Cei*(2._wp/RMiMe*PsiArg*Uz/U3)
                   ColDxxi= Cei*((U2*PhiArg-RViVe**2*PsiArg/2._wp)*Uz*Uz/U5 &
                        + RViVe**2*PsiArg*Ux*Ux/U5)
                   ColDxzi= - Cei*((U2*PhiArg-RViVe**2*PsiArg/2._wp)*Ux*Uz/U5 &
                        - RViVe**2*PsiArg*Ux*Uz/U5)
                   ColDzzi= Cei*((U2*PhiArg-RViVe**2*PsiArg/2._wp)*Ux*Ux/U5 &
                        + RViVe**2*PsiArg*Uz*Uz/U5)

                   ColAx(l,m)= Aux*(ColAxe+ColAxi)
                   ColAz(l,m)= Aux*(ColAze+ColAzi)
                   ColDxx(l,m)= Aux*(ColDxxe+ColDxxi)
                   ColDxz(l,m)= Aux*(ColDxze+ColDxzi)
                   ColDzz(l,m)= Aux*(ColDzze+ColDzzi)
                end do
             end do

          case DEFAULT
             open(98,FILE='Warning_Coef_Coll.wt')
             write(98,*) ' CollTerm= ',CollTerm
             write(98,*) ' CollTerm must be (Yes) or (No ) !!'
             close(98)
             stop

          end select

       case("No ")
          do m= 1,nuz
             do l= 1,nux
                ColAx(l,m)= 0._wp
                ColAz(l,m)= 0._wp
                ColDxx(l,m)= 0._wp
                ColDxz(l,m)= 0._wp
                ColDzz(l,m)= 0._wp
             end do
          end do

       case DEFAULT
          open(98,FILE='Warning_Coef_Coll.wt')
          write(98,*) ' CollTerm= ',CollTerm
          write(98,*) ' CollTerm must be (Yes) or (No ) !!'
          close(98)
          stop

       end select

       return
     end subroutine Coef_Coll


     subroutine Coll_Damping
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       use Common_Params
       use Common_Arrays
       use Math_Constants
       use Phys_Constants
       use DQsimp_m
       use zetaFunc_m
       implicit none
       real(wp) :: Qx,Qz,Q2,Q
       real(wp) :: Zlq,Zsq,Mu
       real(wp) :: RTemp0,Res
       real(wp) :: Res1L,Res2L,Res1S,Res2S
       real(wp) :: DQ,Dmu
       real(wp) :: Aux
       real(wp) :: RTempf
       integer :: ires,m,sigma
       integer :: i,j,it,sigpm
       integer, parameter :: nqp= 1024,nmu=64 !64 !odd
       real(wp), dimension(nmu) :: Vmu,VintGcL,VintGcS
       real(wp) :: AuxSig
       real(wp), parameter :: Qi= 1.E-4_wp, Qf= 10.E0_wp

       m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
       RTemp0= VRTeTs(m)
       RTempf= VRTfTs(m)

       DQ= (Qf-Qi)/(nqcd-1)
       do i= 1,nqcd
          VQQ(i)= Qi+(i-1)*DQ
       end do
       VQQ(nqcd)= Qf
       
       Dmu= (1._wp-(-1._wp))/(nmu-1)
       do i= 1,nmu
          Vmu(i)= -1._wp+(i-1)*Dmu
       end do
       Vmu(nmu)= 1

       select case(Gcoll)
       case("Yes")
          open(1,file="GcollL") !! one column file
          ! Evaluation of the 1D expression:
          do i= 1,nqcd
             Q= VQQ(i)
             Q2= Q**2
             Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
             Zsq= Q*AA*sqrt(VRTeTs(m))/sqrt(1._wp+Q2/2._wp*VRTeTs(m)/VRNeNs(m))
             sigma= 1   ! Evaluated only for sigma=1, due to the symmetry
             ! ! Contribution of the background distribution (aproximated, U0= 0)
             Aux2_Gcoll(1)=i
             Aux1_Gcoll(1)= Rtemp0
             Aux1_Gcoll(2)= Q
             Aux1_Gcoll(3)= Q2
             Aux1_Gcoll(4)= Zlq
             Aux1_Gcoll(5)= Zsq   
             Aux2_Gcoll(2)= sigma
             ! do it= 1,nmu
             !    Aux2_Gcoll(3)= it
             !    Mu= Vmu(it)
             !    Aux1_Gcoll(6)= Mu
             !    call DQsimp(Aux_GcollL,1.E-4_wp,0.28E0_wp,Res1L)
             !    call DQsimpb(Aux_GcollL,0.28E0_wp,1.E30_wp,Res2L)
             !    VintGcL(it)= Res1L+Res2L
             ! end do
             ! call Simpson(Vmu,VintGcL,nmu,Res)
             ! GcollL1D(i)=-2.0_wp*Pi**(2.5)*VRNeNs(m)**4/VRTeTs(m)**3 &
             !      * Geff/Q2/Q*Res*(1._wp-RatioNf)/Zlq**3
             read(1,*) GcollL1D(i)
             !! For S waves !!
             call DQsimp(Aux_GcollS,1.E-4_wp,1._wp,Res1S)
             call DQsimpb(Aux_GcollS,1._wp,1.E30_wp,Res2S)
             Res= Res1S+Res2S
             GcollS1D(i)= -Zlq*Zsq*(16._wp*sqrt(Pi))/Q2 * (Q**3*AA/2._wp) &
                  * (VRNeNs(m)/VRTeTs(m))**(5./2.)*Geff * Res * (1._wp-RatioNf)
             !! GcollL1D(i)= -Zlq**2*(16._wp*sqrt(Pi))/(sqrt(VRTeTs(m)))**3/Q2 &
             !!      * (VRNeNs(m)/VRTeTs(m))**2*Geff * Res * (1._wp-RatioNf)            
             !! call DQsimp(Aux_GcollL,1.E-4_wp,4._wp,Res1L)
             !! call DQsimpb(Aux_GcollL,4._wp,1.E30_wp,Res2L)
             !! Res= Res1L+Res2L
             ! Evaluates the quasilinear damping:
             GqlL1D(i)= -sqrt(Pi)/(sqrt(VRTeTs(m)))**3*Zlq**2/Q**3 &
                  * exp(-Zlq**2/Q2/VRTeTs(m)) * (1._wp-RatioNf)
             GqlS1D(i)= -sqrt(Pi)*(Q**3*AA/2._wp)*Zlq*Zsq/Q**3 &
                  * (exp(-Zsq**2/Q2/VRTeTs(m))/(sqrt(VRTeTs(m)))**3 &
                  + exp(-Zsq**2/Q2/VRTeTs(m)*RTeTi*RMiMe)*sqrt(RMiMe*RTeTi/VRTeTs(m)) &
                  * RTeTi/VRTeTs(m) ) * (1._wp-RatioNf)
          end do     ! Q
          close(1)
          
       case("No ")
          do i= 1,nuz
             GcollL1D(i)= 0._wp
             GcollS1D(i)= 0._wp
             GqlL1D(i)= 0._wp
             GqlS1D(i)= 0._wp
          end do
       end select

       ! Evaluation of the 2D expression:
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
             call Aitp1d2(nqcd,VQQ,GcollL1D,Q,Aux,ires)
             GcollLm(i,j)= Aux
             call Aitp1d2(nqcd,VQQ,GqlL1D,Q,Aux,ires)
             GqlLp(i,j)= Aux
             GqlLm(i,j)= Aux
             call Aitp1d2(nqcd,VQQ,GcollS1D,Q,Aux,ires)
             GcollSp(i,j)= Aux
             GcollSm(i,j)= Aux
             call Aitp1d2(nqcd,VQQ,GqlS1D,Q,Aux,ires)
             GqlSp(i,j)= Aux
             GqlSm(i,j)= Aux
          end do    ! Qz
       end do     ! Qx

       return
     end subroutine Coll_Damping

     real(wp) function Aux_GcollL(Qp)
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       use Common_Params
       use Common_Arrays
       use Math_Constants
       use DQsimp_m
       use zetaFunc_m
       implicit none
       real(wp), intent(in) :: Qp
       real(wp) :: Q,Q2
       real(wp) :: Qp2
       real(wp) :: Zlq,Zsq,Rtemp
       real(wp) :: Qp2EpsQpSq
       real(wp) :: AuxA,AuxB,AuxA2,AuxB2
       real(wp) :: Aux,Aux0,Aux02
       integer :: m
       complex(wp) :: Qp2EpsQp
       complex(wp) :: Zeta
       m= 1
       
       Rtemp= Aux1_Gcoll(1)
       Q= Aux1_Gcoll(2)
       Q2= Aux1_Gcoll(3)
       Zlq= Aux1_Gcoll(4) 
       Zsq= Aux1_Gcoll(5)
 
       if(Q > 2.e-4_wp) then
          if(Qp > 2.e-4_wp) then
             Qp2= Qp**2
             Zeta= Zlq/Q/Qp/sqrt(VRTeTs(m))
             Qp2EpsQp= 1._wp+2._wp*VRNeNs(m)/VRTeTs(m)/Q2/Qp2*(1._wp+Zeta*Zfn(Zeta))
             !EpsQp2= ((Zlq**2-VRNeNs(m)-1.5_wp*VRTeTs(m)*Qp2)**2+4._wp*Pi*(VRNeNs(m))**5 &
             !  /(VRTeTs(m))**3/Qp**6*EXP(-2._wp*Zlq**2/Qp2))/(Zlq)**4
             Qp2EpsQpSq= (abs(Qp2EpsQp))**2
             ! AuxA= -2._wp*Q*Qp
             ! AuxB= 2._wp*(1._wp+RTeTi)*VRNeNs(m)/VRTeTs(m)+Q2+Qp2
             ! AuxA2= AuxA**2
             ! AuxB2= AuxB**2
             ! Aux= Qp*(Qp*(-AuxA2+2._wp*AuxB2)/(AuxB2-AuxA2)+(AuxB/2._wp/Q) &
             !      * log(1._wp+2._wp*AuxA/(AuxB-AuxA)))
             Aux= ((1._wp+Qp2*Qp2)/(1._wp-Qp2)**2+(1._wp+Qp2)/4._wp/Qp &
                  * log((1._wp-Qp)**2/(1._wp+Qp)**2))/Qp2
             Aux_GcollL= exp(-Zlq**2/Qp2/Q2/VRTeTs(m))*Aux/Qp2EpsQpSq
          else
             Aux_GcollL= 0.E0_wp
          end if
       else
          Aux_GcollL= 0.E0_wp
       end if
       return
     end function Aux_GcollL

     real(wp) function Aux_GcollS(Qp)
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       use Common_Params
       use Common_Arrays
       use Math_Constants
       use DQsimp_m
       use zetaFunc_m
       implicit none
       real(wp), intent(in) :: Qp
       real(wp) :: Q,Q2,Qp2
       real(wp) :: Zlq,Zsq,Rtemp
       real(wp) :: EpsQp2,Mu,Aux0
       real(wp) :: AuxA,AuxB,AuxA2,AuxB2,Aux
       integer :: sigma,m,i,it
       complex(wp) :: EpsQp
       complex(wp) :: Zetae,Zetai

       m= 1
       i= Aux2_Gcoll(1) 
       sigma= Aux2_Gcoll(2)
       Rtemp= Aux1_Gcoll(1)
       Q= Aux1_Gcoll(2)
       Q2= Aux1_Gcoll(3)
       Zlq= Aux1_Gcoll(4) 
       Zsq= Aux1_Gcoll(5) 

       if(Q > 2.e-4_wp) then
          if(Qp > 2.e-4_wp) then
             Qp2= Qp**2
             Zetae= sigma*Zsq/Qp/sqrt(VRTeTs(m))
             Zetai= sigma*Zsq/Qp/sqrt(VRTeTs(m))*sqrt(RTeTi*RMiMe)
             EpsQp= 1._wp+2._wp*VRNeNs(m)/VRTeTs(m)/Qp2 *(1._wp+Zetae*Zfn(Zetae)) &
                  + 2._wp*VRNeNs(m)/VRTeTs(m)*RTeTi/Qp2 *(1._wp+Zetai*Zfn(Zetai))
             EpsQp2= (CDABS(EpsQp))**2

             AuxA= -2._wp*Q*Qp
             AuxB= 2._wp*(1.e0_wp+RTeTi)*VRNeNs(m)/VRTeTs(m)+Q2+Qp2
             AuxA2= AuxA**2
             AuxB2= AuxB**2
             Aux= (4._wp/(AuxB2-AuxA2)+RTeTi*Qp/Q/Q2/Qp2*(-2._wp*AuxA*AuxB/(AuxB2-AuxA2) &
                  + log((AuxA+AuxB)/(AuxB-AuxA))))
             Aux_GcollS= (exp(-Zsq**2/Qp2/VRTeTs(m))/VRTeTs(m)/sqrt(VRTeTs(m)) &
                  + exp(-Zsq**2/Qp2/VRTeTs(m)*RTeTi*RMiMe)*RTeTi/VRTeTs(m) &
                  * sqrt(RTeTi/VRTeTs(m)*RMiMe)) * Aux/Qp**3/EpsQp2

             ! Aux_GcollS= (GAMMA(Kappae+1._wp)/GAMMA(Kappae-0.5_wp)/(SQRT(VRTeTs(m)*(Kappae-1.5_wp)))**3 &
             !      * (1._wp+Zsq**2/Qp2/(Kappae-1.5_wp)/VRTeTs(m))**(-(Kappae+1._wp)) &
             !      + (1._wp+Zsq**2/Qp2/(Kappai-1.5_wp)/VRTeTs(m)*RTeTi*RMiMe)**(-(Kappai+1._wp)) &
             !      * RTeTi/VRTeTs(m)/(Kappai-1.5_wp)*SQRT(RTeTi/VRTeTs(m)/(Kappai-1.5_wp)*RMiMe) &
             !      * GAMMA(Kappai+1._wp)/GAMMA(Kappai-0.5_wp))*Aux/Qp**3/EpsQp2

             ! Aux_GcollS= (GAMMA(Kappae+1._wp)/GAMMA(Kappae-0.5_wp)/(SQRT(VRTeTs(m)*(Kappae-1.5_wp)))**3 &
             !      * (1._wp+Zsq**2/Qp2/(Kappae-1.5_wp)/VRTeTs(m))**(-(Kappae+1._wp)) &
             !      + EXP(-Zsq**2/Qp2/VRTeTs(m)*RTeTi*RMiMe)*RTeTi/VRTeTs(m) &
             !      * SQRT(RTeTi/VRTeTs(m)*RMiMe) ) * Aux/Qp**3/EpsQp2
          else
             Aux_GcollS= 0.E0_wp
          end if
       else
          Aux_GcollS= 0.E0_wp
       end if
       return
     end function Aux_GcollS

     subroutine Bremsstrahlung
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       use Common_Params
       use Common_Arrays
       use Math_Constants
       use Phys_Constants
       !use zetaFunc_m
       use DQsimp_m
       implicit none
       real(wp) :: Qx,Qz,Q2,Q
       real(wp) :: DQ,Dmu
       real(wp) :: Mu,Zlq,Zsq
       real(wp) :: RTemp0,ResL,ResS
       real(wp) :: Res1L,Res1S,Res2L,Res2S
       real(wp) :: Aux
       real(wp) :: time
       integer :: i,j,it
       integer :: m,sigma,ires
       integer, parameter :: nqp= 1024,nmu= 64 !odd
       real(wp), parameter :: Qi= 1.E-4_wp, Qf= 10.E0_wp
       !real(wp), parameter :: Qi= 1.E0_wp, Qf= 1.E2_wp
       real(wp), dimension(nmu) :: Vmu,VintL,VintS
       real(wp), dimension(nqx,nqz) :: BremssL,BremssS


       m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
       RTemp0= VRTeTs(m)

       DQ= (Qf-Qi)/(nqcd-1)
       do i= 1,nqcd
          VQQ(i)= Qi+(i-1)*DQ
       end do
       VQQ(nqcd)= Qf

       Dmu= (1._wp-(-1._wp))/(nmu-1)
       do i= 1,nmu
          Vmu(i)= -1._wp+(i-1)*Dmu
       end do
       Vmu(nmu)= 1._wp

       select case(Bremss)
       case("Yes")
          open(1,file="BremL")  !! one column file
          ! Evaluation of the 1D expression:
          sigma=1
          do i= 1,nqcd
             Q= VQQ(i)
             Q2= Q**2
             Zlq= sqrt(VRNeNs(m))*sqrt(1._wp+1.5_wp*Q2*VRTeTs(m)/VRNeNs(m))
             Zsq= Q*AA*sqrt(VRTeTs(m))/sqrt(1._wp+Q2/2._wp*VRTeTs(m)/VRNeNs(m))
             Aux1_Bremss(1)= Rtemp0
             Aux1_Bremss(2)= Q
             Aux1_Bremss(3)= Q2
             Aux1_Bremss(4)= Zlq
             Aux1_Bremss(5)= Zsq
             do it= 1,nmu
                Mu= Vmu(it)               
                Aux1_Bremss(6)= Mu
                ! call DQsimp(Aux_BremL,1.E-4_wp,4._wp,Res1L)
                ! call DQsimpb(Aux_BremL,4._wp,Infinity,Res2L)
                ! VintL(it)= Res1L+Res2L
                call DQsimp(Aux_BremS,1.E-4_wp,4._wp,Res1s)
                call DQsimpb(Aux_BremS,4._wp,Infinity,Res2S)
                VintS(it)= Res1S+Res2S
             end do
             ! call Simpson(Vmu,VintL,nmu,ResL)
             call Simpson(Vmu,VintS,nmu,ResS)
             ! BremL1D(i)= 6.0_wp*Pi**1.5_wp/Zlq**2*VRNeNs(m)**4/VRTeTs(m)**2.5*Geff**2/Q2/Q2*ResL
             ! BremL1D(i)=384._wp*sqrt(pi)/Zlq**2*(1._wp-1._wp/RMiMe/VRTiTs(m))**2 &
             !      * VRNeNs(m)**4/VRTeTs(m)*Geff**2/Q2*ResL
             ! BremL1D(i)= 96._wp*SQRT(pi)*VRNeNs(m)**5/Rtemp0**4 &
             !      * (1._wp-1._wp/VRTiTs(m)**2)**2*Geff**2/Q2*ResL           
             read(1,*) BremL1D(i)
             BremS1D(i)= 96._wp*sqrt(pi)*(AA/2._wp*Q2*Q)*VRNeNs(m)**5/Rtemp0**4 &
                  * (1._wp-1._wp/VRTiTs(m)**2)**2*Geff**2/Q2*ResS
          end do
          close(1)
       case("No ")
          do i= 1,nuz
             BremL1D(i)= 0._wp
             BremS1D(i)= 0._wp
          end do
       end select

       ! Evaluation of the 2D expression:
       
       BremLm= 0._wp
       BremLp= 0._wp
       BremSm= 0._wp
       BremSp= 0._wp
       do i= 1,nqx
          Qx= VQx(i)
          do j= 1,nqz
             Qz= VQz(j)
             Q2= Qx**2+Qz**2
             Q= sqrt(Q2)
             if(Q<=5.E-3_wp) Q=5.E-3_wp
             call Locate(VQQ,nqcd,Q,ires)
             call Aitp1d2(nqcd,VQQ,BremL1D,Q,Aux,ires)
             BremLp(i,j)= Aux
             BremLm(i,j)= Aux
             BremssL(i,j)= Aux
             call Aitp1d2(nqcd,VQQ,BremS1D,Q,Aux,ires)
             BremSp(i,j)= Aux
             BremSm(i,j)= Aux
             BremssS(i,j)= Aux
             !WRITE(98,*)BremLp(i,j),BremSp(i,j)
          end do
       end do

       return
     end subroutine Bremsstrahlung

     real(wp) function Aux_BremL(Qp)
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       use Common_Params
       use Common_Arrays
       use Math_Constants
       use DQsimp_m
       use zetaFunc_m
       implicit none
       real(wp), intent(in) :: Qp
       real(wp) :: Zlq,Zsq,Rtemp0
       real(wp) :: Q,Q2,Mu
       real(wp) :: Qp2
       real(wp) :: Qp2EpsQpSq
       real(wp) :: Aux,Aux0,Aux02
       real(wp) :: AuxA,AuxB,AuxA2,AuxB2
       integer :: m
       complex(wp) :: Qp2EpsQp
       complex(wp) :: Zeta

       m= 1
       Rtemp0= Aux1_Bremss(1)
       Q= Aux1_Bremss(2)
       Q2= Aux1_Bremss(3)
       Zlq= Aux1_Bremss(4)
       Zsq= Aux1_Bremss(5)
       Mu= Aux1_Bremss(6)
       Aux_BremL= 0._wp

       ! if(Q > 2.e-4_wp) then
       !    if(Qp > 2.e-4_wp) then
       !       Qp2= Qp**2
       !       Zeta= sigma*Zlq/Q/Qp/sqrt(VRTeTs(m))
       !       Qp2EpsQp= 1._wp+2._wp*VRNeNs(m)/VRTeTs(m)/Q2/Qp2*(1._wp+Zeta*Zfn(Zeta))
       !       ! EpsQp2= ((Zlq**2-VRNeNs(m)-1.5_wp*VRTeTs(m)*Qp2)**2+4._wp*Pi*(VRNeNs(m))**5 &
       !       !   /(VRTeTs(m))**3/Qp**6*EXP(-2._wp*Zlq**2/Qp2))/(Zlq)**4
       !       Qp2EpsQpSq= (ABS(Qp2EpsQp))**2
       !       ! AuxA= -2._wp*Q*Qp
       !       ! ! AuxB= 2._wp*(1._wp+RTeTi)*VRNeNs(m)/VRTeTs(m)+Q2+Qp2
       !       ! AuxB= Q2+Qp2
       !       ! AuxA2= AuxA**2
       !       ! AuxB2= AuxB**2
       !       ! ! Aux0= AuxB-AuxA
       !       ! Aux= (Qp*(-AuxA2+2._wp*AuxB2)/(AuxB2-AuxA2)+(AuxB/2._wp/Q) &
       !       !      * log((AuxB+AuxA)/(AuxB-AuxA)))/Qp2
       !       ! ! Aux= (Qp*(-AuxA2+2._wp*AuxB2)/(AuxB2-AuxA2)+(AuxB/2._wp/Q) &
       !       ! !      * log(1._wp+2._wp*AuxA/(AuxB-AuxA)))
       !       Aux= ((1._wp+Qp2*Qp2)/(1._wp-Qp2)**2+(1._wp+Qp2)/4._wp/Qp &
       !            * log((1._wp-Qp)**2/(1._wp+Qp)**2))/Qp2/Qp
       !       Aux_BremL= exp(-Zlq**2/Qp2/Q2/VRTeTs(m))*Aux/Qp2EpsQpSq
       !    else
       !       Aux_BremL= 0.E0_wp
       !    end if
       ! else
       !    Aux_BremL= 0.E0_wp
       ! end if
       
       Aux0= 0._wp
       Aux0= Qp**4*(Q**2+Qp**2-2._wp*Q*Qp*Mu) &
            / (2._wp/Rtemp0*VRNeNs(m)*(1._wp+1._wp/VRTiTs(m))+Qp**2)**2 &
            / (2._wp/Rtemp0*VRNeNs(m)*(1._wp+1._wp/VRTiTs(m))+Q**2+Qp**2-2._wp*Q*Qp*Mu)**2
       Aux_BremL= Aux0*(sqrt(1._wp/(Rtemp0*(1._wp*Qp**2 &   ! b=e e a=e
            + 1._wp*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &    
            * exp(-Zlq**2/(Rtemp0*(1._wp*Qp**2 &
            + 1._wp*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            + sqrt(1._wp/(Rtemp0*(1._wp*Qp**2 &   ! b=e e a=i
            + VRTiTs(m)/RMiMe*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            * exp(-Zlq**2/(Rtemp0*(1._wp*Qp**2 &
            + VRTiTs(m)/RMiMe*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            + sqrt(1._wp/(Rtemp0*(VRTiTs(m)/RMiMe*Qp**2 &  !b=i e a=e
            + 1._wp*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            * exp(-Zlq**2/(Rtemp0*(VRTiTs(m)/RMiMe*Qp**2 &
            + 1._wp*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            + sqrt(1._wp/(Rtemp0*(VRTiTs(m)/RMiMe*Qp**2 &  !b=i e a=i
            + VRTiTs(m)/RMiMe*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            * exp(-Zlq**2/(Rtemp0*(VRTiTs(m)/RMiMe*Qp**2 &
            + VRTiTs(m)/RMiMe*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))))

       return
     end function Aux_BremL

     real(wp) function Aux_BremS(Qp)
       use,intrinsic :: iso_fortran_env, only: wp=>real64
       use Common_Params
       use Common_Arrays
       use Math_Constants
       use DQsimp_m
       implicit none
       real(wp), intent(in) :: Qp
       real(wp) :: Zsq,Zlq,Rtemp0
       real(wp) :: Q,Q2,Mu
       real(wp) :: Aux0
       integer :: m

       m =1
       Rtemp0= Aux1_Bremss(1)
       Q= Aux1_Bremss(2)
       Q2= Aux1_Bremss(3)
       Zlq= Aux1_Bremss(4)
       Zsq= Aux1_Bremss(5)
       Mu= Aux1_Bremss(6)
       
       Aux_BremS= 0._wp
       Aux0= 0._wp

       ! Aux0= Qp/(2._wp*Rtemp0*VRNeNs(m)*(1._wp+1._wp/VRTiTs(m))+Qp**2)**2 &
       !      * 1._wp/(2._wp*Rtemp0*VRNeNs(m)*(1._wp+1._wp/VRTiTs(m))+Q**2+Qp**2 &
       !      - 2._wp*Q*Qp*Mu)**2

       Aux0= Qp**2/(2._wp*Rtemp0*VRNeNs(m)*(1._wp+1._wp/VRTiTs(m))+Qp**2)**2 &
            * 1._wp/(2._wp*Rtemp0*VRNeNs(m)*(1._wp+1._wp/VRTiTs(m))+Q**2+Qp**2 &
            - 2._wp*Q*Qp*Mu)**2
       ! Aux0= Qp**4*(Q**2+Qp**2-2._wp*Q*Qp*Mu) &
       !      / (2._wp*Rtemp0*VRNeNs(m)*(1._wp+1._wp/VRTiTs(m))+Qp**2)**2 &
       !      / (2._wp*Rtemp0*VRNeNs(m)*(1._wp+1._wp/VRTiTs(m))+Q**2+Qp**2-2._wp*Q*Qp*Mu)**2
       Aux_BremS= Aux0*(sqrt(1._wp/(Rtemp0*(1._wp*Qp**2 & ! b=e e a=e 
            + 1._wp*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            * exp(-Zsq**2/(Rtemp0*(1._wp*Qp**2 &
            + 1._wp*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            + sqrt(1._wp/(Rtemp0*(1._wp*Qp**2 &   !b=e e a=i 
            + VRTiTs(m)/RMiMe*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            * exp(-Zsq**2/(Rtemp0*(1._wp*Qp**2 &
            + VRTiTs(m)/RMiMe*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            + sqrt(1._wp/(Rtemp0*(VRTiTs(m)/RMiMe*Qp**2 &   !b=i e a=e
            + 1._wp*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            * exp(-Zsq**2/(Rtemp0*(VRTiTs(m)/RMiMe*Qp**2 &
            + 1._wp*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            + sqrt(1._wp/(Rtemp0*(VRTiTs(m)/RMiMe*Qp**2 &   !b=i e a=i
            + VRTiTs(m)/RMiMe*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))) &
            * exp(-Zsq**2/(Rtemp0*(VRTiTs(m)/RMiMe*Qp**2 &
            + VRTiTs(m)/RMiMe*(Q**2+Qp**2-2._wp*Q*Qp*Mu)))))

       ! veja que quando temos a=e ou b=e, eu deixo o 1. multiplicando
       ! explicitamente, para facilitar a identificação

       return
     end function Aux_BremS

     ! subroutine Rebuild
     !   ! Developed by J. Pavan, April 2009.

     !   use Common_Params
     !   use Common_Arrays
     !   use Math_Constants
     !   use Phys_Constants
     !   implicit none
     !   real(wp) :: P0,P1,P2,P3,t
     !   real(wp), dimension(nqx,nqz) :: ILp1,ILm1
     !   integer :: i,k,j,nuz2,isqr,ksqr

     !   call Locate(VQx,nqx,QxSqr,isqr)
     !   call Locate(VQz,nqz,QzSqr,ksqr)
     !   !
     !   ! Uses Four-Points Bezier Method.
     !   !==============================
     !   ! ISp,ISm: z=0 axis.
     !   !
     !   ! DO i=1,nqx
     !   ! P0=ISm(i,4)
     !   ! P1=ISm(i,3)
     !   ! P2=ISp(i,3)
     !   ! P3=ISp(i,4)
     !   ! k=4
     !   ! t=0.0
     !   ! DO j=1,3
     !   ! k=k-1
     !   ! t=t+1./6.
     !   ! ISm(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! k=4
     !   ! t=1.0
     !   ! DO j=1,3
     !   ! k=k-1
     !   ! t=t-1./6.
     !   ! ISp(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! END DO
     !   ! !=========================
     !   ! ! ISp,ISm: x=0 axis. 
     !   ! !
     !   ! DO k=1,nqz
     !   ! P0=ISm(4,k)
     !   ! P1=ISm(3,k)
     !   ! P2=ISm(3,k)
     !   ! P3=ISm(4,k)
     !   ! i=4
     !   ! t=1.0
     !   ! DO j=1,3
     !   ! i=i-1
     !   ! t=t-1./6.
     !   ! ISm(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! P0=ISp(4,k)
     !   ! P1=ISp(3,k)
     !   ! P2=ISp(3,k)
     !   ! P3=ISp(4,k)
     !   ! i=4
     !   ! t=1.0
     !   ! DO j=1,3
     !   ! i=i-1
     !   ! t=t-1./6.
     !   ! ISp(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! END DO
     !   !
     !   ! IT waves are not smoothed out with this subroutine.
     !   !
     !   !==============================
     !   ! ILp,ILm: z=0 axis.
     !   !
     !   ! DO i=1,isqr-1	! Until the beginning of the 'squared region'; 
     !   !                 ! see 'SUBROUTINE Coef_Lwave'.
     !   ! P0=ILm(i,4)
     !   ! P1=ILm(i,3)
     !   ! P2=ILp(i,3)
     !   ! P3=ILp(i,4)
     !   ! k=4
     !   ! t=0.0
     !   ! DO j=1,3
     !   ! k=k-1
     !   ! t=t+1./6.
     !   ! ILm(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! k=4
     !   ! t=1.0
     !   ! DO j=1,3
     !   ! k=k-1
     !   ! t=t-1./6.
     !   ! ILp(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! END DO
     !   ! !=========================
     !   ! ! ILp,ILm: x=0 axis. 
     !   ! !
     !   ! DO k=1,nqz
     !   ! P0=ILm(4,k)
     !   ! P1=ILm(3,k)
     !   ! P2=ILm(3,k)
     !   ! P3=ILm(4,k)
     !   ! i=4
     !   ! t=1.0
     !   ! DO j=1,3
     !   ! i=i-1
     !   ! t=t-1./6.
     !   ! ILm(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! P0=ILp(4,k)
     !   ! P1=ILp(3,k)
     !   ! P2=ILp(3,k)
     !   ! P3=ILp(4,k)
     !   ! i=4
     !   ! t=1.0
     !   ! DO j=1,3
     !   ! i=i-1
     !   ! t=t-1./6.
     !   ! ILp(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! END DO
     !   ! !==============================
     !   ! ! L spectrum 'squared region'.
     !   ! ! Arbitrary weights.
     !   ! !
     !   ! ILp1=ILp
     !   ! DO i=isqr,nqx
     !   ! DO k=ksqr+1,1,-1
     !   ! ILp1(i,k)=(0.5*ILp1(i-1,k)+1.5_wp*ILp1(i,k+1))/2._wp
     !   ! END DO
     !   ! END DO
     !   ! !
     !   ! ILm1=ILm
     !   ! DO i=isqr,nqx
     !   ! DO k=ksqr+1,1,-1
     !   ! ILm1(i,k)=(0.5*ILm1(i-1,k)+1.5_wp*ILm1(i,k+1))/2._wp
     !   ! END DO
     !   ! END DO
     !   ! !
     !   ! ILp=ILp1
     !   ! ILm=ILm1
     !   ! !==============================
     !   ! ! ILp,ILm: z=0 axis, 'squared region'.
     !   ! !
     !   ! DO i=isqr,nqx
     !   ! P0=ILm(i,4)
     !   ! P1=ILm(i,3)
     !   ! P2=ILp(i,3)
     !   ! P3=ILp(i,4)
     !   ! k=4
     !   ! t=0.0
     !   ! DO j=1,3
     !   ! k=k-1
     !   ! t=t+1./6.
     !   ! ILm(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! k=4
     !   ! t=1.0
     !   ! DO j=1,3
     !   ! k=k-1
     !   ! t=t-1./6.
     !   ! ILp(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! END DO
     !   !=================================
     !   ! Fe: z=0 axis.
     !   !
     !   ! nuz2= (nuz-1)/2+1
     !   ! DO i=1,nux
     !   ! P0=Fe(i,nuz2-3)
     !   ! P1=Fe(i,nuz2-2)
     !   ! P2=Fe(i,nuz2+2)
     !   ! P3=Fe(i,nuz2+3)
     !   ! k=nuz2-3
     !   ! t=0.0
     !   ! DO j=1,5
     !   ! k=k+1
     !   ! t=t+1./6.
     !   ! Fe(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! END DO
     !   ! !
     !   ! ! Fe: x=0 axis.
     !   ! !
     !   ! DO k=1,nuz
     !   ! P0=Fe(4,k)
     !   ! P1=Fe(3,k)
     !   ! P2=Fe(3,k)
     !   ! P3=Fe(4,k)
     !   ! i=4
     !   ! t=0.0
     !   ! DO j=1,3
     !   ! i=i-1
     !   ! t=t+1./6.
     !   ! Fe(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
     !   ! END DO
     !   ! END DO
     !   !
     !   return
     ! end subroutine Rebuild

     ! !MODULE Zfn_MODULE
     ! !CONTAINS
     ! !
     ! !Função Z para plasma Maxwelliano inserido ao programa akfvenanisobikappapro1 em 29/08/2013

     ! !c
     ! !c
     ! !c--------------------First Z-function code-----------------------------
     ! !c
     ! complex*16 function zfn(z)
     !   !c
     !   !c     Evaluates the plasma dispersion function (Fried and Conte 
     !   !c     function) of complex argument with a relative error of 1.e-6_wp.
     !   !c
     !   !c     Algorithm: based closely on that described in Piero Barberio-
     !   !c                Corsetti 'Calculation of the Plasma Dispersion 
     !   !c                Function'. 
     !   !c
     !   !c     Precision: Double
     !   !c
     !   !c     Author: R. L. Mace, Plasma Physics Research Institute
     !   !c               University of Natal, Durban
     !   !c
     !   !c     Modificacoes no programa original:
     !   !c        1)Precisao default: tol= 1.0d-6, dlim= 4.0d+00.  Para aumentar
     !   !c          precisao basta colocar tol= 1.0d-14, dlim= 6.0d+00.
     !   !c
     !   implicit none
     !   complex*16 z
     !   !c
     !   !c     constants
     !   !c
     !   real*8 tol,zero,half,one,dlim,thrhlf,pid4
     !   parameter( tol=1.0d-14, zero=0.d00, half=0.5d00, one=1.0d00 )
     !   parameter( dlim=6.0d00, thrhlf=1.5d00, pid4=0.785398163397448d00 )
     !   complex*16 czero,chalf,cone,ctwo,irtpi,i2rtpi
     !   parameter( czero=(0.d00,0.d00), chalf=(0.5d00,0.d00) )
     !   parameter( cone=(1.d00,0.d00), ctwo=(2.d00,0.d00) )
     !   parameter( irtpi=(0.d00,1.772453850905516d00) )
     !   parameter( i2rtpi=(0.d00,3.544907701811032d00) )
     !   !c
     !   !c     local variables
     !   !c
     !   real*8 x,y,abx,aby,xymax,fn,cn,aslim,yasm
     !   complex*16 errz,an,anm1,bn,bnm1,anp1,bnp1
     !   complex*16 z2,zinv,aa,bb,sum,term,pterm
     !   !c
     !   x=dreal(z)
     !   y=dimag(z)
     !   abx=dabs(x)
     !   aby=dabs(y)
     !   if (aby.gt.abx) then
     !      xymax=aby
     !   else
     !      xymax=abx
     !   endif
     !   fn=zero
     !   !c
     !   !c     based on the magnitude of the real and imaginary parts of z, 
     !   !c     determine which of power series, continued fraction, or 
     !   !c     asymptotic forms to use
     !   !c
     !   if (aby.gt.one) then
     !      !c
     !      !c       **********************************
     !      !c       employ the continued fraction form
     !      !c       **********************************
     !      !c
     !      z2=half-z*z
     !      an=z
     !      anm1=czero
     !      bn=z2
     !      bnm1=cone
     !      xymax=one/xymax
     !      !c
     !      !c       compute the continued fraction
     !      !c
     !      zfn=an/bn
     !      errz=zfn-anm1/bnm1
     !      do while (dabs(dreal(errz)).gt.tol*dabs(dreal(zfn)) .or.&
     !           dabs(dimag(errz)).gt.tol*dabs(dimag(zfn)))
     !         fn=fn+one
     !         cn=xymax/fn
     !         aa=-fn*(fn-half)*cn
     !         bb=(z2+fn+fn)*cn
     !         anp1=bb*an+aa*anm1
     !         bnp1=bb*bn+aa*bnm1
     !         anm1=an*cn
     !         an=anp1
     !         bnm1=bn*cn
     !         bn=bnp1
     !         zfn=an/bn
     !         errz=zfn-anm1/bnm1
     !      end do
     !      !c
     !      !c        add the contribution from the pole if Im(z) .le. 0
     !      !c
     !      if (y.le.zero) then
     !         zfn=zfn+i2rtpi*cdexp(-z*z)
     !      end if
     !   else if (abx.gt.dlim)  then
     !      !c
     !      !c        ****************************
     !      !c        use the asmyptotic expansion
     !      !c        ****************************
     !      !c
     !      zinv=cone/z
     !      z2=chalf*zinv*zinv
     !      sum=cone
     !      term=cone
     !      aslim=x*x+y*y-one
     !      do while (&
     !           (dabs(dreal(term)).gt.tol*dabs(dreal(sum)) .or.&
     !           dabs(dimag(term)).gt.tol*dabs(dimag(sum))) .and.&
     !           fn.le.aslim &
     !           )
     !         fn=fn+one
     !         term=term*(fn+fn-one)*z2
     !         sum=sum+term
     !      end do
     !      zfn=-zinv*sum
     !      yasm=pid4/abx
     !      if (y.lt.-yasm) then
     !         zfn=zfn+i2rtpi*cdexp(-z*z)
     !      else if (y.le.yasm) then
     !         zfn=zfn+irtpi*cdexp(-z*z)
     !      end if
     !   else
     !      !c
     !      !c        *************************
     !      !c        use the power series form
     !      !c        *************************
     !      !c
     !      z2=z*z
     !      sum=cone
     !      term=cone
     !      pterm=cone
     !      do while (&
     !           dabs(dreal(term)).gt.tol*dabs(dreal(sum)) .or.&
     !           dabs(dimag(term)).gt.tol*dabs(dimag(sum))&
     !           )
     !         fn=fn+one
     !         pterm=pterm*z2/fn
     !         term=pterm/(fn+fn+one)
     !         sum=sum+term
     !      end do
     !      zfn=(irtpi-ctwo*z*sum)*cdexp(-z2)
     !   end if
     !   return
     ! end function zfn
     ! !
     ! !END MODULE Zfn_MODULE
     
   end module Sub_Prog


   !---------------------------------------------------------------------
   ! Main program:
   !---------------------------------------------------------------------

   program WT_LST
     use,intrinsic :: iso_fortran_env, only: wp=>real64
     use Common_Params
     use Common_Arrays
     use Math_Constants
     use Phys_Constants
     use Sub_Prog
     implicit none
     real(wp) :: DTau,Tau,Tau1,Tau2,TauAdd
     real(wp) :: Dqx,Dqz,Dux,Duz
     real(wp) :: Ewave0,Ewave,Epart0,Epart,EppEw0,EppEw
     real(wp) :: EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
     real(wp) :: Rn,Rp,Rw,Rs
     real(wp) :: Qx,Qz
     real(wp) :: Anorm0,Anorm
     integer :: Iflag,i,j
     integer :: it,Nitera 
     character(LEN=6) :: WaveType
     character(LEN=3) :: TimeEvol
     character(LEN=3) :: OneDfiles,TwoDfiles,Onecolumnfiles,ADcoefficients

     open(1,FILE='Start.wt')
     read(1,*) DTau
     read(1,*) TauAdd
     read(1,*) TimeEvol
     read(1,*) OneDfiles
     read(1,*) TwoDfiles
     read(1,*) Onecolumnfiles
     read(1,*) ADcoefficients
     close(1)

     open(2,FILE='Ratios.wt')
     open(3,FILE='Ewave.wt')
     open(1,FILE='Ini.wt')
     read(1,*) Iflag
     read(1,*) nux,nuz
     read(1,*) nqx,nqz
     read(1,*) nqx2,nqz2
     read(1,*) nrz
     read(1,*) Ulim
     read(1,*) Ucrit
     read(1,*) Qxi 
     read(1,*) Qxf 
     read(1,*) Qzi 
     read(1,*) Qzf 
     read(1,*) InitialLevel 
     read(1,*) Iw0 
     read(1,*) AsympT
     read(1,*) Kappae
     read(1,*) Kappai
     read(1,*) RatioNf 
     read(1,*) Uf 
     read(1,*) RTfTe
     read(1,*) RatioNb 
     read(1,*) Ub 
     read(1,*) RTbTe
     read(1,*) RTeTi 
     read(1,*) G
     read(1,*) Ve2C2
     read(1,*) Lemis
     read(1,*) LdecayLS
     read(1,*) LdecayLT
     read(1,*) LdecayST
     read(1,*) LdecayTT
     read(1,*) LscatLL
     read(1,*) LscatLT
     read(1,*) Semis
     read(1,*) SdecayLL
     read(1,*) SdecayLT
     read(1,*) Sscat
     read(1,*) TdecayLL
     read(1,*) TdecayLS
     read(1,*) TdecayTL
     read(1,*) TscatLT
     read(1,*) CollTerm
     read(1,*) CollTermForm
     read(1,*) Gcoll
     read(1,*) GcollEvol
     read(1,*) Bremss
     read(1,*) BremssEvol
     read(1,*) NewEff_Init
     read(1,*) SpontEmis
     read(1,*) ScatElSpo
     read(1,*) RenormFe
     read(1,*) DerivLn

     if(Qxi<Qmin .or. Qzi<Qmin) then
        open(98,FILE='Warning_Main_Qi.wt')
        write(98,*) ' Qxi= ',Qxi,'  Qzi= ',Qzi
        write(98,*) ' Qxi and Qzi can not be too small !!'
        close(98)
        if(Qxi<Qmin) then
           Qxi= Qmin
        else
        end if
        if(Qzi<Qmin) then
           Qzi= Qmin
        else
        end if
     else
     end if

     if(nqx*nqz>6561) then
        open(98,FILE='Warning_Main_Size_RK4_Derivs.wt')
        write(98,*) ' nqx= ',nqx,'  nqz= ',nqz,'  nqx*nqz= ',nqx*nqz
        write(98,*) ' nqx*nqz can not be larger than NMAX in RK4 and DERIVS !!'
        close(98)
        stop
     else
     end if
     if (nqz2>nqz) then
        nph= nqz2
     else
        nph= nqz 
     end if

     call Allocate_Arrays
     call Definitions
     call Space_Profiles ! Just to generate profiles, which may be useful for
     ! eventual extension to the case of inhomogeneous medium.

     select case(Iflag)

     case(0)
        close(1)
        call Init_Wave(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
             EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
        call Output("Fe0")
        ! call Output("Fi ")
        call Output("IL0")
        call Output("IS0")
        ! call Output("IT0")
        call Output("GcL")
        call Output("GqL")
        call Output("GcS")
        call Output("GqS")
        call Output("BrL")
        call Output("BrS")
        it= 0
        Rn= 1.     ! Anorm/Anorm0
        Rp= 1.     ! Epart/Epart0
        Rw= 1.     ! Ewave/Ewave0
        Rs= 1.     ! EppEw/EppEw0
        !  write(2,2001) it,Tau,Rn,Rp,Rw,Rs
2001    format(1x,i5,7(1x,e13.5e3))
        !  write(3,2001) it,Tau,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
! 500     format(1x,5(e13.5e3,1x))
!         do i=1,nqx
!            Qx=VQx(i)
!            write(11,*)'  '
!            write(12,*)'  '
!            do j=1,nqz
!               Qz=VQz(j)
!               write(11,500) Qx,Qz,GcollLp(i,j),GcollLm(i,j)
!               write(12,500) Qx,Qz,GcollL(i,j),GcollLp(i,j)-GcollLm(i,j)
!            end do
!         end do
        

     case(1)
        call Read_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
             Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
        it= 0
        !  write(2,2001) it,Tau,Rn,Rp,Rw,Rs
        !  write(3,2001) it,Tau,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
        call Coll_Damping
        call Bremsstrahlung

     case DEFAULT
        open(98,FILE='Warning_Main.wt')
        write(98,*) ' Iflag= ',Iflag
        write(98,*) ' Iflag must be 0 or 1 !!'
        close(98)
        stop

     end select

     select case(TimeEvol)

     case("Yes")
        ! Start time evolution:

        call Res_Cond
        call Coef_A
        call Coef_Coll

        !call Output("IL0")
        !call Output("IS0")

        Auxinterp= 1


        ! call Output("GcL")
        ! call Output("GqL")
        ! call Output("GcS")
        ! call Output("GqS")
        ! call Output("BrL")
        ! call Output("BrS")
        !stop

        Tau1= Tau
        Tau2= Tau1+TauAdd
        Nitera= nint((Tau2-Tau1)/DTau)

        !***     Time evolution:    ***
        do it= 1,Nitera
           Tau= Tau+DTau
           ! Evolution of ILp:
           WaveType= "Lwavep"
           call Evol_Iwave(DTau,Tau,WaveType)
           ! Evolution of ILm:
           WaveType= "Lwavem"
           call Evol_Iwave(DTau,Tau,WaveType)
           ! Evolution of ISp:
           WaveType= "Swavep"
           call Evol_Iwave(DTau,Tau,WaveType)
           ! Evolution of ISm:
           WaveType= "Swavem"
           call Evol_Iwave(DTau,Tau,WaveType)
           ! Evolution of ITp:
           WaveType= "Twavep"
           call Evol_Iwave(DTau,Tau,WaveType)
           ! Evolution of ITm:
           WaveType= "Twavem"
           call Evol_Iwave(DTau,Tau,WaveType)

           call Split(Dux,Duz,DTau)

           Auxinterp= 0
           call Fnorm(Anorm)
           call Energy(Epart,Ewave,EppEw)
           Rn= Anorm/Anorm0
           Rp= Epart/Epart0
           Rw= Ewave/Ewave0
           Rs= EppEw/EppEw0
           !write (7,*) Rn, Anorm, Anorm0
           write(2,2001) it,Tau,Rn,Rp,Rw,Rs
           call Energy2(EwaveL,EWaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
           write(3,2001) it,Tau,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
           open(98,FILE='status.wt')
           write(98,*) ' it= ',it,'   Tau= ',Tau
           close(98)
        end do

        close(2)
        close(3)
        close(99)

        call Save_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
             Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)

        Iflag= 1

     case("No ")

        ! Nothing to be done at this point. Proceed to the next option.

     case DEFAULT
        open(98,FILE='Warning_WT_LST.wt')
        write(98,*) ' TimeEvol= ',TimeEvol
        write(98,*) ' TimeEvol must be (Yes) or (No ) !! '
        close(98)
        stop

     end select

     if (Iflag==1) then
        !call Rebuild

        select case(TwoDfiles)
        case("Yes")
           call Output("Fe ")
           call Output("IL ")
           call Output("IS ")
           call Output("Feu")
           !call Output("IT ")
        case("No ")
        case DEFAULT
           open(98,FILE='Warning_WT_LST.wt')
           write(98,*) ' TwoDfiles= ',TwoDfiles
           write(98,*) ' TwoDfiles must be (Yes) or (No ) !! '
           close(98)
           stop
        end select
        select case(OneDfiles)
        case("Yes")
           call Output("Fe1")
           call Output("IL1")
           call Output("IS1")
           !call Output("IT1")
        case("No ")
        case DEFAULT
           open(98,FILE='Warning_WT_LST.wt')
           write(98,*) ' OneDfiles= ',OneDfiles
           write(98,*) ' OneDfiles must be (Yes) or (No ) !! '
           close(98)
           stop
        end select
        select case(Onecolumnfiles)
        case("Yes")
           call Output2("Ui ")
           call Output2("Qi ")
           call Output2("Q  ")
           !call Output2("ZTQ")
           call Output2("Fe ")
           call Output2("IL ")
           call Output2("IS ")
           !call Output2("IT ")
           !call Output2("ITQ")
           call Output2("Fe1")
           call Output2("IL1")
           call Output2("IS1")
           !call output2("IT1")
        case("No ")
        case DEFAULT
           open(98,FILE='Warning_WT_LST.wt')
           write(98,*) ' Onecolumnfiles= ',Onecolumnfiles
           write(98,*) ' Onecolumnfiles must be (Yes) or (No ) !! '
           close(98)
           stop
        end select
        select case(ADcoefficients)
        case("Yes")
           call Output("Ai ")
           call Output("Dij")
        case("No ")
        case DEFAULT
           open(98,FILE='Warning_WT_LST.wt')
           write(98,*) ' ADCoefficients= ',ADcoefficients
           write(98,*) ' ADcoefficients must be (Yes) or (No ) !! '
           close(98)
           stop
        end select
     else
        open(98,FILE='Warning_WT_LST.wt')
        write(98,*) ' Iflag= ',Iflag
        write(98,*) ' Iflag must be (1) in order to generate output files!! '
        close(98)
        stop
     end if

   end program WT_LST

