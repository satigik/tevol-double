module zetaFunc_m
  ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
  ! !   This is a re-written version of the F77 version, originally written   ! !
  ! !    by "R. L. Mace". Runs slightly faster than the original version.     ! ! 
  ! !                                                                         ! !
  ! !     Evaluates the plasma dispersion function (Fried and Conte           ! !
  ! !     function) of complex argument with a relative error of 1.e-6.       ! !
  ! !                                                                         ! !
  ! !     Algorithm: based closely on that described in Piero Barberio-       ! !
  ! !                Corsetti 'Calculation of the Plasma Dispersion           ! !  
  ! !                Function'.                                               ! ! 
  ! !                                                                         ! ! 
  ! !     Precision: Double                                                   ! ! 
  ! !                                                                         ! !
  ! !     Author: R. L. Mace, Plasma Physics Research Institute               ! !
  ! !               University of Natal, Durban                               ! !
  ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
contains
  complex*16 function zfn(z)
    implicit none
    complex*16 :: z
    !! Constants
    real*8, parameter :: tol=1.D-14,zero=0.D0,half=0.5D00,one=1.D0
    real*8, parameter :: dlim=6.0D0,thrhlf=1.5D0,pid4=0.785398163397448D0
    complex*16 :: czero, chalf
    parameter (czero=(0.D0,0.D0), chalf=(0.5D0,0.D0))
    complex*16 :: cone, ctwo
    parameter (cone=(1.D0,0.D0), ctwo=(2.D0,0.D0))
    complex*16 :: irtpi,i2rtpi
    parameter (irtpi=(0.D0,1.772453850905516D0))
    parameter (i2rtpi=(0.D0,3.544907701811032D0))
    !! Local variables
    real*8 :: x,y,abx,aby,xymax
    real*8 :: fn,cn,aslim,yasm
    complex*16 :: an,anm1,bn,bnm1,errz
    complex*16 :: anp1,bnp1,aa,bb
    complex*16 :: z2,zinv,summ,term,pterm

    x = dreal(z)
    y = dimag(z)
    abx = dabs(x)
    aby = dabs(y)
    if(aby > abx) then
       xymax = aby
    else
       xymax = abx
    end if
    fn = zero
    !!     based on the magnitude of the real and imaginary parts of z, 
    !!     determine which of power series, continued fraction, or 
    !!     asymptotic forms to use
    if(aby > one) then
       !! **********************************
       !! employ the continued fraction form
       !! **********************************
       z2 = half-z*z
       an = z
       anm1 = czero
       bn = z2
       bnm1 = cone
       xymax = one/xymax
       !!       compute the continued fraction
       zfn = an/bn
       errz = zfn-anm1/bnm1
       do while(dabs(dreal(errz)) > tol*dabs(dreal(zfn)) .or. &
            dabs(dimag(errz)) > tol*dabs(dimag(zfn)))
          fn = fn+one
          cn = xymax/fn
          aa = -fn*(fn-half)*cn
          bb = (z2+fn+fn)*cn
          anp1 = bb*an+aa*anm1
          bnp1 = bb*an+aa*bnm1
          anm1 = an*cn
          an = anp1
          bnm1 = bn*cn
          bn = bnp1
          zfn = an/bn
          errz = zfn-anm1/bnm1
       end do
       
       !!  add the contribution from the pole if Im(z) .le. 0

       if(y < zero) zfn = zfn+i2rtpi*cdexp(-z*z)
    else if(abx > dlim) then
       !!  ****************************
       !!  use the asmyptotic expansion
       !!  ****************************
       zinv = cone/z
       z2 = chalf*zinv*zinv
       summ = cone
       term = cone
       aslim = x*x+y*y-one
       do while((dabs(dreal(term)) > tol*dabs(dreal(summ)) .or. &
            dabs(dimag(term)) > tol*dabs(dimag(summ))) .and. &
            fn <= aslim)
          fn = fn+one
          term = term*(fn+fn-one)*z2
          summ = summ+term
       end do
       zfn = -zinv*summ
       yasm = pid4/abx
       if(y < -yasm) then
          zfn = zfn+i2rtpi*cdexp(-z*z)
       else if(y <= yasm) then
          zfn = zfn+irtpi*cdexp(-z*z)
       end if
    else
       !!  *************************
       !!  use the power series form
       !!  ************************* 
       z2 = z*z
       summ = cone
       term = cone
       pterm = cone
       do while(dabs(dreal(term)) > tol*dabs(dreal(summ)) .or. &
            dabs(dimag(term)) > tol*dabs(dimag(summ)))
          fn = fn+one
          pterm = pterm*z2/fn
          term = pterm/(fn+fn+one)
          summ = summ+term
       end do
       zfn = (irtpi-ctwo*z*summ)*cdexp(-z2)
    end if
    return
  end function zfn
end module zetaFunc_m
