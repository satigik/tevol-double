module DQsimp_m
contains
  subroutine DQsimp(func,a,b,s)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    ! Version 22/Aug/2015
    real(wp) :: a,b,s
    real(wp) :: EPS= 1.E-8_wp
    integer :: j
    integer :: jmax= 10
    real(wp) :: os,ost,st
    !      EXTERNAL func
    interface 
       real(wp) function FUNC(X)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: X
       end function FUNC
    end interface
    ost=-1.e+30_wp
    os= -1.e+30_wp
    do j=1,JMAX
       call trapzd(func,a,b,st,j)
       s=(4._wp*st-ost)/3._wp
       if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or. &
               (s.eq.0._wp.and.os.eq.0._wp)) return
       end if
       os=s
       ost=st
    end do
    open(98,file="WARNING_qsimp")
    write(98,*) " Too many steps in qsimp! "
    close(98)
    !      STOP
  end subroutine DQsimp


  subroutine Trapzd(func,a,b,s,n)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    ! Version 22/Aug/2015
    integer :: n
    integer :: it,j
    real(wp) :: a,b,s
    real(wp) :: del,suma,tnm,x
    !      external func
    interface 
       real(wp) function FUNC(X)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: X
       end function FUNC
    end interface
    if (n.eq.1) then
       s=0.5_wp*(b-a)*(func(a)+func(b))
    else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5_wp*del
       suma=0._wp
       do j=1,it
          suma=suma+func(x)
          x=x+del
       end do
       s=0.5_wp*(s+(b-a)*suma/tnm)
    endif
    return
  end subroutine Trapzd

  subroutine DQsimpb(func,a,b,s)
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    ! Version 22/Aug/2015
    ! To be used if "a' is finite, positive, and "b" goes to infinity, or if
    ! "b" is finite, negative, and "a" goes to -infinity.
    real(wp) :: a,b,s
    real(wp) :: EPS= 1.E-8_wp
    integer :: j
    integer :: jmax= 10
    real(wp) :: os,ost,st
    !      EXTERNAL func
    interface 
       real(wp) function func(x)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: x
       end function func
    end interface
    ost=-1.e+30_wp
    os= -1.e+30_wp
    do j=1,JMAX
       call midinf(func,a,b,st,j)
       s=(4._wp*st-ost)/3._wp
       if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or. &
               (s.eq.0._wp.and.os.eq.0._wp)) return
       end if
       os=s
       ost=st
    end do
    open(98,file="WARNING_qsimpb")
    write(98,*) " Too many steps in qsimpb! "
    close(98)
    !      STOP
  end subroutine DQsimpb

  subroutine MIDINF(FUNK,AA,BB,S,N)
    !	Numerical Recipes, Cambridge, 1989, p. 118:
    !	This routine is an exact replacement for MIDPNT, i.e. return as S 
    !	the Nth stage of refinement of the integral of FUNK from AA to BB,
    !	except that the function is evaluated at evenly spaced points in 1/x
    !	rather than in x. This allows the upper limit BB to be as large and
    !	positive as the computer allows, or the lower limit AA to be as large
    !	and negative, but not both. AA and BB must have the same sign.
    ! Version 22/Aug/2015
    !! use precision_m
    ! use cudafor
    ! use openacc
    use,intrinsic:: iso_fortran_env, only: wp=>real64
    implicit none
    integer :: n
    integer :: it,j
    real(wp) :: aa,bb,s,a,b
    real(wp) :: del,summ,tnm,x,ddel,func
    !      EXTERNAL func
    interface 
       real(wp) function FUNK(X)
         use,intrinsic:: iso_fortran_env, only: wp=>real64
         real(wp), intent(in) :: X
       end function FUNK
    end interface
    FUNC(X)=FUNK(1._wp/X)/X**2   ! This is a statement function which effects
    ! the change of variable. 
    B=1._wp/AA 		! These two statements change the limits of
    ! integration accordingly.
    A=1._wp/BB 
    if (N.eq.1) then 		! From this point on, the routine is exactly
       ! identical to MIDPNT.
       S=(B-A)*FUNC(0.5_wp*(A+B)) 
       IT=1 
    else 
       TNM=IT 
       DEL=(B-A)/(3._wp*TNM) 
       DDEL=DEL+DEL 
       X=A+0.5_wp*DEL 
       SUMM=0._wp 
       do 11 J=1,IT
       !do J=1,IT   
          SUMM=SUMM+FUNC(X) 
          X=X+DDEL 
          SUMM=SUMM+FUNC(X) 
          X=X+DEL
       !end do
11        continue 
       S=(S+(B-A)*SUMM/TNM)/3._wp 
       IT=3*IT 
    end if
    return 
  end subroutine MIDINF
end module DQsimp_m
