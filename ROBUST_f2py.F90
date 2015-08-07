!-----------------------------------------------------------------------------------------------
! CREATED BY RUSSELL JOHNSTON (rwi.johnston at gmail.com)
!               and
!            LUIS TEODORO     (luis.f.teodoro@nasa.gov)
!               and 
!            MARTIN HENDRY    (martin.hendry@glasgow.ac.uk)
!##################################################################
! This is a modified version of ROBUST.F90 for use of f2py conversion
! to import in python.
!##################################################################
!
! THIS MODULE CONTAINS VARIOUS SUBROUTINES TO ESTIMATE MAGNITUDE COMPLETENESS FOR
! A MAGNITUDE-REDSHIFT SURVEY.
!
! tctv_faint():
! COMPUTES  (i) THE TRADITIONAL RAUZY Tc COMPLETENESS TEST (http://arxiv.org/abs/astro-ph/0010357 )
!                which samples the CLF
!           (ii) THE JOHNSTON, TEODORO & HENDRY (2007) [JTH07- http://arxiv.org/abs/astro-ph/0703040] 
!                VARIANT Tv statistic which samples the CDF
!
! tctv_bright():
! COMPUTES   (i) THE JTH07 Tc and Tv statisitcs taking into account a bright apparent
!                magnitude limit of a given survey
!
! BOTH SUBROUTINES USE SORT3.F90 based on  Numerical Recipes, Created by 
! http://www.science-softcon.de/autochem/box-dir/box_html/998.html
!
!
! BOTH SUBROUTINES REQUIRE AS INPUT:
!
! Ngal = Number of galaxies in sample
! mag_min = brightest galaxy in sample (apparent magnitude)
! mag_max = faintest galaxy in sample (apparent magnitude)
! mu      = distance modulus,   [real array, size=Ngal]
! am      = absolute magnitude, [real array, size=Ngal]
! mt      = apparent magnitude, [real array, size=Ngal]
!
! for tctv_bright, it requires a further 2 parameters
! delta_mu  = the width of the boxes for Tc
! delta_am  = the width of the boxes for Tv
! (requires a little trial and error for the choice of size. default should be delta_mu=delta_am=0.5)
!
!----------------------------------------------------------------------------------------------
module ROBUST
  implicit none
contains

   SUBROUTINE tctv_R01(t_c,t_v,appml,tc_inc,mag_min,mag_max,mur,mtr,amr,bin,Ngal)
    implicit none
    integer, intent(in) :: Ngal
    integer, intent(out) :: tc_inc
    integer mstar_max
    PARAMETER (mstar_max=200)
    integer m,igal,jgal
    real(kind=4), DIMENSION(ngal), intent(in)  :: mtr,mur,amr
    real(kind=4), intent(in) :: mag_max,mag_min,bin
    real(kind=4), intent(out) ::  t_v(mstar_max),t_c(mstar_max),appml(mstar_max)
    real(kind=4) :: zmax,am(Ngal),mu(Ngal),mt(Ngal)
    real(kind=4) :: absmlim,maxmaglim,n(mstar_max,Ngal),p(mstar_max,Ngal)
    real(kind=4) :: varsum,s1,s2,s3,s4,r(mstar_max,Ngal)
    real(kind=4) :: zeta(mstar_max,Ngal),var(mstar_max,Ngal)
    real(kind=4) :: tau(mstar_max,Ngal)

    am = amr
    mu = mur
    mt = mtr


    print *,'CALLING: SUBROUTINE TCTV_R01'
    print *,'-----------------------------------'
    print *, 'faint app mag limit of sample  = ',mag_min
    print *, 'bright app mag limit of sample = ',mag_max
    print *, 'Tc(m*),Tv(m*) bin incrememnt   = ',bin
    print *,'-----------------------------------'

    maxmaglim = mag_max + 1.0 !making sure m* for Tc,Tv go beyond the survey limit
    tc_inc = floor((maxmaglim-mag_min)/bin)
    if(tc_inc.gt.mstar_max)then
       print*, 'too many steps in tc_inc (',tc_inc,'). nominaly set to 200,  increase mstar_max'
       stop
    endif

    print *, 'computing Tc first'
    r	  = 0.0
    n	  = 0.0
    zeta  = 0.0
    appml = 0.0
    t_c   = 0.0
    call SORT3(Ngal,mu,mt,am) ! sort in terms of Mu.
    do m = 1,tc_inc
       appml(m) = mag_min + bin*(m-1)
       varsum = 0.0
       do igal = 1,Ngal
          s1			= 0.0
          s2			= 0.0

          if(mt(igal) .le. appml(m)) then
             absmlim = appml(m)- mu(igal)
             if(am(igal) .le. absmlim) then
                do jgal = igal,1,-1
                   if(mt(jgal) .le. appml(m)) then
                      if(am(jgal) .le. absmlim) then
                         if(mu(jgal) .le. mu(igal)) then
                            if(am(jgal) .le. am(igal))s1 = s1 + 1.0
                            if(am(jgal) .gt. am(igal))s2 = s2 + 1.0
                         endif
                      endif
                   endif
                enddo

                r(m,igal) = s1
                n(m,igal) = s1+s2          
                zeta(m,igal)  = r(m,igal)/(n(m,igal)+1.0)
                var(m,igal)   = (n(m,igal)-1.0)/(12.0*(n(m,igal)+1.0))
                if(var(m,igal).ne.0.0)varsum = varsum + var(m,igal)
                t_c(m) = t_c(m) + zeta(m,igal) - 0.50
             endif
          endif
       enddo
       t_c(m) = t_c(m)/sqrt(varsum)
    enddo

    print *,'Now computing Tv'
    print *,''
!    print *,'# m*    Tc(m*)        Tv(m*)'
    call SORT3(Ngal,am,mt,mu)

    appml = 0.0
    t_v   = 0.0
    r     = 0.0
    p     = 0.0
    tau   = 0.0
    do m = 1,tc_inc
       appml(m) = mag_min + bin*(m-1) 
       varsum   = 0.0

       do igal = 1,Ngal
          s3		= 0.0
          s4		= 0.0        

          if(mt(igal) .le. appml(m)) then
             zmax= appml(m)-am(igal)
             if(mu(igal) .le. zmax) then

                do jgal = igal,1,-1
                   if(mt(jgal) .le. appml(m)) then
                      if(mu(jgal) .le. zmax) then
                         if(am(jgal) .le. am(igal)) then
                            if(mu(jgal) .le. mu(igal))s3 = s3 + 1.0
                            if(mu(jgal) .gt. mu(igal))s4 = s4 + 1.0
                         endif
                      endif
                   endif
                enddo
                r(m,igal) = s3
                p(m,igal) = s3+s4 

                tau(m,igal)= r(m,igal)/(p(m,igal)+1.0)                 
                var(m,igal) = (p(m,igal)-1.0)/(12.0*(p(m,igal)+1.0))  
                if(var(m,igal) .ne. 0.0)varsum = varsum +var(m,igal)                          
                t_v(m) = t_v(m) + tau(m,igal) -0.50
             endif
          endif
       enddo
       t_v(m) = t_v(m)/sqrt(varsum)
    enddo

  END SUBROUTINE tctv_R01



  SUBROUTINE tctv_JTH07(t_c,t_v,appml,tc_inc,mag_min,mag_max,mur,mtr,amr,bin,delta_am,delta_mu,ngal)

    implicit none
    integer :: mstar_max
    integer, intent(in) :: ngal
    integer, intent(out) :: tc_inc
    PARAMETER (mstar_max=200)
    integer       :: m,igal,jgal

    real(kind=4), DIMENSION(ngal), intent(in)  :: mtr,mur,amr
    real(kind=4), intent(in) :: mag_max,mag_min,bin,delta_mu,delta_am
    real(kind=4), intent(out) ::  t_v(mstar_max),t_c(mstar_max),appml(mstar_max)
    real(kind=4)  :: zmax,zmin,nn(mstar_max,ngal),mt(ngal),mu(ngal),am(ngal)
    real(kind=4)  :: varsum,s1,s2,s3,s4,r(mstar_max,ngal),n(mstar_max,ngal)
    real(kind=4)  :: absmlim_f,tau(mstar_max,ngal),rr(mstar_max,ngal)
    real(kind=4)  :: zeta(mstar_max,ngal),var(mstar_max,ngal),varsumzeta
    real(kind=4)  :: absmlim_b,varzeta(mstar_max,ngal),maxmaglim

    print *,'---------------------------------'
    print *,'CALLING SUBROUTINE: TCTV_JTH07'
    print *,'binsize  = ',bin
    print *,'delta_mu = ',delta_mu
    print *,'delta_am = ',delta_am
    print *,'---------------------------------'

    am = amr
    mu = mur
    mt = mtr

    maxmaglim = mag_max + 1.0 !making sure m* for Tc,Tv go beyond the survey limit
    tc_inc = floor((maxmaglim-mag_min)/bin)
    if(tc_inc.gt.mstar_max)then
       print*, 'too many steps in tc_inc (',tc_inc,'). nominaly set to 200,  increase mstar_max'
       stop
    endif
    print *, 'computing Tc first'
    call SORT3(Ngal,mu,mt,am)

    r     = 0.0
    n     = 0.0
    zeta  = 0.0
    do m = 1,tc_inc
       appml(m) = mag_min + bin*(m-1)
       t_c(m)   = 0.0        
       varsum   = 0.0       
       varsumzeta=0.0
       do igal =1,Ngal 
          s1  = 0.0
          s2  = 0.0

          if(mt(igal).le.appml(m).and. mt(igal).ge.mag_min) then
             absmlim_f = appml(m)- mu(igal)
             absmlim_b = mag_min - mu(igal)+delta_mu
             if(am(igal).ge.absmlim_b .and. am(igal).le.absmlim_f) then

                do jgal = igal,1,-1
                   if(am(jgal)  .ge. absmlim_b)  then
                      if(am(jgal)  .le. absmlim_f)  then
                         if(mu(jgal)  .lt. mu(igal) - delta_mu) goto 601

                         if(am(jgal) .le. am(igal))s1 = s1+1.0
                         if(am(jgal) .gt. am(igal))s2 = s2 + 1.0
                      endif
                   endif

                enddo
601             continue
                r(m,igal) = s1
                n(m,igal) = s1+s2  
                zeta(m,igal)  = r(m,igal)/(n(m,igal)+1.0)           
                var(m,igal)   = (n(m,igal)-1.0)/(12.0*(n(m,igal)+1.0))
                varzeta(m,igal)=((zeta(m,igal)-0.50)**2.0)

                if(var(m,igal).ne.0.0)varsum = varsum + var(m,igal)
                if(varzeta(m,igal).ne.0.0)varsumzeta = varsumzeta + varzeta(m,igal)
                t_c(m) = t_c(m) + zeta(m,igal) - 0.50
             endif
          endif
       enddo

       t_c(m) = t_c(m)/sqrt(varsum)
    enddo

    print *,'Now computing Tv'
    print *,''
    !print *,'# m*    Tc(m*)        Tv(m*)'
    call SORT3(Ngal,am,mt,mu)
    rr  = 0.0
    nn  = 0.0
    tau = 0.0
    t_v   = 0.0   
    do m = 1,tc_inc
       appml(m) = mag_min + bin*(m-1)
       varsum   = 0.0
       do igal =1,Ngal 
          s3           = 0.0
          s4           = 0.0
          if(mt(igal).le.appml(m).and.mt(igal).ge.mag_min) then
             zmax = appml(m)- am(igal)
             zmin = mag_min - am(igal)+delta_am
             if(mu(igal) .ge. zmin .and. mu(igal) .le. zmax) then

                   do jgal = igal,1,-1
                      if(mu(jgal)  .ge. zmin)     then
                         if(mu(jgal)  .le. zmax)     then
                            if(am(jgal)  .le. am(igal)) then
                               if(am(jgal)  .lt. am(igal) - delta_am) goto 602

                               if(mu(jgal) .le. mu(igal))s3 = s3+1.0
                               if(mu(jgal) .gt. mu(igal))s4 = s4+1.0

                            endif
                         endif
                      endif

                   enddo
602                continue
                   rr(m,igal) = s3
                   nn(m,igal) = s3+s4  
                   tau(m,igal)  = rr(m,igal)/(nn(m,igal)+1.0) 

                   var(m,igal)   = (nn(m,igal)-1.0)/(12.0*(nn(m,igal)+1.0))
                   if(var(m,igal).ne.0.0)varsum = varsum + var(m,igal)
                   t_v(m) = t_v(m) + tau(m,igal) - 0.50
                !endif
             endif
          endif

       enddo
       t_v(m) = t_v(m)/sqrt(varsum)
    enddo
  END SUBROUTINE tctv_JTH07


  SUBROUTINE SORT3(N,Ra,Rb,Rc)
    IMPLICIT NONE
    INTEGER :: N
    REAL , DIMENSION(N) :: Ra , Rb , Rc
    INTENT (IN) N
    INTENT (INOUT) Ra , Rb , Rc
    INTEGER :: i , ir , j , l
    REAL :: rra , rrb , rrc
    !
    !-----------------------------------------------------------------------
    !
    !     argrument list.
    !
    !     name      type              description.
    !     n         integer           length of arrays ra, rb, rc & rd.
    !     (unchanged on exit).
    !
    !     ra        array of real     array of length n to be sorted into
    !     ascending numerical order.
    !     (contains result on exit).
    !
    !     rb        array of real     array of length n to be sorted into
    !     an order corresponding to that of ra.
    !     (contains result on exit).
    !
    !     rc        array of real     array of length n to be sorted into
    !     an order corresponding to that of ra.
    !     (contains result on exit).
    !
    !-----------------------------------------------------------------------
    !
    !     looks at array ra and sorts it into ascending order, as well as
    !     making the corresponding rearrangements to rb, rc, & rd. the
    !     heapsort algorithm is used.
    !
    !     based on numerical recipes, the art of scientific computing,
    !     section 8.2, by press, flannery, teukolsky & vetterling,
    !     cambridge university press, 1987.
    !
    !     modified by: david lary
    !     ----------
    !
    !     date started : 28/9/1991
    !
    !     last modified: 28/9/1991
    !
    !-----------------------------------------------------------------------
    !
    l = N/2 + 1
    ir = N
    DO WHILE ( .TRUE. )
       IF ( l>1 ) THEN
          l = l - 1
          rra = Ra(l)
          rrb = Rb(l)
          rrc = Rc(l)
       ELSE
          rra = Ra(ir)
          rrb = Rb(ir)
          rrc = Rc(ir)
          Ra(ir) = Ra(1)
          Rb(ir) = Rb(1)
          Rc(ir) = Rc(1)
          ir = ir - 1
          IF ( ir<=1 ) THEN
             Ra(1) = rra
             Rb(1) = rrb
             Rc(1) = rrc
             EXIT
          ENDIF
       ENDIF
       i = l
       j = l + l
       DO WHILE ( .TRUE. )
          IF ( .NOT..TRUE. ) THEN
             RETURN
          ELSEIF ( j<=ir ) THEN
             IF ( j<ir ) THEN
                IF ( Ra(j)<Ra(j+1) ) j = j + 1
             ENDIF
             IF ( rra<Ra(j) ) THEN
                Ra(i) = Ra(j)
                Rb(i) = Rb(j)
                Rc(i) = Rc(j)
                i = j
                j = j + j
             ELSE
                j = ir + 1
             ENDIF
             CYCLE
          ENDIF
          Ra(i) = rra
          Rb(i) = rrb
          Rc(i) = rrc
          GOTO 100
       ENDDO
       EXIT
100 ENDDO
  END SUBROUTINE SORT3




end module ROBUST
