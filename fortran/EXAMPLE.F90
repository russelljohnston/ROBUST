  include 'ROBUST.F90'
  program example

    use robust

    implicit none
    integer, parameter :: ngal=7877
    integer            :: i,answer
    real(kind=4)       :: mt(ngal),am(ngal),mu(ngal),dmu
    real(kind=4)       :: bin,mag_min,mag_max,delta_mu,delta_am
    character(len=12)  :: filein
    character(len=23)  :: fileout


    if(iargc().eq.2) then
        call getarg(1,filein)
        call getarg(2,fileout)
    else
       stop 'YOU PROBABLY HAVENT GOT THE RIGHT NUMBER OF INPUT PARAMS FROM SCRIPT'
    end if


    print *,'Reading in test data'
    filein='testdata.txt'
    open(unit=14,file='../'//filein,status='old')
    do i = 1,ngal
       read(14,*,err=998,end=999)mt(i),am(i),mu(i)
998 enddo
999 close(14)


    print *,"TEST DATA READ IN SUCCESSFULLY"
    print *,'min/max dist mod in sample is .........',minval(mu(1:ngal)),maxval(mu(1:ngal))
    print *,'min/max abs mag in sample is ..........',minval(am(1:ngal)),maxval(am(1:ngal))
    print *,'min/max app mag in sample is ..........',minval(mt(1:ngal)),maxval(mt(1:ngal))


    ! Run a Tc and Tv assuming a faint limit only
    mag_min = minval(mt(1:ngal))
    mag_max = maxval(mt(1:ngal))
    bin=0.1


    print *,'--------------------------------------------'
111 print *,'Which ROBUST version would you like to run?'
    print *,'1 - Rauzy 2001 (faint app mag lim only) - includes JTH07 Tv estimator'
    print *,'2 - Johnston et al 2007  - faint and bright lim - Tc & Tv estimator'
    read *, answer
    if (answer.eq.1)then
       open(unit=90,file=trim(fileout),status='unknown')
       rewind(90)
       call tctv_R01(ngal,mag_min,mag_max,mu,mt,am,bin)
    elseif (answer.eq.2)then
       open(unit=90,file=trim(fileout),status='unknown')
       rewind(90)
       print *,'input delta_mu width (e.g. 0.5)'
       read *, dmu
       delta_mu=dmu
       delta_am=dmu
       call tctv_JTH07(ngal,mag_min,mag_max,mu,mt,am,bin,delta_am,delta_mu)
    elseif (answer.ne.1 .or. answer.ne.2)then
       print*, 'you must enter 1 or 2'
       goto 111
    endif

  end program example
