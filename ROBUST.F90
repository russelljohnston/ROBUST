!-----------------------------------------------------------------------------------------------
! CREATED BY RUSSELL JOHNSTON (rwi.johnston at gmail.com)
!               and
!            LUIS TEODORO     (luis.f.teodoro@nasa.gov)
!               and 
!            MARTIN HENDRY    (martin.hendry@glasgow.ac.uk)
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
! BOTH SUBROUTINES USE sort3.f and indexx.f from Numerical Recipes. 
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

  SUBROUTINE tctv_R01(Ngal,mag_min,mag_max,mu,mt,am,bin)
    implicit none
    integer Ngal,mstar_max
    PARAMETER (mstar_max=200)
    integer m,igal,jgal,iwksp(Ngal),ii,tc_inc
    real(kind=4) :: t_v(mstar_max),zmax,am(Ngal),mu(Ngal),mt(Ngal),mag_max,bin
    real(kind=4) :: appml(mstar_max),mag_min,absmlim,t_c(mstar_max),maxmaglim
    real(kind=4) :: varsum,s1,s2,s3,s4,r(mstar_max,Ngal),n(mstar_max,Ngal),p(mstar_max,Ngal)
    real(kind=4) :: zeta(mstar_max,Ngal),var(mstar_max,Ngal)
    real(kind=4) :: wksp(Ngal),tau(mstar_max,Ngal)


    print *,'CALLING: SUBROUTINE TCTV_FAINT'
    print *,'-----------------------------------'
    print *, 'faint app mag limit  = ',mag_min
    print *, 'bright app mag limit = ',mag_max
    print *, 'Tc,Tv bin incrememnt = ',bin
    print *,'-----------------------------------'

    maxmaglim = mag_max + 1.0 !making sure m* for Tc,Tv go beyond the survey limit
    tc_inc = floor((maxmaglim-mag_min)/bin)
    if(tc_inc.gt.mstar_max)then
       print*, 'too many steps in tc_inc (',tc_inc,'). nominaly set to 200,  increase mstar_max'
       stop
    endif
    write(90,*)'# m*        Tc(m*)      Tv(m*)'

    print *, 'computing Tc first'
    call sort3(Ngal,mu,mt,am,wksp,iwksp)
    do m = 1,tc_inc
       appml(m) = mag_min + bin*(m-1)        
       t_c(m)   = 0.0        
       varsum = 0.0
       do igal = 1,Ngal     
          s1			= 0.0
          s2			= 0.0
          r(m,igal)	= 0.0
          n(m,igal)	= 0.0
          zeta(m,igal)= 0.0

          if(mt(igal) .gt. appml(m)) goto 556 !REMOVE ALL GALS > m*
          absmlim = appml(m)- mu(igal) 
          if(am(igal) .gt. absmlim) goto 556 
          do jgal = igal,1,-1

             if(mt(jgal) .gt. appml(m))   goto 555 
             if(am(jgal) .gt. absmlim)    goto 555      
             if(mu(jgal) .gt. mu(igal))   goto 555	
             if(am(jgal) .le. am(igal))s1 = s1 + 1.0             
             if(am(jgal) .gt. am(igal))s2 = s2 + 1.0                
555       enddo
          r(m,igal) = s1
          n(m,igal) = s1+s2          
          zeta(m,igal)  = r(m,igal)/(n(m,igal)+1.0)                         
          var(m,igal)   = (n(m,igal)-1.0)/(12.0*(n(m,igal)+1.0))
          if(var(m,igal).ne.0.0)varsum = varsum + var(m,igal)
          t_c(m) = t_c(m) + zeta(m,igal) - 0.50        
556    enddo
       t_c(m) = t_c(m)/sqrt(varsum)
       !     print*,m,appml(m),t_c(m)
111    format(i4,1x,f6.2,1x,f11.6)
    enddo

    print *,'Now computing Tv'
    print *,''
    print *,'# m*    Tc(m*)        Tv(m*)'
    call sort3(Ngal,am,mt,mu,wksp,iwksp)
    do m = 1,tc_inc
       appml(m) = mag_min + bin*(m-1) 

       t_v(m)   = 0.0
       varsum   = 0.0 

       do igal = 1,Ngal
          s3		= 0.0
          s4		= 0.0        
          r(m,igal)	= 0.0
          p(m,igal)	= 0.0
          tau(m,igal)	= 0.0

          if(mt(igal) .gt. appml(m)) goto 558 !REMOVE ALL GALS > m*                 
          zmax= appml(m)-am(igal)
          if(mu(igal) .gt. zmax) goto 558        

          do jgal = igal,1,-1
             if(mt(jgal) .gt. appml(m))  goto 557
             if(mu(jgal) .gt. zmax)      goto 557
             if(am(jgal) .gt. am(igal))  goto 557   
             if(mu(jgal) .le. mu(igal))s3 = s3 + 1.0 					
             if(mu(jgal) .gt. mu(igal))s4 = s4 + 1.0
557       enddo
          r(m,igal) = s3
          p(m,igal) = s3+s4 

          tau(m,igal)= r(m,igal)/(p(m,igal)+1.0)                 
          var(m,igal) = (p(m,igal)-1.0)/(12.0*(p(m,igal)+1.0))  
          if(var(m,igal) .ne. 0.0)varsum = varsum +var(m,igal)                          
          t_v(m) = t_v(m) + tau(m,igal) -0.50              
558    enddo
       t_v(m) = t_v(m)/sqrt(varsum)
       write(6,112)appml(m),t_c(m),t_v(m)
       write(90,112)appml(m),t_c(m),t_v(m)
112    format(f6.2,1x,f12.4,1x,f12.4)

    enddo
    close(90)

  END SUBROUTINE tctv_R01



  SUBROUTINE tctv_JTH07(Ngal,mag_min,mag_max,mu,mt,am,bin,delta_am,delta_mu)

    implicit none
    integer ngal,mstar_max
    PARAMETER (mstar_max=200)
    integer       :: m,igal,jgal,j,tc_inc,iwksp(ngal),ii
    real(kind=4)  :: mt(ngal),mu(ngal),am(ngal),zmax,zmin,nn(mstar_max,ngal)
    real(kind=4)  :: appml(mstar_max),appm_min,t_c(mstar_max),t_v(mstar_max)
    real(kind=4)  :: varsum,s1,s2,s3,s4,r(mstar_max,ngal),n(mstar_max,ngal),bin
    real(kind=4)  :: absmlim_f,tau(mstar_max,ngal),rr(mstar_max,ngal),delta_mu,mag_min
    real(kind=4)  :: zeta(mstar_max,ngal),var(mstar_max,ngal),varsumzeta,delta_am
    real(kind=4)  :: wksp(ngal),absmlim_b,varzeta(mstar_max,ngal),maxmaglim,mag_max

    print *,'---------------------------------'
    print *,'CALLING SUBROUTINE: tctv_bright'
    print *,'binsize  = ',bin
    print *,'delta_mu = ',delta_mu
    print *,'delta_am = ',delta_am
    print *,'---------------------------------'


    maxmaglim = mag_max + 1.0 !making sure m* for Tc,Tv go beyond the survey limit
    tc_inc = floor((maxmaglim-mag_min)/bin)
    if(tc_inc.gt.mstar_max)then
       print*, 'too many steps in tc_inc (',tc_inc,'). nominaly set to 200,  increase mstar_max'
       stop
    endif
    write(91,*)'# m*        Tc(m*)      Tv(m*)'
    print *, 'computing Tc first'
    call sort3(Ngal,mu,mt,am,wksp,iwksp)
    do m = 1,tc_inc
       appml(m) = mag_min + bin*(m-1)
       t_c(m)   = 0.0        
       varsum   = 0.0       
       varsumzeta=0.0
       do igal =1,Ngal 
          s1            = 0.0
          s2            = 0.0
          r(m,igal)     = 0.0
          n(m,igal)     = 0.0        
          zeta(m,igal)  = 0.0 

          if(mt(igal).gt.appml(m)) goto 556	!REMOVE ALL GALS > m*
          if(mt(igal).lt.mag_min)  goto 556
          absmlim_f = appml(m)- mu(igal) !DEFINE ABSOLUTE MAG FAINT LIM
          absmlim_b = mag_min - mu(igal)+delta_mu !DEFINE ABSOLUTE MAG BRIGHT LIM 
          if(am(igal).lt.absmlim_b) goto 556
          if(am(igal).gt.absmlim_f) goto 556  

          do jgal = igal,1,-1
             if(am(jgal)  .lt. absmlim_b)  goto 555
             if(am(jgal)  .gt. absmlim_f)  goto 555  
             if(mu(jgal)  .lt. mu(igal) - delta_mu) goto 601

             if(am(jgal) .le. am(igal))s1 = s1+1.0										
             if(am(jgal) .gt. am(igal))s2 = s2 + 1.0  

555	  enddo
601	  continue
          r(m,igal) = s1
          n(m,igal) = s1+s2  
          zeta(m,igal)  = r(m,igal)/(n(m,igal)+1.0)           
          var(m,igal)   = (n(m,igal)-1.0)/(12.0*(n(m,igal)+1.0))
          varzeta(m,igal)=((zeta(m,igal)-0.50)**2.0)

          if(var(m,igal).ne.0.0)varsum = varsum + var(m,igal)
          if(varzeta(m,igal).ne.0.0)varsumzeta = varsumzeta + varzeta(m,igal)
          t_c(m) = t_c(m) + zeta(m,igal) - 0.50   
556    enddo

       t_c(m) = t_c(m)/sqrt(varsum)
       !	   write(6,111)m,appml(m),t_c(m)
       ! 111	   format(i4,1x,f6.2,1x,f11.6)
    enddo

    print *,'Now computing Tv'
    print *,''
    print *,'# m*    Tc(m*)        Tv(m*)'
    call sort3(Ngal,am,mt,mu,wksp,iwksp)
    do m = 1,tc_inc
       appml(m) = mag_min + bin*(m-1)
       t_v(m)   = 0.0        
       varsum   = 0.0       

       do igal =1,Ngal 
          s3           = 0.0
          s4           = 0.0
          rr(m,igal)   = 0.0
          nn(m,igal)   = 0.0
          tau(m,igal)  = 0.0 

          if(mt(igal).gt.appml(m)) goto 559	!REMOVE ALL GALS > m*
          if(mt(igal).lt.mag_min)  goto 559
          zmax = appml(m)- am(igal)	!DEFINE ABSOLUTE MAG FAINT LIM
          zmin = mag_min - am(igal)+delta_am !DEFINE ABSOLUTE MAG BRIGHT LIM 
          if(mu(igal) .lt. zmin) goto 559
          if(mu(igal) .gt. zmax) goto 559 
          if(zmax .lt. zmin) goto 559 

          do jgal = igal,1,-1
             if(mu(jgal)  .lt. zmin)     goto 558
             if(mu(jgal)  .gt. zmax)     goto 558	
             if(am(jgal)  .gt. am(igal)) goto 558			
             if(am(jgal)  .lt. am(igal) - delta_am) goto 602							

             if(mu(jgal) .le. mu(igal))s3 = s3+1.0
             if(mu(jgal) .gt. mu(igal))s4 = s4+1.0

558	  enddo
602	  continue
          rr(m,igal) = s3
          nn(m,igal) = s3+s4  
          tau(m,igal)  = rr(m,igal)/(nn(m,igal)+1.0) 

          var(m,igal)   = (nn(m,igal)-1.0)/(12.0*(nn(m,igal)+1.0))
          if(var(m,igal).ne.0.0)varsum = varsum + var(m,igal)
          t_v(m) = t_v(m) + tau(m,igal) - 0.50   

559    enddo
       t_v(m) = t_v(m)/sqrt(varsum)
112    format(f6.2,1x,f12.4,1x,f12.4)
       write(6,112)appml(m),t_c(m),t_v(m)
       write(90,112)appml(m),t_c(m),t_v(m)
    enddo
    close(90)
  END SUBROUTINE tctv_JTH07

!#### NUMERICAL RECIPES SUBROUTINES ####!!

    SUBROUTINE indexx(n,arr,indx)
    implicit none
    INTEGER n,indx(n),M,NSTACK
    real :: arr(n)
    PARAMETER (M=7,NSTACK=50)
    INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
    real(kind=4) :: a
    do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END 


    SUBROUTINE sort3(n,ra,rb,rc,wksp,iwksp)
	implicit none
	INTEGER n,iwksp(n)
	 real ::ra(n),rb(n),rc(n),wksp(n)
!CU    USES indexx
	INTEGER j
	call indexx(n,ra,iwksp)
	do 11 j=1,n
		wksp(j)=ra(j)
11    continue
	do 12 j=1,n
		ra(j)=wksp(iwksp(j))
12    continue
	do 13 j=1,n
		wksp(j)=rb(j)
13    continue
	do 14 j=1,n
		rb(j)=wksp(iwksp(j))
14    continue
	do 15 j=1,n
		wksp(j)=rc(j)
15    continue
	do 16 j=1,n
		rc(j)=wksp(iwksp(j))
16    continue

      return
      END





end module ROBUST
