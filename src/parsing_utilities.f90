module parsing_utilities
use parameters

contains 

!==============================================================================

subroutine readgendf(fname,energy,sigt,sigel,chi,sigf,vsigf,signxn,sigdis,sigs,sigup)

    implicit none

    !inputs
    character(len=*), intent(in) :: fname 
 
    !outputs
    real(kind=8), dimension(0:groups), intent(out), optional :: energy
    real(kind=8), dimension(groups), intent(out), optional :: sigt, sigel
    real(kind=8), dimension(groups), intent(out), optional :: sigdis, signxn
    real(kind=8), dimension(groups), intent(out), optional :: chi, vsigf, sigf
    real(kind=8), dimension(maxscat), intent(out), optional :: sigs
    real(kind=8), dimension(maxup, maxup), intent(out), optional :: sigup
    real(kind=8), dimension(groups) :: sig_2, chi_d, fiss, vsigf_d, chi_numer

    !counter variables
    integer :: i,j,ii,jj,n,nj,njjs

    !variables for reading in each line
    integer :: mat, mf, mt, nsp
    real(kind=8), dimension(6) :: a

    !other temporary variables
    real(kind=8), dimension(:), allocatable :: sig
    real(kind=8) :: fissflux
    real(kind=8), dimension(groups) :: fissvec
    integer :: ming, maxg, numsig,maxjj

    !variables to store data
    integer :: matl, nl
    real(kind=8) :: temp

    !--------------------------------------------------------------------------
    !executable statments
    
    !open GENDF file, move past header row
    open(unit=8,file=fname,action='read')
    read(8,*)

    !read in line
    read(8,100) a,mat,mf,mt,nsp
    100 format (6E11.0,I4,I2,I3,I5)

    matl = int(a(1))

    read(8,100) a,mat,mf,mt,nsp
    temp = a(1)

    !prepare lower triangular matrix for scattering 
    !    (upscattering corrections later)
    allocate(sig(groups))
   
    !read in energy structure
    if (present(energy)) then
        do i=1,ceiling(a(5)/6)
            read(8,100) a,mat,mf,mt,nsp
            !first energy at third position of first line
            if (i==1) then
                do ii=0,3
                    energy(groups-ii) = a(ii+3)
                enddo
                j=4
            !six energies per line
            elseif (j+5<groups) then
                do ii=1,6
                    energy(groups-ii-j+1) = a(ii)
                enddo
                j=j+6
            !last line may not be full
            elseif (j<=groups) then
                do ii=1,1+groups-j
                    energy(groups-ii-j+1) = a(ii)
                enddo
                j=groups+1
            endif
        enddo
    endif

    !get total cross section
    if (present(sigt)) then
        call find(3,1,a,mat,mf,mt,nsp)
        call vectorxs(groups,sig,ming,maxg)
        sigt=0.0
        sigt(ming:maxg) = sig(ming:maxg)
    endif
    
    !get total elastic scattering cross section
    if (present(sigel) .or. present(sigup)) then
        call find(3,2,a,mat,mf,mt,nsp)
        call vectorxs(groups,sig,ming,maxg)
        if (present(sigup)) then
            sig_2 = 0.0
            sig_2(ming:maxg) = sig(ming:maxg)
        endif
        if (present(sigel)) then
            sigel=0.0
            sigel(ming:maxg) = sig(ming:maxg)
        endif
    endif

    if (present(signxn)) then
        signxn = 0.0
        call find(3,16,a,mat,mf,mt,nsp)
        call vectorxs(groups,sig,ming,maxg)
        signxn(ming:maxg) = signxn(ming:maxg)+sig(ming:maxg)
        
        call find(3,17,a,mat,mf,mt,nsp)
        call vectorxs(groups,sig,ming,maxg)
        signxn(ming:maxg) = signxn(ming:maxg)+2*sig(ming:maxg)
        
        call find(3,37,a,mat,mf,mt,nsp)
        call vectorxs(groups,sig,ming,maxg)
        signxn(ming:maxg) = signxn(ming:maxg)+3*sig(ming:maxg)
        
        rewind(unit=8)
        read(8,*)
    endif
    
    !get sigf
    if (present(vsigf) .or. present(sigf)) then
        call find(3,18,a,mat,mf,mt,nsp)
        call vectorxs(groups,sig,ming,maxg)

        fiss=0.0
        fiss(ming:maxg) = sig(ming:maxg)
        if (present(sigf)) then
            sigf = fiss
        endif
    endif
    
    !get disappearance cross section
    if (present(sigdis)) then
        sigdis=0.0
        call find(3,102,a,mat,mf,mt,nsp)
        find_dis: do
            if (mt>110) then
                exit find_dis
            endif
            if (mt/=0) then
                write(*,*) mt
                call vectorxs(groups,sig,ming,maxg)
                sigdis(ming:maxg) = sigdis(ming:maxg)+sig(ming:maxg)
            endif
            read(8,100) a,mat,mf,mt,nsp
        enddo find_dis        
    endif

    ! !upscattering correction
    ! if (present(sigup) .and. present(sigt)) then
        ! findup01: do
            ! if (mt>=220 .and. mt<=250) then
                ! exit findup01
            ! endif
            ! read(8,100) a,mat,mf,mt,nsp
        ! enddo findup01
        ! call vectorxs(groups,sig,ming,maxg)
        ! sigt(ming:maxg)=sigt(ming:maxg)-sig_2(ming:maxg)+sig(ming:maxg)
        ! if (present(sigel)) then
            ! sigel(ming:maxg) = sig(ming:maxg)
        ! endif
    ! endif
    
    !get nu and build vsigf
    if (present(vsigf)) then
        call find(3,452,a,mat,mf,mt,nsp)
        call vectorxs(groups,sig,ming,maxg)
        vsigf=0.0
        vsigf(ming:maxg) = (/(fiss(i)*sig(i),i=ming,maxg)/)
        
        call find(3,455,a,mat,mf,mt,nsp)
        call vectorxs(groups,sig,ming,maxg)
        vsigf_d=0.0
        vsigf_d(ming:maxg) = (/(fiss(i)*sig(i),i=ming,maxg)/)
    endif
    
    deallocate(sig)
    allocate(sig((groups+1)*groups/2))
    
    !delayed chi spectrum
    if (present(chi)) then
        chi_d=0.0
        call find(5,455,a,mat,mf,mt,nsp)
        read(8,100) a,mat,mf,mt,nsp
        ming=int(a(4))
        maxg=int(a(3))-ming
        read(8,100) a,mat,mf,mt,nsp
        do i=ming,maxg
            read(8,100) a,mat,mf,mt,nsp
            chi_d(groups-i+1)=sum(a)
        enddo
    
    endif
        
    !elastic scattering
    if (present(sigs)) then
        sigs=0.0
        call find(6,2,a,mat,mf,mt,nsp)
        call scatxs(groups,maxscat,sigs)
    endif
    
    !check if n,2n is in wrong place
    if (present(sigs) .and. matl/=1001) then
        read(8,100) a,mat,mf,mt,nsp
        inelastic01: do 
            if (mt>37 .or. mt==18 .or. mt==19) then
                exit inelastic01
            endif
            if (mt/=0) then
                call scatxs(groups,maxscat,sigs)
            endif
            read(8,100) a,mat,mf,mt,nsp
        enddo inelastic01
    endif

    !chi distribution, similar to energy read
    if (present(chi)) then
        chi = 0.0
        if (mt/=18) then
            call find(6,18,a,mat,mf,mt,nsp)
        endif
        read(8,100) a,mat,mf,mt,nsp
        numsig = int(a(5))
        maxg = groups-int(a(4))+1
        j=1
        do i=1,ceiling(a(5)/6)
            read(8,100) a,mat,mf,mt,nsp
            if (j+5<numsig) then
                do ii=1,6
                    chi(maxg-ii-j+2) = a(ii)
                enddo
                j=j+6
            elseif (j<=numsig) then
                do ii=1,1+numsig-j
                    chi(maxg-ii-j+2) = a(ii)
                enddo
                j=groups+1
            endif    
        enddo
        
        chi_numer = 0.0
        do i=1,groups
            read(8,100) a,mat,mf,mt,nsp
            if (int(a(6)) /= i) then
                write(*,*) 'chi read error', i
                stop
            endif
            
            if (a(3)==2) then
                read(8,100) a,mat,mf,mt,nsp
                ! fissflux=1
                fissflux=a(1)
                chi_numer=chi_numer+fissflux*a(2)*chi + vsigf_d(groups-i+1)*fissflux*chi_d
            
            else
                !write(*,*) a
                
                maxg = groups-int(a(4))+1
                njjs = int(a(5))-1
                ming = groups-int(a(3))+2

                do j=1,ceiling(a(5)/6)
                    read(8,100) a,mat,mf,mt,nsp
                    if (j==1) then
                    !first line only contains up to 5 values
                        ! fissflux=1
                        fissflux=a(1)
                        nj=0
                        if (njjs>=5) then
                            maxjj=5
                        else
                            maxjj=njjs
                        endif

                        do jj=1,maxjj
                            chi_numer(maxg-nj)=chi_numer(maxg-nj) + &
                                & fissflux*a(jj+1)+fissflux*vsigf_d(groups-i+1)*chi_d(maxg-nj)
                            nj=nj+1
                        enddo
                
                    elseif (nj+5<njjs) then
                    !six values per line after first
                        do jj=1,6
                            chi_numer(maxg-nj)=chi_numer(maxg-nj) + &
                                & fissflux*a(jj)+fissflux*vsigf_d(groups-i+1)*chi_d(maxg-nj)
                            nj=nj+1
                        enddo
                
                    elseif (nj<=njjs) then
                    !last line may not be full
                        do jj=1,njjs-nj
                            chi_numer(maxg-nj)=chi_numer(maxg-nj) + &
                                & fissflux*a(jj)+fissflux*vsigf_d(groups-i+1)*chi_d(maxg-nj)
                            nj=nj+1
                        enddo
                    endif
                enddo

                
            endif
            
        enddo
        chi = chi_numer/sum(chi_numer)
    endif

    ! inelastic scattering
    if (present(sigs) .and. matl/=1001 .and. matl/=2004) then
        !find first inelastic cross section
        call find(6,51,a,mat,mf,mt,nsp)

        !continue for all inelastic cross sections
        !    also include n,xn and miscellaneous x,n
        !    may need to add separate grab for n,xn later
        inelastic: do 
            if (mt>91 .or. mt==18) then
                exit inelastic
            endif
            if (mt/=0) then
                ! write(*,*) mf, mt
                call scatxs(groups,maxscat,sigs)
            endif
            read(8,100) a,mat,mf,mt,nsp
        enddo inelastic
    endif
   
    !new upscatter matrix treatment
    if (present(sigup)) then
        if (.not. present(sigs)) then
            write(*,*) 'No upscattering without downscattering!'
            stop
        endif
        findup: do
            if (mt>=220 .and. mt<=250) then
                exit findup
            endif
            read(8,100) a,mat,mf,mt,nsp
        enddo findup
        
        sigup=0.0
        call matxs(maxup,sigup,maxg)
 
        !add downscatter matrix to empty parts of upscatter matrix
        do i=1,maxup
            n=groups-maxup+i
            if (i<=maxup-maxg) then
                ii = n*(n-1)/2+1
                jj = 0
                do j=i,1,-1
                    sigup(i,j)=sigs(ii+jj)
                    jj=jj+1
                    if (i>maxup .or. i<1 .or. j>maxup .or. j<1) then
                        write(*,*) 'Index Error!', i, j
                    endif
                    if (ii+jj > maxscat+1) then
                        write(*,*) 'Index Error 2!', ii, jj, maxscat
                    endif
                enddo
            elseif (maxup-maxg>0) then
                ii = ii+groups-maxup+i
                jj = 0
                do j=maxup-maxg,1,-1
                    sigup(i,j)=sigs(ii+jj)
                    jj=jj+1
                    if (i>maxup .or. i<1 .or. j>maxup .or. j<1) then
                        write(*,*) 'Index Error!', i, j
                    endif
                    if (ii+jj > maxscat+1) then
                        write(*,*) 'Index Error 2!', ii, jj, maxscat
                    endif
                enddo
            endif
        enddo
        
        ! correct total cross section
        do i=groups-maxup+1,groups
            sigt(i) = sigt(i)-sig_2(i)+sum(sigup(:,i-(groups-maxup)))
        enddo
        
    endif

    close(unit=8)

end subroutine readgendf

!==============================================================================

subroutine find(imf,imt,a,mat,mf,mt,nsp)
!skip over unwanted data until next data block of interest is reached
    implicit none
    integer,intent(in) :: imf, imt
    integer,intent(out) :: mat, mf, mt, nsp
    real(kind=8),dimension(6),intent(out) :: a
    integer :: io
    logical :: eof = .false.
    
    findit: do
        read(8,100, iostat=io) a,mat,mf,mt,nsp
        if (io<0 .and. .not. eof) then
            eof = .true.
            rewind(unit=8)
        elseif (io<0 .and. eof) then
            write(*,*) imf,imt,'not found on GENDF tape'
            stop
        endif
        100 format (6E11.0,I4,I2,I3,I5)
        if (mf==imf .and. mt==imt) then
            exit findit
        endif
    enddo findit    
end subroutine find

!==============================================================================

subroutine vectorxs(groups,sigma,ming,maxg)
!grab a vector cross section
!return the cross section and the minimum and maximum groups it includes
!could add functionality to insert zeros for skipped groups, but I haven't
!    run across this situation in my test cases
    implicit none
    integer, intent(in) :: groups
    integer, intent(out) :: ming, maxg
    real(kind=8),dimension(groups),intent(out) :: sigma
    integer :: mat, mf, mt, nsp, i, j
    real(kind=8), dimension(6) :: a
  
    j=0
    i=2
    loop: do
        read(8,100) a,mat,mf,mt,nsp
        if (mt==0) then
            exit loop
        elseif (int(a(6))==groups .and. i/=2) then
            exit loop
        endif
        i = groups-int(a(6))+1
        if (j==0) then
            maxg=i
            j=1
        endif
        read(8,100) a,mat,mf,mt,nsp
        sigma(i) = a(2)
    enddo loop
    ming=i
    
    100 format (6E11.0,I4,I2,I3,I5) 
end subroutine vectorxs

!==============================================================================

subroutine scatxs(groups,maxscat,sigma)
!grab scattering cross section by adding components to given scattering matrix
!lower triangular scattering matrix, only includes in-group and down-scattering
!thermal groups should be replaced by upscatter matrix
    implicit none
    integer, intent(in) :: groups, maxscat
    real(kind=8),dimension(maxscat),intent(inout) :: sigma
    integer :: njjs, mingroup
    integer :: mat, mf, mt, nsp, i, ii, j, jj, nj, maxjj, fix
    real(kind=8), dimension(6) :: a
    
    loop: do
        read(8,100) a,mat,mf,mt,nsp
        if (mt==0) then
            exit loop
        endif

        i = int(a(6)) !initial energy group
        mingroup = int(a(4)) !minimum energy group scattered to
        njjs = int(a(5))-1 !number of groups to read in
        
        ii = maxscat-(groups-i)
        do nj = 1,mingroup-1 !correct for minimum group
            ii = ii-(groups-nj+1)
        enddo

        do j=1,ceiling(a(5)/6)
            read(8,100) a,mat,mf,mt,nsp

            if (j==1) then
            !first line only contains up to 5 values
                nj=0
                if (njjs>=5) then
                    maxjj=5
                else
                    maxjj=njjs
                endif

                do jj=1,maxjj
                    sigma(ii)=sigma(ii)+a(jj+1)
                    ii=ii-(groups-nj-mingroup+1) 
                    nj=nj+1
                enddo
        
            elseif (nj+5<njjs) then
            !six values per line after first
                do jj=1,6
                    sigma(ii)=sigma(ii)+a(jj)
                    ii=ii-(groups-nj-mingroup+1) 
                    nj=nj+1
                enddo
        
            elseif (nj<=njjs) then
            !last line may not be full
                do jj=1,njjs-nj
                    sigma(ii)=sigma(ii)+a(jj)
                    ii=ii-(groups-nj-mingroup+1) 
                    nj=nj+1
                enddo

            endif
    
        enddo
    enddo loop
    
    !advance to end of block if loop exited early
    !    I don't think this is needed with the current version of the code
    do while (mt/=0)
        read(8,100) a,mat,mf,mt,nsp
        write(*,*) 'hello there!'
    enddo
    
    100 format (6E11.0,I4,I2,I3,I5)
end subroutine scatxs

!==============================================================================

subroutine matxs(size,sigma,maxg)
!grab a matrix cross section and store in two dimensional array
!no compressed storage
    implicit none
    integer, intent(in) :: size
    real(kind=8), dimension(size,size), intent(inout) :: sigma
    integer, intent(out) :: maxg
    integer :: i, j, ii, jj
    integer :: group, max, njjs, maxjj, nj
    integer :: mat, mf, mt, nsp, olda6
    real(kind=8), dimension(6) :: a
    
    olda6=0
    loop: do
        read(8,100) a,mat,mf,mt,nsp
        !write(*,*) a,mat,mf,mt,nsp
        if (int(a(6))>size .or. int(a(6))-olda6>1) then
            exit loop
        endif
        olda6=int(a(6))
        group = size-int(a(6))+1
        max = size-int(a(4))+1
        njjs = int(a(5))-1
        maxg = int(a(6))
            
        do j=1,ceiling(a(5)/6)
            read(8,100) a,mat,mf,mt,nsp
            if (j==1) then
            !first line only contains up to 5 values
                nj=0
                if (njjs>=5) then
                    maxjj=5
                else
                    maxjj=njjs
                endif

                do jj=1,maxjj
                    sigma(max-nj,group)=a(jj+1)
                    nj=nj+1
                enddo
        
            elseif (nj+5<njjs) then
            !six values per line after first
                do jj=1,6
                    sigma(max-nj,group)=a(jj)
                    nj=nj+1
                enddo
        
            elseif (nj<=njjs) then
            !last line may not be full
                do jj=1,njjs-nj
                    sigma(max-nj,group)=a(jj)
                    nj=nj+1
                enddo

            endif

        enddo
    enddo loop
    
    100 format (6E11.0,I4,I2,I3,I5) 
end subroutine matxs

end module parsing_utilities






