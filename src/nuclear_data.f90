module nuclear_data
use parameters
use parsing_utilities
use timer

!permanent nuclear data
type xs_data
    real(kind=8), dimension(:), allocatable :: sigt, vsigf
    real(kind=8), dimension(:), allocatable :: chi
    real(kind=8), dimension(:,:), allocatable :: sig_up
    real(kind=8), dimension(:), allocatable :: sig_dn
    integer :: groups, maxup, maxscat
end type

real(kind=8), dimension(:), allocatable :: phi
real(kind=8), dimension(:), allocatable :: energy
type(xs_data), dimension(:), allocatable :: xs

! real(kind=8), dimension(groups) :: sigt=0.0, vsigf=0.0
! real(kind=8), dimension(groups) :: chi = 0.0
! real(kind=8), dimension(maxup,maxup) :: sig_up=0.0
! real(kind=8), dimension(maxscat) :: sig_dn=0.0
! real(kind=8), dimension(groups) :: phi

contains

!==============================================================================

subroutine allocate_xsdata(this, g, gg)
    implicit none
    type(xs_data), intent(inout) :: this
    integer, intent(in), optional :: g, gg
    
    if (present(g)) then
        this%groups = g
    endif
    if (present(gg)) then
        this%maxup = gg
    endif
    
    this%maxscat = this%groups*(this%groups+1)/2
    allocate(this%sigt(this%groups), this%vsigf(this%groups), this%chi(this%groups))
    allocate(this%sig_up(this%maxup,this%maxup), this%sig_dn(this%maxscat))
    
    this%sigt=0.0
    this%vsigf=0.0
    this%chi=0.0
    this%sig_up=0.0
    this%sig_dn=0.0

end subroutine

!==============================================================================

subroutine output_data_ascii(xsdat)
    implicit none
    type(xs_data), intent(in) :: xsdat
    integer :: i

    open(unit=51,file='data',status='replace',action='write')
    open(unit=52,file='sig_dn',status='replace',action='write')
    open(unit=53,file='sig_up',status='replace',action='write')
    open(unit=54,file='mg_flux',status='replace',action='write')


    151 format (1X, 4E15.8)    
    
    write(51,*) 'energy_f90    chi_f90    sigt_f90    vsigf_f90'
    do i=1,groups
        write(51,151) energy(i), xsdat%chi(i), xsdat%sigt(i), xsdat%vsigf(i)
    enddo

    do i=1,maxscat
        ! write(52,*) xsdat%sig_dn(i)
    enddo

    do i=1,maxup
        ! write(53,103) xsdat%sig_up(i,:)
    enddo
    
    do i=1,groups
        write(54,*) phi(i)
    enddo

    103 format (1X, 1000E16.5)

    close(unit=51)
    close(unit=52)
    close(unit=53)
    close(unit=54)
    
end subroutine output_data_ascii

!==============================================================================

subroutine output_data_binary()
    implicit none
    open(unit=51,file='data_binary',status='replace',action='write',form='unformatted')
    write(51) energy
    write(51) xs(1)%chi
    write(51) xs(1)%sigt
    write(51) xs(1)%vsigf 
    write(51) xs(1)%sig_dn
    write(51) xs(1)%sig_up    
    close(unit=51)
    
end subroutine output_data_binary

!==============================================================================

subroutine input_data_binary()
    implicit none
    
    allocate(xs(1))
    call allocate_xsdata(xs(1),groups,maxup)
    
    open(unit=51,file='data_binary',action='read',form='unformatted')
    read(51) energy
    read(51) xs(1)%chi
    read(51) xs(1)%sigt
    read(51) xs(1)%vsigf 
    read(51) xs(1)%sig_dn
    read(51) xs(1)%sig_up
    close(unit=51)
end subroutine input_data_binary

!==============================================================================

subroutine input_data_ascii()
    implicit none
    integer :: i
    
    open(unit=51,file='data',action='read')
    open(unit=52,file='sig_dn',action='read')
    open(unit=53,file='sig_up',action='read')
    
    read(51,*)
    do i=1,groups
        read(51,*) energy(i), xs(1)%chi(i), xs(1)%sigt(i), xs(1)%vsigf(i)
    enddo
    
    do i=1,maxscat
        read(52,*) xs(1)%sig_dn(i)
    enddo
    
    do i=1,maxup
        read(53,103) xs(1)%sig_up(i,:)
    enddo
    
    103 format (1X, 1000E16.5)
    
    close(unit=51)
    close(unit=52)
    close(unit=53)
    
end subroutine input_data_ascii

!==============================================================================

subroutine input_from_gendf()

    implicit none
    
    !temporary nuclear data
    real(kind=8), dimension(:), allocatable :: sigma_t, vsigma_f, chi0, flux
    real(kind=8), dimension(:,:), allocatable :: sigma_upsc
    real(kind=8), dimension(:), allocatable :: sigma_scat

    !other variables
    integer :: i,maxi,mat,maxmat,ifiss,j
    real(kind=8) :: numdens
    real(kind=8), dimension(2) :: sigf
    character(len=20) :: fname
    
    !start executable statements
    allocate(sigma_t(groups),vsigma_f(groups),chi0(groups),flux(groups))
    allocate(sigma_upsc(maxup,maxup),sigma_scat(maxscat))
    
    open(unit=11,file='inputs/material_input',action='read')
    read(11,'(I3)') maxmat
    allocate(xs(maxmat))
    do mat=1,maxmat
        call allocate_xsdata(xs(mat),groups,maxup)
        read(11,'(I3)') maxi
        do i=1,maxi
            read(11,'(I3,E13.0,A20)') ifiss, numdens, fname
            write(*,*) mat, fname, read_timer()
            if (i==1) then
                call readgendf('inputs/'//fname,energy=energy)
                do j=1,groups
                    flux(j)=1/energy(j)
                enddo
            endif
            if (ifiss==1) then

                call readgendf('inputs/'//fname,sigt=sigma_t, chi=chi0, &
                    & vsigf=vsigma_f,sigs=sigma_scat,sigup=sigma_upsc)
                
                ! if (i==1) then
                    ! write(*,*) fname, numdens
                    ! open(unit=32,file='dil_238',action='read')
                    ! do j=1,maxup
                        ! read(32,*)
                        ! read(32,*)
                        ! read(32,*)
                        ! read(32,*)
                    ! enddo
                    ! do j=groups-maxup,1,-1
                        ! read(32,'(11X,1E9.0)') sigma_t(j)
                        ! read(32,*)
                        ! read(32,*)
                        ! read(32,*)
                    ! enddo
                    ! read(32,*)
                    ! do j=groups,1,-1
                        ! read(32,'(2E10.0)') sigf
                        ! vsigma_f(j) = vsigma_f(j)*sigf(2)/sigf(1)
                    ! enddo
                    ! close(unit=32)
                                          
                ! elseif (i==2) then
                    ! write(*,*) fname, numdens
                    ! open(unit=32,file='dil_235',action='read')
                    ! do j=1,maxup
                        ! read(32,*)
                        ! read(32,*)
                        ! read(32,*)
                        ! read(32,*)
                    ! enddo
                    ! do j=groups-maxup,1,-1
                        ! read(32,'(11X,1E9.0)') sigma_t(j)
                        ! read(32,*)
                        ! read(32,*)
                        ! read(32,*)
                    ! enddo
                    ! read(32,*)
                    ! do j=groups,1,-1
                        ! read(32,'(2E10.0)') sigf
                        ! vsigma_f(j) = vsigma_f(j)*sigf(2)/sigf(1)
                    ! enddo
                    ! close(unit=32)
                ! endif
                
                
                xs(mat)%sigt = xs(mat)%sigt+numdens*sigma_t
                xs(mat)%vsigf = xs(mat)%vsigf+numdens*vsigma_f
                xs(mat)%sig_up = xs(mat)%sig_up+numdens*sigma_upsc
                xs(mat)%sig_dn = xs(mat)%sig_dn+numdens*sigma_scat
                xs(mat)%chi = xs(mat)%chi+numdens*dot_product(vsigma_f,flux)*chi0
            else
                call readgendf('inputs/'//fname,sigt=sigma_t, &
                    & sigs=sigma_scat,sigup=sigma_upsc)
                xs(mat)%sigt = xs(mat)%sigt+numdens*sigma_t
                xs(mat)%sig_up = xs(mat)%sig_up+numdens*sigma_upsc
                xs(mat)%sig_dn = xs(mat)%sig_dn+numdens*sigma_scat
            endif  
        enddo
        
        ! open(unit=12,file='fission',action='read')
        ! do i=1,groups
            ! read(12,*) xs(mat)%chi(i)
        ! enddo
        ! close(unit=12)
        
        if (sum(xs(mat)%chi)>0) then
            xs(mat)%chi=xs(mat)%chi/sum(xs(mat)%chi)  
        endif
    enddo
    
    deallocate(sigma_t,vsigma_f,chi0,sigma_upsc,sigma_scat)
    close(unit=11)
    

end subroutine input_from_gendf

!==============================================================================

subroutine grab_xs()

    implicit none
    
    !temporary nuclear data
    real(kind=8), dimension(:), allocatable :: sigma_t, vsigma_f, chi0
    real(kind=8), dimension(:,:), allocatable :: sigma_f, sigma_nxn, sigma_dis
    real(kind=8), dimension(:,:), allocatable :: sigma_upsc
    real(kind=8), dimension(:), allocatable :: sigma_scat

    !other variables
    integer :: i,maxi,mat,maxmat,ifiss,j
    real(kind=8) :: numdens
    character(len=20) :: fname
    
    !start executable statements
    allocate(sigma_t(groups),vsigma_f(groups),chi0(groups))
    allocate(sigma_upsc(maxup,maxup),sigma_scat(maxscat))
    allocate(sigma_f(groups,2),sigma_nxn(groups,2),sigma_dis(groups,2))
    sigma_f = 0.0
    sigma_nxn = 0.0
    sigma_dis = 0.0
    
    open(unit=11,file='inputs/energy_input',action='read')
    read(11,'(I3)') maxmat
    allocate(xs(1))
    
    call allocate_xsdata(xs(1),groups,maxup)
    read(11,'(I3)') maxi
    do i=1,maxi
        read(11,'(I3,E13.0,A20)') ifiss, numdens, fname
        if (i==1) then
            call readgendf('inputs/'//fname,energy=energy)
        endif
        if (ifiss==1) then
            
            call readgendf('inputs/'//fname,sigt=sigma_t, chi=chi0, &
                & vsigf=vsigma_f,sigs=sigma_scat,sigup=sigma_upsc, &
                sigdis=sigma_dis(:,1), signxn=sigma_nxn(:,1), sigf=sigma_f(:,1))
                        
            xs(1)%sigt = xs(1)%sigt+numdens*sigma_t
            xs(1)%vsigf = xs(1)%vsigf+numdens*vsigma_f
            xs(1)%sig_up = xs(1)%sig_up+numdens*sigma_upsc
            xs(1)%sig_dn = xs(1)%sig_dn+numdens*sigma_scat
            sigma_f(:,2) = sigma_f(:,2)+numdens*sigma_f(:,1)
            sigma_nxn(:,2) = sigma_nxn(:,2)+numdens*sigma_nxn(:,1)
            sigma_dis(:,2) = sigma_dis(:,2)+numdens*sigma_dis(:,1)
            if (i==1) then
                xs(1)%chi=xs(1)%chi+0.2828*chi0
            elseif (i==2) then
                xs(1)%chi=xs(1)%chi+107.1805*chi0
            endif
        else
            call readgendf('inputs/'//fname,sigt=sigma_t, &
                & sigs=sigma_scat,sigup=sigma_upsc,sigdis=sigma_dis(:,1))
            xs(1)%sigt = xs(1)%sigt+numdens*sigma_t
            xs(1)%sig_up = xs(1)%sig_up+numdens*sigma_upsc
            xs(1)%sig_dn = xs(1)%sig_dn+numdens*sigma_scat
            sigma_dis(:,2) = sigma_dis(:,2)+numdens*sigma_dis(:,1)
        endif  
    enddo
            
    if (sum(xs(1)%chi)>0) then
        xs(1)%chi=xs(1)%chi/sum(xs(1)%chi)  
    endif
    
    open(unit=51,file='data',status='replace',action='write')
    open(unit=52,file='sig_dn',status='replace',action='write')
    open(unit=53,file='sig_up',status='replace',action='write')


    151 format (1X, 7E15.8)    
    
    write(51,*) 'energy_f90  chi_f90  sigt_f90  vsigf_f90  sigf_f90  signxn_f90  sigdis_f90'
    do i=1,groups
        write(51,151) energy(i), xs(1)%chi(i), xs(1)%sigt(i), xs(1)%vsigf(i), &
            & sigma_f(i,2), sigma_nxn(i,2), sigma_dis(i,2)
    enddo

    do i=1,maxscat
        write(52,*) xs(1)%sig_dn(i)
    enddo

    do i=1,maxup
        write(53,103) xs(1)%sig_up(i,:)
    enddo
    
    103 format (1X, 1000E16.5)

    close(unit=51)
    close(unit=52)
    close(unit=53)
    
    
    deallocate(sigma_t,vsigma_f,chi0,sigma_upsc,sigma_scat)
    close(unit=11)



end subroutine grab_xs

!==============================================================================
    
end module nuclear_data