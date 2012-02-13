module coarse_group_type
use nuclear_data
use DCT_type

type coarse_group
    integer :: N, mingroup, index, up_index
    type(discreteLegendre) :: DLP
    real(kind=8), dimension(:), allocatable :: scal_flux            ! (region)
    real(kind=8), dimension(:,:,:), allocatable :: ang_flux         ! (region, moment, angle)
    real(kind=8), dimension(:), allocatable :: vsigf, sigt0         ! (region)
    real(kind=8), dimension(:,:,:), allocatable :: delta            ! (region, moment, angle)
    real(kind=8), dimension(:,:), allocatable :: chi                ! (region, moment)
    real(kind=8), dimension(:,:,:), allocatable :: sig_dn, sig_up   ! (region, moment, origin)
end type coarse_group

type(coarse_group), dimension(:), allocatable :: CG                 ! (group)

contains

!==============================================================================

subroutine input_coarse_groups()

    implicit none
    integer :: i, j
    integer :: mingroup
    logical :: is_unique

    open(unit=17,file='inputs/energy_input',action='read')
    read(17,*) groups, coarse_groups, upsc_groups
    maxup=0
    maxscat=groups*(groups+1)/2

    ! make a set of DLPs for each coarse group
    allocate(CG(coarse_groups))
    mingroup=1
    do i=1,coarse_groups
        is_unique = .true.
        read(17,*) CG(i)%N
        CG(i)%mingroup=mingroup
        mingroup=mingroup+CG(i)%N
        
        ! check if same length DLPs have already been generated
        ! if so, point to DLP rather than generate redundancy
        unique: do j=1,i-1
            if (CG(i)%N == CG(j)%N) then
                CG(i)%DLP%values => CG(j)%DLP%values
                CG(i)%DLP%a => CG(j)%DLP%a
                is_unique = .false.
                exit unique
            endif
        enddo unique
        if (is_unique) then
            call generate_DLPs(CG(i)%DLP,CG(i)%N)
        endif
        
        ! label scattering indices for groups
        ! build maxup value from upscatter coarse groups
        CG(i)%index=i
        if (i>coarse_groups-upsc_groups) then
            CG(i)%up_index=i-(coarse_groups-upsc_groups)
            maxup=maxup+CG(i)%N
        else
            CG(i)%up_index=0
        endif
    enddo

    ! error if inconsistent group mapping
    if (sum(CG%N) /= groups) then
        write(*,*) 'Invalid coarse group mapping.'
        stop
    endif
    
end subroutine input_coarse_groups

!==============================================================================

subroutine flux_moments(cgroup, flux, aflux)
    implicit none

    type(coarse_group), intent(inout) :: cgroup
    real(kind=8), dimension(:,:), intent(in) :: flux            ! (region, group)
    real(kind=8), dimension(:,:,:), intent(in) :: aflux         ! (region, group, angle)
    
    integer :: i, K, N, angle, region

    N = cgroup%N-1
    do region=1,size(aflux,1)
    
        ! condense coarse group scalar flux
        if (region < size(aflux,1)) then
            cgroup%scal_flux(region)=0.0
            do K=0,N
                cgroup%scal_flux(region) = cgroup%scal_flux(region)+flux(region,cgroup%mingroup+K)
            enddo
        endif
        
        ! expand angular flux
        cgroup%ang_flux(region,:,:)=0.0
        do angle=1,size(aflux,3)
            do i=0,N
                do K=0,N
                    cgroup%ang_flux(region,i,angle) = cgroup%ang_flux(region,i,angle) + &
                        & cgroup%DLP%values(i,K)*aflux(region,cgroup%mingroup+K,angle)
                enddo    
            enddo
        enddo
        
    enddo

end subroutine flux_moments

!==============================================================================

subroutine condense_flux_moments(cgroup, aflux)
    implicit none
    
    type(coarse_group), intent(in) :: cgroup
    real(kind=8), dimension(:,:,:), intent(out) :: aflux        ! (region, group, angle)
    
    integer :: i,K,N,angle,region
    
    N=cgroup%N-1
    
    do region=1,size(aflux,1)
    
        aflux(region,cgroup%mingroup:cgroup%mingroup+N,:)=0.0
        do angle=1,size(aflux,3)
            do K=0,N
                do i=0,N
                    aflux(region,cgroup%mingroup+K,angle) = aflux(region,cgroup%mingroup+K,angle) + &
                        & cgroup%DLP%a(i)*cgroup%DLP%values(i,K)*cgroup%ang_flux(region,i,angle)
                enddo
            enddo   
        enddo
        
    enddo
        
end subroutine condense_flux_moments

!==============================================================================

subroutine condense_flux_moments2(cgroup, dgm_aflux, aflux)
    implicit none
    
    type(coarse_group), intent(in) :: cgroup
    real(kind=8), dimension(:,:,:,:), intent(in) :: dgm_aflux    ! (region, group, moment, angle)
    real(kind=8), dimension(:,:,:), intent(out) :: aflux        ! (region, group, angle)
    
    integer :: i,K,N,angle,region
    
    N=cgroup%N-1
    
    do region=1,size(aflux,1)
    
        aflux(region,cgroup%mingroup:cgroup%mingroup+N,:)=0.0
        do angle=1,size(aflux,3)
            do K=0,N
                do i=0,N
                    aflux(region,cgroup%mingroup+K,angle) = aflux(region,cgroup%mingroup+K,angle) + &
                        & cgroup%DLP%a(i)*cgroup%DLP%values(i,K)*dgm_aflux(region,cgroup%index,i,angle)
                enddo
            enddo   
        enddo
        
    enddo
        
end subroutine condense_flux_moments2

!==============================================================================

subroutine data_moments(cgroup, xsdat, flux, aflux)
    implicit none

    type(coarse_group), intent(inout) :: cgroup
    type(xs_data), dimension(:), intent(in) :: xsdat            ! (region)
    real(kind=8), dimension(:,:), intent(in) :: flux            ! (region, group)
    real(kind=8), dimension(:,:,:), intent(in) :: aflux         ! (region, group, angle)
    
    integer :: i,j,K,N,region,regions,angle
    real(kind=8), dimension(:,:), allocatable :: delta_g        ! (region, group)
    
    N=cgroup%N-1
    regions=size(aflux,1)-1
    allocate(delta_g(regions,0:N))
        
    do region=1,regions
    
        ! chi moments
        cgroup%chi(region,:) = 0.0
        do K=0,N
            do i=0,N
                cgroup%chi(region,i)=cgroup%chi(region,i)+ &
                    & cgroup%DLP%values(i,K)*xsdat(region)%chi(cgroup%mingroup+K)
            enddo
        enddo
        
        ! sigt0
        cgroup%sigt0(region) = 0.0
        do K=0,N
            cgroup%sigt0(region)=cgroup%sigt0(region)+ &
                & xsdat(region)%sigt(cgroup%mingroup+K)*flux(region,cgroup%mingroup+K)
        enddo
        cgroup%sigt0(region)=cgroup%sigt0(region)/cgroup%scal_flux(region)
        
        ! delta_g
        delta_g(region,:) = 0.0
        do K=0,N
            delta_g(region,K)=xsdat(region)%sigt(cgroup%mingroup+K)-cgroup%sigt0(region)
        enddo
        
        ! delta moments
        cgroup%delta(region,:,:) = 0.0
        do angle=1,size(aflux,3)
            do K=0,N
                do i=0,N
                    cgroup%delta(region,i,angle)=cgroup%delta(region,i,angle) + &
                        & cgroup%DLP%values(i,K)*delta_g(region,K)* &
                        & aflux(region,cgroup%mingroup+K,angle)
                enddo
            enddo
            cgroup%delta(region,:,angle)=cgroup%delta(region,:,angle)/cgroup%ang_flux(region,0,angle)
        enddo

        
        ! vsigf
        cgroup%vsigf(region) = 0.0
        do K=0,N
            cgroup%vsigf(region)=cgroup%vsigf(region)+ & 
                & xsdat(region)%vsigf(cgroup%mingroup+K)*flux(region,cgroup%mingroup+K)
        enddo
        cgroup%vsigf(region)=cgroup%vsigf(region)/cgroup%scal_flux(region)
   
   
    enddo

end subroutine data_moments

!==============================================================================

subroutine scat_moments_down(origin, destination, xsdat, flux)
    implicit none
    
    type(coarse_group), intent(in) :: origin
    type(coarse_group), intent(inout) :: destination
    type(xs_data), dimension(:), intent(in) :: xsdat                ! (region)
    real(kind=8), dimension(:,:), intent(in) :: flux                ! (region, group)                        
    integer :: i, K, L, N, M, region
    integer :: dest_fg, scat_diag, scat_corner, maxL
    real(kind=8) :: temp
    
    M=origin%N-1
    N=destination%N-1
    
    do region=1,size(flux,1)
    
        destination%sig_dn(region,:,origin%index)=0.0
        do K=0,N
            dest_fg = destination%mingroup+K
            scat_diag = dest_fg*(dest_fg-1)/2+1
            scat_corner = scat_diag + destination%mingroup - origin%mingroup + K
            temp=0.0
            
            ! loop over less elements in block along diagonal
            if (origin%index==destination%index) then
                maxL=K
            else
                maxL=M
            endif
            
            ! condense reaction rate of origin groups
            do L=0,maxL
                temp=temp+xsdat(region)%sig_dn(scat_corner-L)* &
                    & flux(region,origin%mingroup+L)
            enddo
            
            ! DLP expansion of reaction rates in destination groups
            do i=0,N
                destination%sig_dn(region,i,origin%index) = &
                    & destination%sig_dn(region,i,origin%index) + &
                    & destination%DLP%values(i,K)*temp
            enddo
           
        enddo
        destination%sig_dn(region,:,origin%index)= &
            & destination%sig_dn(region,:,origin%index)/origin%scal_flux(region)

        !destination%sig_dn(region,1:N,origin%index) = 0.0


    enddo
    
end subroutine scat_moments_down

!==============================================================================

subroutine scat_moments_up(origin, destination, xsdat, flux)
    implicit none
    
    type(coarse_group), intent(in) :: origin
    type(coarse_group), intent(inout) :: destination
    type(xs_data), dimension(:), intent(in) :: xsdat                ! (region)
    real(kind=8), dimension(:,:), intent(in) :: flux                ! (region, group) 
    integer :: i,K,L,N,M,region
    integer :: orig_up, dest_up
    real(kind=8) :: temp
    
    M=origin%N-1
    N=destination%N-1
    
    orig_up = origin%mingroup-(groups-maxup)
    dest_up = destination%mingroup-(groups-maxup)
    
    do region=1,size(flux,1)
    
        destination%sig_up(region,:,origin%up_index)=0.0
        do K=0,N
            temp = 0.0
            
            ! condense reaction rate of origin groups
            do L=0,M
                temp=temp+xsdat(region)%sig_up(dest_up+K,orig_up+L)* &
                    & flux(region,origin%mingroup+L)
            enddo
            
            ! DLP expansion of reaction rates in destination groups
            do i=0,N
                destination%sig_up(region,i,origin%up_index) = &
                    & destination%sig_up(region,i,origin%up_index) + &
                    & destination%DLP%values(i,K)*temp
            enddo
        enddo
        destination%sig_up(region,:,origin%up_index)= &
            & destination%sig_up(region,:,origin%up_index)/origin%scal_flux(region)
            
        !destination%sig_up(region,1:N,origin%up_index) = 0.0
            
    enddo

end subroutine scat_moments_up

!==============================================================================

end module coarse_group_type
