module dgm_solvers
use nuclear_data
use coarse_group_type
use solvers

contains

!==============================================================================

subroutine dgm_infinite_medium(xsdat,matl)
    implicit none

    integer :: i, g, gg, j, ii
    
    real(kind=8), dimension(1,groups) :: flux
    real(kind=8), dimension(1,groups,1) :: aflux
    
    
    type(xs_data), dimension(:), intent(in) :: xsdat
    integer, intent(in) :: matl
    
    type(xs_data) :: xsDGM

    integer, parameter :: max_DGM_steps = 100
    integer :: DGM_steps
    real(kind=8), dimension(0:max_DGM_steps) :: k_stored
    real(kind=8), dimension(coarse_groups) :: f_0
    real(kind=8) :: rms,avg,mre

    call allocate_xsdata(xsDGM, coarse_groups, upsc_groups)
    
    do g=1,coarse_groups
        allocate(CG(g)%scal_flux(1),CG(g)%ang_flux(1,0:CG(g)%N-1,1))
        allocate(CG(g)%vsigf(1), CG(g)%sigt0(1), CG(g)%delta(1,0:CG(g)%N-1,1))
        allocate(CG(g)%chi(1,0:CG(g)%N-1),CG(g)%sig_dn(1,0:CG(g)%N-1,CG(g)%index))
        if (g>coarse_groups-upsc_groups) then
            allocate(CG(g)%sig_up(1,0:CG(g)%N-1,upsc_groups))
        endif
    enddo
    
    
    ! k_inf = 1.0
    ! flux = 1.0
    ! flux = flux/sum(flux)
    
    k_inf = 1.0
    do i=1,groups
        flux(1,i)=1/(energy(i)*xsdat(matl)%sigt(i))
    enddo
    flux=flux/sum(flux)

    k_stored(0) = k_inf

    call start_timer()
    DGM_loop: do DGM_steps=1,max_DGM_steps

        write(*,*) 'Starting DGM step ', DGM_steps, read_timer()

        aflux(1,:,1) = flux(1,:)
        
        do i=1,coarse_groups
            call flux_moments(CG(i),flux,aflux)
            call data_moments(CG(i),(/xsdat(matl)/),flux,aflux)
        enddo

        do i=1,coarse_groups
            do j=1,i
                call scat_moments_down(CG(j),CG(i),xsdat,flux)
            enddo
            if (i>coarse_groups-upsc_groups) then
                do j=coarse_groups-upsc_groups+1,coarse_groups
                    call scat_moments_up(CG(j),CG(i),xsdat,flux)
                enddo
            endif
        enddo
        
        !------------------------!
        ! 0th order solution     !
        !------------------------!
        
        ii=1
        do g=1,coarse_groups

            
            xsDGM%sigt(g) = CG(g)%sigt0(1) + CG(g)%delta(1,0,1)
            xsDGM%vsigf(g) = CG(g)%vsigf(1)
            
            xsDGM%chi(g) = CG(g)%chi(1,0)
            
            ! downscatter matrix
            do j=g,1,-1
                xsDGM%sig_dn(ii) = CG(g)%sig_dn(1,0,j)
                ii=ii+1
            enddo

            ! upscatter matrix
            if (g>coarse_groups-upsc_groups) then
                xsDGM%sig_up(g-(coarse_groups-upsc_groups),:) = CG(g)%sig_up(1,0,:)
            endif
        enddo
        
        call inf_medium_power(xsDGM, f_0)

        !------------------------!
        ! higher order solutions !
        !------------------------!

        ! scattering source with downscatter only
        do g=1,coarse_groups
            do i=1,CG(g)%N-1
                if (g>coarse_groups-upsc_groups) then
                    gg = coarse_groups-upsc_groups
                else
                    gg = g
                endif
                
                ! downscatter source
                CG(g)%ang_flux(1,i,1) = dot_product(CG(g)%sig_dn(1,i,1:gg),f_0(1:gg))
                
                ! upscatter source
                if (g>coarse_groups-upsc_groups) then
                    CG(g)%ang_flux(1,i,1) = CG(g)%ang_flux(1,i,1) + dot_product(CG(g)%sig_up(1,i,:), & 
                        & f_0(coarse_groups-upsc_groups+1:coarse_groups))
                endif   
                
                ! add chi*vsigf/k
                CG(g)%ang_flux(1,i,1)=CG(g)%ang_flux(1,i,1)+ &
                    & CG(g)%chi(1,i)*dot_product(xsDGM%vsigf,f_0)/k_inf
                    
                ! subtract delta_i
                CG(g)%ang_flux(1,i,1)=CG(g)%ang_flux(1,i,1)-CG(g)%delta(1,i,1)*f_0(g)    
                
                ! divide by sigt for final answer
                CG(g)%ang_flux(1,i,1)=CG(g)%ang_flux(1,i,1)/CG(g)%sigt0(1)
                
                ! if (g==1 .and. i==1) then
                    ! write(*,*) 'ang_flux', CG(g)%ang_flux(i,1)
                    ! write(*,*) 'delta_i', CG(g)%delta_i(1)
                    ! write(*,*) 'chi_i', CG(g)%chi_i
                    ! write(*,*) 'sig_dn', CG(g)%sig_dn(g)
                ! endif
                
            enddo
        enddo


        !------------------------!
        ! condense and update    !
        !------------------------!    
        
        do i=1,coarse_groups
            CG(i)%ang_flux(1,0,1) = f_0(i)
        enddo
        
        ! condense to fine flux
        do i=1,coarse_groups
            call condense_flux_moments(CG(i),aflux)
        enddo
        flux(1,:)=aflux(1,:,1)/sum(aflux(1,:,1))
              
        ! flux updates
        call inf_medium_source(xsdat(matl), 3, .false., flux(1,:))
        write(*,*) '    k_update = ', k_inf
        
        ! calculate rms error
        rms=0.0
        do g=1,groups
            rms=rms+((flux(1,g)-phi(g))/phi(g))**2
        enddo
        rms=sqrt(rms/groups)
        write(*,*) '    rms = ', rms
        
        ! calculate mean relative error
        avg=0.0
        do g=1,groups
            avg=avg+phi(g)
        enddo
        avg=avg/groups
        
        mre=0.0
        do g=1,groups
            mre=mre+abs(flux(1,g)-phi(g))/avg
        enddo
        mre=mre/groups
        write(*,*) '    mre = ', mre
        
        
        k_stored(DGM_steps) = k_inf
        if (abs(k_stored(DGM_steps)-k_stored(DGM_steps-1))<1e-6) then
            write(*,*) 'k_inf converged at DGM step', DGM_steps
            exit DGM_loop
        endif
        
    enddo DGM_loop


    ! open(unit=21, file='phi', action='write', status='replace')
    ! do i=1,groups
        ! write(21,*) flux(1,i)
    ! enddo
    ! close(unit=21)


end subroutine dgm_infinite_medium

!==============================================================================

subroutine dgm_sn(xsdat,meshes,width,matl,angles,BC,diff_scheme)
    implicit none

    !------------------------!
    ! variable declaration   !
    !------------------------!
    
    type(xs_data), dimension(:), intent(in) :: xsdat
    integer, dimension(:), intent(in) :: meshes, matl
    real(kind=8), dimension(:), intent(in) :: width
    integer, intent(in) :: angles, diff_scheme
    real(kind=8), dimension(2), intent(in) :: BC
    
    integer :: regions
    real(kind=8) :: k_eff
    type(xs_data), dimension(:), allocatable :: xsDGM, xs_reg
    real(kind=8), dimension(:,:,:), allocatable :: delta            ! (region,group,angle)
    real(kind=8), dimension(:,:), allocatable :: flux, dgm_flux     ! (mesh,group)
    real(kind=8), dimension(:,:,:), allocatable :: aflux, dgm_aflux ! (mesh,group,angle)
    real(kind=8), dimension(:,:), allocatable :: reg_flux           ! (region,group)
    real(kind=8), dimension(:,:,:), allocatable :: reg_aflux        ! (region,group,angle)    
    integer, dimension(:), allocatable :: regID, matID
    
    real(kind=8), dimension(:), allocatable :: Q_f0, Q_f, Q_ds, Q_us, Q
    real(kind=8), dimension(:), allocatable :: mu, wgt
    real(kind=8), dimension(:), allocatable :: alpha, dx
    real(kind=8), dimension(:,:), allocatable :: A, B
    real(kind=8) :: denom
    integer :: smu
    
    integer :: i,ii,j,g,region,angle
    
    !------------------------!
    ! allocate data          !
    !------------------------!
    
    regions = size(meshes)
    allocate(xsDGM(regions),xs_reg(regions))
    do i=1,regions
        call allocate_xsdata(xsDGM(i), coarse_groups, upsc_groups)
    enddo
    allocate(delta(regions,coarse_groups,angles))
    allocate(flux(sum(meshes),groups), aflux(sum(meshes)+1,groups,angles))
    allocate(reg_flux(regions,groups), reg_aflux(regions,groups,angles))
    allocate(dgm_flux(sum(meshes),coarse_groups), dgm_aflux(sum(meshes)+1,coarse_groups,angles))
    allocate(regID(sum(meshes)),matID(sum(meshes)))
    allocate(Q_f0(sum(meshes)),Q_f(sum(meshes)),Q_ds(sum(meshes)),Q_us(sum(meshes)),Q(sum(meshes)))
    
    do g=1,coarse_groups
        allocate(CG(g)%scal_flux(regions),CG(g)%ang_flux(regions,0:CG(g)%N-1,angles))
        allocate(CG(g)%vsigf(regions), CG(g)%sigt0(regions), CG(g)%delta(regions,0:CG(g)%N-1,angles))
        allocate(CG(g)%chi(regions,0:CG(g)%N-1),CG(g)%sig_dn(regions,0:CG(g)%N-1,CG(g)%index))
        if (g>coarse_groups-upsc_groups) then
            allocate(CG(g)%sig_up(regions,0:CG(g)%N-1,upsc_groups))
        endif
    enddo
    
    allocate(alpha(angles), A(sum(meshes),angles),B(sum(meshes),angles))
    allocate(dx(size(meshes)))

    !------------------------!
    ! weights and angles     !
    !------------------------!
    
    if (angles==2) then
        allocate(mu(2), wgt(2))
        mu = (/-0.5773502691, 0.5773502691/)
        wgt = (/1.0, 1.0/)
    elseif (angles==4) then
        allocate(mu(4), wgt(4))
        mu = (/-0.8611363115, -0.3399810435, 0.3399810435, 0.8611363115/)
        wgt = (/0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451/)
    elseif (angles==8) then
        allocate(mu(8), wgt(8))
        mu = (/-0.9602898564,-0.7966664774,-0.5255324099,-0.1834346424, &
            & 0.1834346424,0.5255324099,0.7966664774,0.9602898564 /)
        wgt = (/ 0.1012285363,0.2223810344,0.3137066459,0.3626837834, &
            & 0.3626837834,0.3137066459,0.2223810344,0.1012285363 /)
    else
        write(*,*) 'Invalid N in S_N solver.'
        stop
    endif
    
    !------------------------!
    ! geometry setup         !
    !------------------------!
    
    ii=1
    do i=1,sum(meshes)
        if (i>sum(meshes(1:ii))) then
            ii=ii+1
        endif
        regID(i)=ii
        matID(i)=matl(ii)
    enddo
    
    do i=1,regions
        ! xs_reg(i) = xsdat(matID(i))
        xs_reg(i) = xsdat(matl(i))
    enddo
    
    do i=1,size(dx)
        dx(i) = (width(i+1)-width(i))/real(meshes(i),kind=8)
    enddo
    
    !------------------------!
    ! initialize DGM data    !
    !------------------------!
    
    flux=1.0
    aflux=1.0
    
    write(*,*) 'starting sn_solver... ', read_timer()
    call sn_solver(xs,meshes,width,matl,angles,BC,diff_scheme,aflux=aflux,flux=flux)
    write(*,*) 'back to dgm...', read_timer()
    
    reg_flux=0.0
    reg_aflux=0.0
    do i=1,sum(meshes)
        reg_flux(regID(i),:)=reg_flux(regID(i),:)+flux(i,:)
        if (diff_scheme==1) then
            reg_aflux(regID(i),:,:)=reg_aflux(regID(i),:,:)+.5*(aflux(i,:,:)+aflux(i+1,:,:))
        else
            do angle=1,angles
                if (angle>angles/2) then
                    reg_aflux(regID(i),:,angle)=reg_aflux(regID(i),:,angle)+aflux(i+1,:,angle)
                else
                    reg_aflux(regID(i),:,angle)=reg_aflux(regID(i),:,angle)+aflux(i,:,angle)
                endif
            enddo
        endif
    enddo

    do g=1,coarse_groups
        call flux_moments(CG(g),reg_flux,reg_aflux)
        call data_moments(CG(g),xs_reg,reg_flux,reg_aflux)
    enddo
    
    do i=1,coarse_groups
        do j=1,i
            call scat_moments_down(CG(j),CG(i),xs_reg,reg_flux)
        enddo
        if (i>coarse_groups-upsc_groups) then
            do j=coarse_groups-upsc_groups+1,coarse_groups
                call scat_moments_up(CG(j),CG(i),xs_reg,reg_flux)
            enddo
        endif
    enddo
    
    !------------------------!
    ! 0th order solution     !
    !------------------------!
   
    do region=1,regions
        ii=1
        do g=1,coarse_groups
            xsDGM(region)%sigt(g) = CG(g)%sigt0(region)
            delta(region,g,:) = CG(g)%delta(region,0,:)
            xsDGM(region)%vsigf(g) = CG(g)%vsigf(region)
            xsDGM(region)%chi(g) = CG(g)%chi(region,0)
            
            ! downscatter matrix
            do j=g,1,-1
                xsDGM(region)%sig_dn(ii) = CG(g)%sig_dn(region,0,j)
                ii=ii+1
            enddo

            ! upscatter matrix
            if (g>coarse_groups-upsc_groups) then
                xsDGM(region)%sig_up(g-(coarse_groups-upsc_groups),:) = CG(g)%sig_up(region,0,:)
            endif
        
        
        enddo
    enddo
    
    call sn_solver(xsDGM,meshes,width,(/(i,i=1,regions)/),angles,BC,diff_scheme,k_eff,dgm_aflux,dgm_flux,delta)
    
    
    
    !------------------------!
    ! higher order solutions !
    !------------------------!

    Q_f0 = 0.0
    do i=1,sum(meshes)
        do g=1,coarse_groups
            Q_f0(i)=Q_f0(i)+CG(g)%vsigf(regID(i))*dgm_flux(i,g)
        enddo
    enddo
 
    dgm_aflux=0.0
    group_sweep: do g=1,coarse_groups
        moment_sweep: do j=1,CG(g)%N-1

            ! build source            
            do i=1,sum(meshes)
                Q_f(i) = CG(g)%chi(regID(i),j)*Q_f0(i)/k_eff
                if (g<=coarse_groups-upsc_groups) then
                    Q_us=0.0
                    Q_ds=0.0
                    do ii=1,g
                        Q_ds(i)=Q_ds(i)+CG(g)%sig_dn(regID(i),j,ii)*dgm_flux(i,ii)
                    enddo
                else
                    Q_us=0.0
                    Q_ds=0.0
                    do ii=1,coarse_groups-upsc_groups
                        Q_ds(i)=Q_ds(i)+CG(g)%sig_dn(regID(i),j,ii)*dgm_flux(i,ii)
                    enddo
                    do ii=1,upsc_groups
                        Q_us(i)=Q_us(i)+CG(g)%sig_up(regID(i),j,ii)* & 
                            & dgm_flux(i,coarse_groups-upsc_groups+ii)
                    enddo
                endif
            enddo
            Q = .5*(Q_f + Q_ds + Q_us)
            
            ! compute coefficients
            do angle=1,angles
                alpha(angle) = mu(angle)/abs(mu(angle))*(-diff_scheme+1)
                if (angle<=angles/2) then
                    smu=-1
                else
                    smu=1
                endif                
                
                do i=1,sum(meshes)
                    denom = 2*mu(angle)+smu*(1+smu*alpha(angle))*CG(g)%sigt0(regID(i))*dx(regID(i))
                    A(i,angle) = (2*mu(angle)-smu*(1-smu*alpha(angle))*CG(g)%sigt0(regID(i))*dx(regID(i)))/denom
                    B(i,angle) = smu*2*dx(regID(i))/denom
                enddo               
            enddo
            
            ! loop through negative angles
            do angle=1,angles/2
                ! boundary flux
                dgm_aflux(sum(meshes)+1,g,angle)=BC(2)*dgm_aflux(sum(meshes)+1,g,angles-angle+1)
                
                ! interior fluxes
                do i=sum(meshes),1,-1                    
                    ! solve for next mesh point
                    dgm_aflux(i,g,angle)=A(i,angle)*dgm_aflux(i+1,g,angle) & 
                        & + B(i,angle)*(Q(i)-CG(g)%delta(regID(i),j,angle)*dgm_aflux(i,g,angle))
                enddo
                
            enddo
            
            ! loop through positive angles
            do angle=angles/2+1,angles
                ! boundary flux
                dgm_aflux(1,g,angle)=BC(1)*dgm_aflux(1,g,angles-angle+1)
                
                ! interior fluxes
                do i=1,sum(meshes)                    
                    ! solve for next mesh point
                    dgm_aflux(i+1,g,angle)=A(i,angle)*dgm_aflux(i,g,angle) & 
                        & + B(i,angle)*(Q(i)-CG(g)%delta(regID(i),j,angle)*dgm_aflux(i,g,angle))                       
                enddo                
            enddo
            
        enddo moment_sweep
    enddo group_sweep
    
end subroutine dgm_sn

!==============================================================================

end module dgm_solvers