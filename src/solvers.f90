module solvers
use nuclear_data
use timer

real(kind=8) :: k_inf

contains

!==============================================================================

subroutine inf_medium_power(xsdat, flux)
    implicit none 
    
    type(xs_data), intent(in) :: xsdat
    real(kind=8), dimension(xsdat%groups), intent(out) :: flux
    
    real(kind=8), dimension(xsdat%groups) :: flux0
    real(kind=8), dimension(xsdat%maxup) :: flux_up
    real(kind=8) :: k_old
    
    integer :: i,j,ii,jj, iter
    real(kind=8), dimension(xsdat%maxscat) :: A
    real(kind=8), dimension(xsdat%groups) :: b
    real(kind=8), dimension(xsdat%maxup,xsdat%maxup) :: AA
    
    k_old=1.0
    flux0=1.0
    flux0=flux0/sum(flux0)

    ! build LHS
    A = -xsdat%sig_dn
    j=1
    do i=1,xsdat%groups-xsdat%maxup
        A(j) = A(j)+xsdat%sigt(i)
        j=j+i
    enddo
    AA = -xsdat%sig_up
    do i=1,xsdat%maxup
        AA(i,i)=AA(i,i)+xsdat%sigt(xsdat%groups-xsdat%maxup+i)
    enddo

    ! build RHS
    b = dot_product(xsdat%vsigf,flux0)*xsdat%chi/k_old
    
    !! begin direct Gauss-Seidel solve
    ! downscatter block
    j=1
    do i=1,xsdat%groups-xsdat%maxup
        flux(i) = b(i)
        do ii=1,i-1
            flux(i)=flux(i)-A(j+ii)*flux(i-ii)
        enddo
        flux(i)=flux(i)/A(j)
        j=j+i
    enddo

    ! build RHS for upscatter block
    do i=1,xsdat%maxup
        ii=xsdat%groups-xsdat%maxup+i
        do jj=1,xsdat%groups-xsdat%maxup
            b(ii)=b(ii)-A(j+i+jj-1)*flux(ii-i-jj+1)
        enddo
        j=j+ii
    enddo

    ! upscatter block
    iter=0
    flux_up = 0.0
    flux(xsdat%groups-xsdat%maxup+1:xsdat%groups)=flux0(xsdat%groups-xsdat%maxup+1:xsdat%groups)
    do while (maxval(abs(flux(xsdat%groups-xsdat%maxup+1:xsdat%groups)-flux_up))>1e-10)
        flux_up=flux(xsdat%groups-xsdat%maxup+1:xsdat%groups)
        do i=1,xsdat%maxup
            j=xsdat%groups-xsdat%maxup+i
            flux(j) = b(j)
            do ii=1,xsdat%maxup
                jj=xsdat%groups-xsdat%maxup+ii
                if (ii/=i) then
                    flux(j)=flux(j)-AA(i,ii)*flux(jj)
                endif
            enddo
            flux(j)=flux(j)/AA(i,i)
        enddo
        iter=iter+1
        if (iter==1000) then
            write(*,*) 'Exceeded Maximum GS Iterations!'
            stop
        endif
    enddo
    !write(*,*) '    iter =', iter

    k_inf = k_old*dot_product(xsdat%vsigf,flux)/dot_product(xsdat%vsigf,flux0)
    flux = flux/sum(flux)

    write(*,*) '    k_inf = ', k_inf

end subroutine inf_medium_power

!==============================================================================

subroutine inf_medium_source(xsdat, maxiter, check_conv, flux)
    implicit none
    
    integer :: iter,i,ii,j,jj
    
    type(xs_data), intent(in) :: xsdat
    real(kind=8), dimension(xsdat%groups), intent(inout) :: flux
    
    integer, intent(in) :: maxiter
    logical, intent(in) :: check_conv
    real(kind=8), dimension(:), allocatable :: flux0
    real(kind=8), dimension(:), allocatable :: flux_up
    real(kind=8) :: k_old
    
    allocate(flux0(xsdat%groups), flux_up(xsdat%maxup))
        
    iter_loop: do iter=1,maxiter
        k_old=k_inf
        flux0=flux
        flux=0.0
        
        ! S*flux
        do i=1,xsdat%groups-xsdat%maxup
            j=i*(i-1)/2+1
            do ii=0,i-1
                flux(i)=flux(i)+xsdat%sig_dn(j+ii)*flux0(i-ii)                
            enddo
        enddo
        flux(xsdat%groups-xsdat%maxup+1:xsdat%groups)= & 
            & matmul(xsdat%sig_up,flux0(xsdat%groups-xsdat%maxup+1:xsdat%groups))
        
        ! build RHS for upscatter block
        j=(xsdat%groups-xsdat%maxup+1)*(xsdat%groups-xsdat%maxup)/2+1
        do i=1,xsdat%maxup
            ii=xsdat%groups-xsdat%maxup+i
            do jj=1,xsdat%groups-xsdat%maxup
                flux(ii)=flux(ii)+xsdat%sig_dn(j+i+jj-1)*flux0(ii-i-jj+1)
            enddo
            j=j+ii
        enddo
        
        ! add chi*vsigf*flux/k
        flux = flux + dot_product(xsdat%vsigf,flux0)*xsdat%chi/k_old
        
        ! diagonal matrix solve
        do i=1,xsdat%groups
            flux(i) = flux(i)/xsdat%sigt(i)
        enddo
        
        ! update eigenvalue
        k_inf = k_old*dot_product(xsdat%vsigf,flux)/dot_product(xsdat%vsigf,flux0)
        
        flux=flux/sum(flux)
        
        if (abs(k_inf-k_old) < 1e-6 .and. check_conv) then
            write(*,*) '    iter = ', iter
            write(*,*) '    k_inf = ', k_inf
            exit iter_loop
        endif
       
    enddo iter_loop
    
end subroutine inf_medium_source


!==============================================================================

subroutine sn_solver(xsdat, meshes, width, matl, N, BC, diff_scheme, k_out, aflux, flux, delta)
    implicit none
    
    integer, dimension(:), intent(in) :: meshes, matl
    real(kind=8), dimension(:), intent(in) :: width
    type(xs_data), dimension(:), intent(in) :: xsdat
    real(kind=8), dimension(:,:,:), intent(in), optional :: delta
    
    integer, intent(in) :: N, diff_scheme ! diff_scheme: 0-SD, 1-DD
    real(kind=8), dimension(2), intent(in) :: BC
    real(kind=8), dimension(:), allocatable :: mu, wgt
    real(kind=8), dimension(:), allocatable :: dx
    
    real(kind=8), dimension(sum(meshes)+1,xsdat(1)%groups,N) :: ang_flux, ang_flux_old
    real(kind=8), dimension(sum(meshes)+1,xsdat(1)%groups,N), intent(out), optional :: aflux
    real(kind=8), dimension(sum(meshes),xsdat(1)%groups):: scal_flux, scal_flux_old
    real(kind=8), dimension(sum(meshes),xsdat(1)%groups), intent(out), optional :: flux
    real(kind=8) :: k_eff, k_old
    real(kind=8), intent(out), optional :: k_out
    
    real(kind=8), dimension(sum(meshes)) :: S_f0
    integer, dimension(sum(meshes)) :: regID, matID
    real(kind=8), dimension(sum(meshes),xsdat(1)%groups) :: S_f, S_ds
    real(kind=8), dimension(sum(meshes),xsdat(1)%maxup) :: S_us
    real(kind=8) :: RR_f, Q, smu, denom
    real(kind=8), dimension(sum(meshes),xsdat(1)%groups,N) :: A, B
    real(kind=8), dimension(N) :: alpha
    integer :: power_iter, source_iter, upsc_iter, g, g_up
    integer :: i, j, ii, gg, maxgg
    
    ! weights and angles
    if (N==2) then
        allocate(mu(2), wgt(2))
        mu = (/-0.5773502691, 0.5773502691/)
        wgt = (/1.0, 1.0/)
    elseif (N==4) then
        allocate(mu(4), wgt(4))
        mu = (/-0.8611363115, -0.3399810435, 0.3399810435, 0.8611363115/)
        wgt = (/0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451/)
    elseif (N==8) then
        allocate(mu(8), wgt(8))
        mu = (/-0.9602898564,-0.7966664774,-0.5255324099,-0.1834346424, &
            & 0.1834346424,0.5255324099,0.7966664774,0.9602898564 /)
        wgt = (/ 0.1012285363,0.2223810344,0.3137066459,0.3626837834, &
            & 0.3626837834,0.3137066459,0.2223810344,0.1012285363 /)
    else
        write(*,*) 'Invalid N in S_N solver.'
        stop
    endif
    
    ! geometry
    allocate(dx(size(meshes)))
    do i=1,size(dx)
        dx(i) = (width(i+1)-width(i))/real(meshes(i),kind=8)
    enddo
    
    ii=1
    do i=1,sum(meshes)
        if (i>sum(meshes(1:ii))) then
            ii=ii+1
        endif
        regID(i)=ii
        matID(i)=matl(ii)
    enddo
        
    ! guess
    k_old=0.0
    k_eff=1.0
    ang_flux=1.0
    scal_flux=1.0
    
    ! normalize by fission reaction rate
    RR_f = 0.0
    do g=1,xsdat(1)%groups
        do i=1,sum(meshes)
            RR_f = RR_f+xsdat(matID(i))%vsigf(g)*scal_flux(i,g)
        enddo
    enddo
    scal_flux=scal_flux/RR_f
    
    ! precompute coefficients
    do g=1,xsdat(1)%groups
        do j=1,N
            alpha(j) = mu(j)/abs(mu(j))*(-diff_scheme+1)
            if (j<=N/2) then
                smu=-1
            else
                smu=1
            endif
            
            if (present(delta)) then
                do i=1,sum(meshes)
                    denom = 2*mu(j)+smu*(1+smu*alpha(j))*(xsdat(matID(i))%sigt(g)+delta(regID(i),g,j))*dx(regID(i))
                    A(i,g,j) = (2*mu(j)-smu*(1-smu*alpha(j))*(xsdat(matID(i))%sigt(g)+delta(regID(i),g,j))*dx(regID(i)))/denom
                    B(i,g,j) = smu*2*dx(regID(i))/denom
                enddo
            else
                do i=1,sum(meshes)
                    denom = 2*mu(j)+smu*(1+smu*alpha(j))*xsdat(matID(i))%sigt(g)*dx(regID(i))
                    A(i,g,j) = (2*mu(j)-smu*(1-smu*alpha(j))*xsdat(matID(i))%sigt(g)*dx(regID(i)))/denom
                    B(i,g,j) = smu*2*dx(regID(i))/denom
                enddo
            endif
        enddo
    enddo
    
    power_iter_loop: do power_iter=1,100
        
        write(*,*) '      k_eff = ', k_eff
        write(*,*) '   starting power iteration         ', power_iter, read_timer()
        
        ! update stored values
        k_old = k_eff
        scal_flux_old = scal_flux
        
        ! build fission source
        S_f0=0.0
        do i=1,sum(meshes)
            do g=1,xsdat(1)%groups
                S_f0(i) = S_f0(i) + xsdat(matID(i))%vsigf(g)*scal_flux(i,g)/k_eff
            enddo
            do g=1,xsdat(1)%groups
                S_f(i,g) = xsdat(matID(i))%chi(g)*S_f0(i)/2
            enddo
        enddo
        
        ! write(*,*) '    starting downscatter block        ', read_timer()       
        group_sweep_down: do g=1,xsdat(1)%groups-xsdat(1)%maxup
            
            ! build downscatter scource
            S_ds(:,g)=0.0
            do i=1,sum(meshes)
                do gg=1,g-1
                    ii=g*(g+1)/2
                    S_ds(i,g)=S_ds(i,g)+xsdat(matID(i))%sig_dn(ii-gg+1)*scal_flux(i,gg)/2
                enddo
            enddo
        
            source_iter_loop1: do source_iter=1,100
                ! update stored value
                ang_flux_old(:,g,:) = ang_flux(:,g,:)
                
                ! loop through negative angles
                do j=1,N/2
                    ! boundary flux
                    ang_flux(sum(meshes)+1,g,j)=BC(2)*ang_flux(sum(meshes)+1,g,N-j+1)
                    
                    ! interior fluxes
                    do i=sum(meshes),1,-1
                        ! build source
                        Q = S_f(i,g)+S_ds(i,g)+xsdat(matID(i))%sig_dn(g*(g-1)/2+1)*scal_flux(i,g)/2
                        
                        ! solve for next mesh point
                        ang_flux(i,g,j)=A(i,g,j)*ang_flux(i+1,g,j)+B(i,g,j)*Q
                        
                    enddo
                    
                enddo
                
                ! loop through positive angles
                do j=N/2+1,N
                    ! boundary flux
                    ang_flux(1,g,j)=BC(1)*ang_flux(1,g,N-j+1)
                    
                    ! interior fluxes
                    do i=1,sum(meshes)
                        ! build source
                        Q = S_f(i,g)+S_ds(i,g)+xsdat(matID(i))%sig_dn(g*(g-1)/2+1)*scal_flux(i,g)/2
                        
                        ! solve for next mesh point
                        ang_flux(i+1,g,j)=A(i,g,j)*ang_flux(i,g,j)+B(i,g,j)*Q                        
                    enddo                
                enddo
                
                ! condense to scalar flux
                scal_flux(:,g)=0.0
                do j=1,N
                    scal_flux(:,g)=scal_flux(:,g)+ang_flux(:,g,j)*wgt(j)
                enddo
                
                if (maxval(abs(ang_flux(:,g,:)-ang_flux_old(:,g,:)))<1e-12) then
                    !if (g==1) write(*,*) source_iter 
                    exit source_iter_loop1
                endif
        
            enddo source_iter_loop1
      
        enddo group_sweep_down
        
        !write(*,*) '    starting upscatter block        ', read_timer()
        ! build downscatter scource
        do g=xsdat(1)%groups-xsdat(1)%maxup+1,xsdat(1)%groups
            S_ds(:,g)=0.0
            do i=1,sum(meshes)
                do gg=1,xsdat(1)%groups-xsdat(1)%maxup
                    ii=g*(g+1)/2
                    S_ds(i,g)=S_ds(i,g)+xsdat(matID(i))%sig_dn(ii-gg+1)*scal_flux(i,gg)/2
                enddo
            enddo
        enddo
        
        upsc_iter_loop: do upsc_iter=1,1

            group_sweep_up: do g_up=1,xsdat(1)%maxup
                do i=1,sum(meshes)
                    S_us(i,g_up) = dot_product(xsdat(matID(i))%sig_up(g_up,:), &
                        & scal_flux(i,xsdat(1)%groups-xsdat(1)%maxup+1:xsdat(1)%groups))/2 &
                        & - xsdat(matID(i))%sig_up(g_up,g_up)*scal_flux(i,xsdat(1)%groups-xsdat(1)%maxup+g_up)/2
                enddo
                g = xsdat(1)%groups-xsdat(1)%maxup+g_up
                
                source_iter_loop2: do source_iter=1,100
                
                    ! update stored value
                    ang_flux_old(:,g,:) = ang_flux(:,g,:)
                
                    ! loop through negative angles
                    do j=1,N/2
                        ! boundary flux
                        ang_flux(sum(meshes)+1,g,j)=BC(2)*ang_flux(sum(meshes)+1,g,N-j+1)
                        
                        ! interior fluxes
                        do i=sum(meshes),1,-1
                            ! build source
                            Q = S_f(i,g)+S_ds(i,g)+S_us(i,g_up)
                            Q=Q+xsdat(matID(i))%sig_up(g_up,g_up)*scal_flux(i,g)/2
                            
                            ! solve for next mesh point
                            ang_flux(i,g,j)=A(i,g,j)*ang_flux(i+1,g,j)+B(i,g,j)*Q                      
                        enddo
                        
                    enddo
                    
                    ! loop through positive angles
                    do j=N/2+1,N
                        ! boundary flux
                        ang_flux(1,g,j)=BC(1)*ang_flux(1,g,N-j+1)
                        
                        ! interior fluxes
                        do i=1,sum(meshes)
                            ! build source
                            Q = S_f(i,g)+S_ds(i,g)+S_us(i,g_up)
                            Q=Q+xsdat(matID(i))%sig_up(g_up,g_up)*scal_flux(i,g)/2
                            
                            ! solve for next mesh point
                            ang_flux(i+1,g,j)=A(i,g,j)*ang_flux(i,g,j)+B(i,g,j)*Q                       
                        enddo 
                    enddo
                    
                    ! condense to scalar flux
                    scal_flux(:,g)=0.0
                    if (diff_scheme==1) then
                        do j=1,N
                            scal_flux(:,g)=scal_flux(:,g)+.5*(ang_flux(1:sum(meshes),g,j)+ & 
                                & ang_flux(2:sum(meshes)+1,g,j))*wgt(j)
                        enddo
                    else
                        do j=1,N
                            if (j>N/2) then
                                scal_flux(:,g)=scal_flux(:,g)+ang_flux(2:sum(meshes)+1,g,j)*wgt(j)
                            else
                                scal_flux(:,g)=scal_flux(:,g)+ang_flux(1:sum(meshes),g,j)*wgt(j)
                            endif
                        enddo
                    endif
                    
                    if (maxval(abs(ang_flux(:,g,:)-ang_flux_old(:,g,:)))<1e-12) then
                        exit source_iter_loop2
                    endif
                    
                enddo source_iter_loop2
            
            enddo group_sweep_up
            
            
        enddo upsc_iter_loop
    
    
        ! update eigenvalue
        RR_f = 0.0
        do g=1,xsdat(1)%groups
            do i=1,sum(meshes)
                RR_f = RR_f+xsdat(matID(i))%vsigf(g)*scal_flux(i,g)
            enddo
        enddo
        k_eff = k_eff*RR_f
        
        ! normalize by fission reaction rate
        scal_flux=scal_flux/RR_f
        
        if (abs(k_eff-k_old)<1e-6) then
            write(*,*) '   k_eff converged at step', power_iter
            write(*,*) '   k_eff = ', k_eff
            exit power_iter_loop
        endif
    
    enddo power_iter_loop
    
    if (present(k_out)) then
        k_out = k_eff
    endif
    if (present(aflux)) then
        aflux = ang_flux
    endif
    if (present(flux)) then
        flux = scal_flux
    endif
    
    
end subroutine sn_solver

!==============================================================================

end module solvers