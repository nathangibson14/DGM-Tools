module DLP_type

type discreteLegendre
    real(kind=8), pointer :: values(:,:) => null()
    real(kind=8), pointer :: a(:) => null()
end type discreteLegendre

contains

!==============================================================================

subroutine generate_DLPs(this,NN)
    implicit none
    type(discreteLegendre), intent(inout) :: this
    integer, intent(in) :: NN
    integer :: i, j, K, N
    real(kind=8) sumsquare
    
    N=NN-1

    ! generate DLPs
    allocate(this%values(0:N,0:N))

    do K=0,N
        this%values(0,K)=1.0
        if (N>0) then
            this%values(1,K)= real(N-2*K,kind=8)/real(N,kind=8)
        endif
    enddo

    do i=1,N-1
        do K=0,N
            this%values(i+1,K) = (real((2*i+1)*(N-2*K),kind=8)*this%values(i,K) - &
                & real(i*(N+i+1),kind=8)*this%values(i-1,K))/real((i+1)*(N-i),kind=8)
        enddo
    enddo
    
    ! generate orthogonality constants
    allocate(this%a(0:N))
    
    do i=0,N
        sumsquare=0.0
        do j=0,N
            sumsquare=sumsquare+this%values(i,j)**2
        enddo
        this%a(i) = 1/sumsquare
    enddo
    
end subroutine generate_DLPs

!==============================================================================

end module DLP_type