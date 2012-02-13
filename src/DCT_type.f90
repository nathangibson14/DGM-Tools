module DCT_type
use parameters

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
    
    do i=0,N
        do K=1,N
            this%values(i,K) = cos((real(i,kind=8)*pi/real(NN,kind=8))*(real(K,kind=8)+.5_8))
        enddo
    enddo

   
    ! generate orthogonality constants
    allocate(this%a(0:N))
    
    this%a(0) = 1._8/real(NN,kind=8)
    do i=1,N
        this%a(i) = 2._8/real(NN,kind=8)
    enddo
    
end subroutine generate_DLPs

!==============================================================================

end module DCT_type
