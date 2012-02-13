module dgm_timer

real(kind=8), dimension(2) :: time

contains

subroutine start_timer()
    call cpu_time(time(1))
end subroutine start_timer

function read_timer()
    real(kind=8) :: read_timer

    call cpu_time(time(2))
    read_timer = time(2)-time(1)
    return 
end function read_timer

end module dgm_timer
