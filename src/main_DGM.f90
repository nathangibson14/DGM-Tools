program main_DGM
use dgm_timer
use DCT_type
use parameters
use nuclear_data
use solvers
use coarse_group_type
use dgm_solvers

implicit none

real(kind=8), dimension(2) :: BC = 1.0
integer :: i

!------------------------------------------------------------------------------

call start_timer()

call input_coarse_groups()

allocate(energy(0:groups))
allocate(phi(groups))

!------------------------------------------------------------------------------

!~ write(*,*) 'starting grab cross sections... ', read_timer()
!~ call grab_xs()

!~ write(*,*) 'starting input my xs... ', read_timer()
!~ call input_my_xs()

write(*,*) 'starting input from gendf... ', read_timer()
call input_from_gendf()

!~ write(*,*) 'starting input data binary... ', read_timer()
!~ call input_data_binary

!~ write(*,*) 'starting input data ascii... ', read_timer()
!~ call input_data_ascii

write(*,*) 'starting inf_medium_power... ', read_timer()
call inf_medium_power(xs(1), phi)

!~ write(*,*) 'starting inf_medium_source... ', read_timer()
!~ k_inf = 1.0
!~ phi = 1.0
!~ phi = phi/sum(phi)
!~ call inf_medium_source(xs(1), 1000, .true., phi)

write(*,*) 'output data ascii... ', read_timer()
call output_data_ascii(xs(1))

!~ write(*,*) 'calc few group... ', read_timer()
!~ call calc_few_group(xs(1))

!~ write(*,*) 'output data binary... ', read_timer()
!~ call output_data_binary

!~ write(*,*) 'starting sn_solver... ', read_timer()
!~ call sn_solver(xs,(/50,50/),(/0._8,1._8,2._8/),(/1,2/),8,BC,0)

!~ write(*,*) 'starting sn_solver... ', read_timer()
!~ call sn_solver(xs,(/10,20/),(/0._8,0.375_8,0.5625_8/),(/1,1/),8,BC,0)

!~ write(*,*) 'starting dgm_sn... ', read_timer()
!~ call dgm_sn(xs,(/(1,i=1,150)/),(/(real(i*1.25/150,kind=8),i=0,150)/), & 
!~     & (/(1,i=1,50),(1,i=1,50),(1,i=1,50)/),8,BC,0)

write(*,*) 'starting dgm_infinite_medium... ', read_timer()
call dgm_infinite_medium(xs,1)

!------------------------------------------------------------------------------

write(*,*) 'program complete... ', read_timer()


end program main_DGM
