program rt_2d ! The program simulates a primitive ray-tracing at an azimuth of 45.
              ! Each ray starts either in the first column j = 1 or in the last row
              ! (i' = i - 1, j' = j + 1) until it stops either at first row i = 1 or
              ! the last column j = n_j.
    implicit none
    real, allocatable :: atmos(:,:)
    integer :: n_i, n_j ! 2D grid size, row 'i', columns 'j'.

    print '(a)', "Grid size of the atmosphere?"
    read *, n_i, n_j
    allocate(atmos(n_i,n_j)) ! THE GRID IS EQUIDISTANT AND NON-PERIODIC
    call init_atmos()
    call print_atmos("Initial")
    call solve()
    call print_atmos("Final")
    deallocate(atmos)


contains
        ! subroutines init_atmos, print_atmos, solve, and propagate
   subroutine init_atmos ! The first subroutine sets up a 2D array atmos,
                         !which we will call the atmosphere for simplicity
       implicit none
       integer :: i, j
       do j = 1, n_j
           do i = 1, n_i
               atmos(i, j) = 1.42857143 * i + 5.8823529 * j
           end do
       end do
   end subroutine init_atmos ! There is no input file, the code 
                             !itself generates the initial data


   subroutine print_atmos(name) ! Prints out the entire atmosphere row-by-row
       implicit none
       character(len=*), intent(in) :: name
       integer :: i
       print '(2a)', name, " values:"
       do i = 1, n_i
           print '(*(es12.6, 1x))', atmos(i,:)
       end do
   end subroutine print_atmos ! It is called two times in the main program, 
                              ! before and after the solution


   subroutine solve() ! Calls the last subroutine propagate for all grid points
                      ! in the first column and in the last row
       implicit none
       integer :: i, j

       ! Rays starting in the 1st column:
       j = 1
       do i = 1, n_i
           call propagate(i, j)
       end do

       ! Rays starting in the last row:
       i = n_i
       do j = 2, n_j
           call propagate(i, j)
       end do
   end subroutine solve


   subroutine propagate(first_i, first_j) ! Traces a ray at an azimuth of 45 and 
                                          ! updates the atmosphere at each step on it
       implicit none
       integer, intent(in) :: first_i, first_j
       integer             :: crnt_i, next_i, crnt_j, next_j

       crnt_i = first_i
       crnt_j = first_j
       do 
           next_i = crnt_i - 1
           next_j = crnt_j + 1
           if (next_i < 1 .or. next_j > n_j) exit
           atmos(next_i, next_j) = &
                   (atmos(crnt_i, crnt_j) + atmos(next_i, next_j)) / 2.0
           crnt_i = next_i
           crnt_j = next_j
       end do
   end subroutine propagate
        


end program rt_2d
