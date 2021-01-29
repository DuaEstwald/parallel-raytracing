! AUTHOR: ELENA ARJONA GALVEZ

program rt_2d ! The program simulates a primitive ray-tracing at an azimuth of 45.
              ! Each ray starts either in the first column j = 1 or in the last row
              ! (i' = i - 1, j' = j + 1) until it stops either at first row i = 1 or
              ! the last column j = n_j.
    
    use mpi_f08                        
    implicit none
    type(MPI_Status) :: status
    integer :: my_rank, n_ranks, bottom_rank, top_rank, left_rank, right_rank
    integer :: n_i, n_j, n_cols, n_rows, remainder_j, remainder_i
    integer :: r, c, tmp
    integer :: height, width, r_i, c_i, r_f, c_f, height_c, ci, ri
    real, allocatable :: atmos_global(:,:), atmos_local(:,:)
    type(MPI_Datatype) :: col_g, row_g, row_r1, row_r2, row__g, row__r2, row__r1 
    integer(kind = MPI_ADDRESS_KIND) :: lb, real_extent

    ! Initialize the parallel mode

    call MPI_Init()
    call MPI_Comm_size(MPI_COMM_WORLD, n_ranks)
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank)
    
    ! Read the data and calculate the initial data.

    if (my_rank == 0) then
        print '(a)', "Grid size of the atmosphere?"
        read *, n_i, n_j ! Read the size of the global matrix.
        allocate(atmos_global(n_i, n_j), source = 0.0)
        print '(a)', "Define a 2D Cartesian topology:"
        read *, n_rows, n_cols ! Read the size of the topology.
    end if

    ! Send the topology to all the ranks with collective calls.
    call MPI_Bcast(n_rows, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
    call MPI_Bcast(n_cols, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)

    call MPI_Bcast(n_j, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
    call MPI_Bcast(n_i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)


    if (n_rows * n_cols /= n_ranks) then
        call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_TOPOLOGY)
        stop "Run this code with '-np 4'."
    end if
    
  
    
    
    remainder_i = mod(n_i, n_rows)
    remainder_j = mod(n_j, n_cols)
        
    height = (n_i - remainder_i)/n_rows
    width  = (n_j - remainder_j)/n_cols

     ! Create the topology. The order is the same as in Fortran arrays, from top
    ! to bottom, from left to right.

    tmp = my_rank
    r = mod(tmp, n_rows)
    tmp = (tmp - r) / n_rows
    c = mod(tmp, n_cols)

    ! Let have an extra column or row
    if (r < remainder_i) height = height + 1
    if (c < remainder_j) width = width + 1
    
    ! Define neighbor rank. If they fall out of the allowed ranges, make them
    ! MPI_PROC_NULL
    top_rank    = get_rank(r - 1, c)
    bottom_rank = get_rank(r + 1, c)
    left_rank   = get_rank(r, c - 1)
    right_rank  = get_rank(r, c + 1)
    
    
    ! Use global indices for the local matrix
    
    call global_index()


    ! Allocate the local matrix with the global indices.    
    allocate(atmos_local(r_i:r_f, c_i:c_f), source = 0.0)

    ! Initializate the atmosphere in each rank
    call init_atmos()

    ! Create the initial derived types
    call create_derived_type()

    ! Print the initial local matrix
    call print_local("Initial")
    ! Send the intial values to the global matrix and print them
    call send_to_global("Initial")
    ! Solve the atmosphere for each rank
    call solve()
    ! Print the final local matrix 
    call print_local("Final")
    ! Send the final values to the global matrix and print them
    call send_to_global("Final")
 

    
    if (my_rank == 0) then
        deallocate(atmos_global)
        call MPI_Type_free(row_r1)
        call MPI_Type_free(row_r2)
    end if
    deallocate(atmos_local)
    call MPI_Type_free(row_g)
    call MPI_Type_free(col_g)
    call MPI_Finalize()


contains
   function get_rank(r, c)
       ! Calculate the rank to the topology
       implicit none
       integer :: get_rank
       integer, intent(in) :: r, c

       if (0 <= c .and. c < n_cols .and. 0 <= r .and. r < n_rows) then
           get_rank = r + c * n_rows
       else
           get_rank = MPI_PROC_NULL
       end if
       return
   end function get_rank


   subroutine create_derived_type()

       implicit none
       integer :: w

       call MPI_Type_get_extent(MPI_REAL, lb, real_extent)

       call MPI_Type_vector(width, 1, height+1, MPI_REAL, row__g)
       call MPI_Type_create_resized(row__g, lb, real_extent, row_g)
       call MPI_Type_commit(row_g)

       call MPI_Type_contiguous(height+1, MPI_REAL, col_g)
       call MPI_Type_commit(col_g)

       if (my_rank == 0) then
           w = (n_j - remainder_j)/n_cols ! Initial width
           call MPI_Type_vector(w+1, 1, n_i, MPI_REAL, row__r2)
           call MPI_Type_create_resized(row__r2, lb, real_extent, row_r2)
           call MPI_Type_commit(row_r2)
           call MPI_Type_vector(w, 1, n_i, MPI_REAL, row__r1)
           call MPI_Type_create_resized(row__r1, lb, real_extent, row_r1)
           call MPI_Type_commit(row_r1)
       end if

   end subroutine create_derived_type



   subroutine global_index()
            
       call MPI_Recv(r_i, 1, MPI_INTEGER, top_rank, 0, MPI_COMM_WORLD, status)
       if (top_rank == MPI_PROC_NULL) r_i = (r*height)+1
       r_f = r_i + height
       call MPI_Send(r_f, 1, MPI_INTEGER, bottom_rank, 0, MPI_COMM_WORLD)
           

       call MPI_Recv(c_i, 1, MPI_INTEGER, left_rank, 0, MPI_COMM_WORLD, status)
       if (left_rank == MPI_PROC_NULL) c_i = (c*width)
       c_f = c_i + width
       call MPI_Send(c_f, 1, MPI_INTEGER, right_rank, 0, MPI_COMM_WORLD)


       if (remainder_i == 0) then
           r_i = (r * height) + 1
           r_f = r_i + height
       else if (remainder_j == 0) then
           c_i = (c * width)
           c_f = c_i + width
       end if  
   
   end subroutine global_index

        ! subroutines init_atmos, print_atmos, solve, and propagate
   subroutine init_atmos ! The first subroutine sets up a 2D array atmos,
                         !which we will call the atmosphere for simplicity
       implicit none
       integer :: i, j


       do j = c_i+1, c_f
           do i = r_i, r_f-1
               atmos_local(i, j) = 1.42857143 * i + 5.8823529 * j
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
           print '(*(es12.6, 1x))', atmos_global(i,:)
       end do
   end subroutine print_atmos ! It is called two times in the main program, 
                              ! before and after the solution

   subroutine print_local(name) ! Prints out the local altmosphere row-by-row
       implicit none
       character(len=*), intent(in) :: name
       integer :: rank, row
       
       do rank = 0, n_ranks - 1
           if (rank == my_rank) then
               if (rank == 0) print '(2a)', name, " rank values:"
               print '(a, i0)', "Rank ", my_rank
               do row = r_i, r_i + height
                   print '(*(es12.6, 1x))', atmos_local(row,:)
               end do
           end if
           call MPI_Barrier(MPI_COMM_WORLD)
       end do
   end subroutine print_local ! It is called two times in the main program,
                              ! before and after the solution

   subroutine send_to_global(name)
       implicit none
       character(len=*), intent(in) :: name
       integer :: rank
       
       if (my_rank == 0) then
         ! Send the local data to the global matrix

           atmos_global(r_i:r_f-1, c_i+1:c_f) = atmos_local(r_i:r_f-1, c_i+1:c_f)
            do rank = 1, n_ranks - 1

                tmp = rank
                r = mod(tmp, n_rows)
                tmp = (tmp - r)/ n_rows
                c = mod(tmp, n_cols)

                call MPI_Recv(ri, 1, MPI_INTEGER, rank, 0, MPI_COMM_WORLD, status)
                call MPI_Recv(ci, 1, MPI_INTEGER, rank, 0, MPI_COMM_WORLD, status)
                call MPI_Recv(height_c, 1, MPI_INTEGER, rank, 0, MPI_COMM_WORLD,&
                        status)

                if (c < remainder_j) then
                    call MPI_Recv(atmos_global(ri, ci+1), height_c, row_r2, rank, 0, &
                       MPI_COMM_WORLD, status)           
                else  
                    call MPI_Recv(atmos_global(ri, ci+1), height_c, row_r1, rank, 0, &
                       MPI_COMM_WORLD, status)
                end if

            
            end do
            call print_atmos(name)
         
        else
         
            call MPI_Send(r_i, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD)
            call MPI_Send(c_i, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD)
            call MPI_Send(height, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD)
            call MPI_Send(atmos_local(r_i,c_i+1), height, row_g, 0, 0,&
                 MPI_COMM_WORLD)

        end if

   end subroutine send_to_global


   subroutine solve() ! Calls the last subroutine propagate for all points of the
                      ! local array in the first column and in the last row
       implicit none
       integer :: i, j

       
       ! Update the ghost layers
       call update_ghost("recv")       
       
       ! Rays starting in the 1st column:
       j = c_i
       if (j == 0) j = j + 1
       do i = r_i, r_f
           if (i > n_i) exit
           call propagate(i, j)
       end do

       ! Rays starting in the last row:
       i = r_f
       if (i > n_i) i = i - 1
       do j = c_i+1, c_f
           if (j == 1) cycle
           call propagate(i, j)
       end do

       ! Update the ghost layers
       call update_ghost("send")

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
           if (next_i < r_i .or. next_j > c_f) exit
           atmos_local(next_i, next_j) = &
                   (atmos_local(crnt_i, crnt_j) + atmos_local(next_i, next_j)) / 2.0
           crnt_i = next_i
           crnt_j = next_j
       end do

   end subroutine propagate
 
   subroutine update_ghost(name) ! Update the ghost layers between ranks
       implicit none
       character(len=*), intent(in) :: name

       if (name=="send") then
           call MPI_Send(atmos_local(r_i, c_i+1), 1, row_g, top_rank, 0, MPI_COMM_WORLD)
           call MPI_Send(atmos_local(r_i, c_f), 1, col_g, right_rank, 1,&
              MPI_COMM_WORLD)
       else if (name=="recv") then
           call MPI_Recv(atmos_local(r_f, c_i+1), 1, row_g, bottom_rank, 0, &
               MPI_COMM_WORLD,status)
           call MPI_Recv(atmos_local(r_i,c_i), 1, col_g, left_rank, 1, &
               MPI_COMM_WORLD,status)   
       end if
   
   end subroutine update_ghost


       
end program rt_2d
