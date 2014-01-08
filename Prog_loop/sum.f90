! This is a simple program to reads in 2 numbers and prints the sum. 
! The sum is recorded in "results_sum.txt", created by initial_input.exe (init.f90),
! and overwrites the 0 with the sum value.

program sum
    implicit none
    integer :: a, b, s
    ! a - one of the two numbers to be added
    ! b - the other number in the sum
    ! s - the sum of a and b
    integer, parameter :: out_unit1 = 20
    integer, parameter :: out_unit2 = 21
    character(*), parameter :: fileplace = './/'

    open (unit=out_unit1, status="old", file=fileplace//'results_prod.txt', action="read")
    read (20, *) a, b
    close (out_unit1)
 
    s = a + b
    
    open (unit=out_unit2, status="old", file=fileplace//"results_sum.txt", action="write")
    write (out_unit2,*) s

end program sum
