! This program asks the user to input 2 values, and records them in
! "results_prod.txt", which will be read by sum.exe for further processing.
! It also creates "results_sum.txt" and writes a 0 in it.

program init
    implicit none
    integer :: x1, x2
    integer, parameter :: out_unit1 = 10
    integer, parameter :: out_unit2 = 11
    character(*), parameter :: fileplace = './/'
    
    print *, 'Type in 2 integers separated by a comma or space:'
    read *, x1, x2
    
    open (unit=out_unit1, file=fileplace//"results_prod.txt", action="write", status="new")
    write (out_unit1,*) x1, x2
    
    open (unit=out_unit2, file=fileplace//"results_sum.txt", action="write", status="new")
    write (out_unit2,*) 0
end program init
