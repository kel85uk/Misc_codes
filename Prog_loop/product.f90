! This program reads in the sum value in "results_sum.txt", and assigns a variable x.
! It then multiplies x by 2, and calls it y. x and y gets recorded in "results_prod.txt"

program product
    implicit none
    integer :: x, y
    integer, parameter :: out_unit1 = 20
    integer, parameter :: out_unit2 = 10
    character(*), parameter :: fileplace = './/'

    open (unit=out_unit1, status="old", file=fileplace//'results_sum.txt', action="read")
    read (20, *) x
    close (out_unit1)
 
    y = x * 2
    
    print *, 'The input number is ', x
    print *, 'Integer y is equal to ', y
    
    open (unit=out_unit2, file="results_prod.txt", action="write", status="replace")
    write (out_unit2,*) x, y
end program product
    
    
