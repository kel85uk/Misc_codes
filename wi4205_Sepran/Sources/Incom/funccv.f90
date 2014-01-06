subroutine funccv ( icurve, t, x, y, z)
implicit none
integer icurve
double precision t, x, y, z, alpha, beta

alpha = 0.7
beta = 0.4

x = t
y = 1 - alpha*exp(-beta*t**2)
if ( t==5 ) then
	y=1
	end if
if ( t==-5 ) then
	y=1
	end if	

end