program bissection_main
    implicit none
    real:: a,b,racine
    double precision:: tolerance
    a=0.0
    b=2.0
    tolerance=1d-5

    call bissection(f,a,b,racine)
    print *, "La racine trouvee est ",racine



 contains
    real function f(x):
        real, intent(in):: x
        
        f=x**2- 2
    end function f


end program bissection_main


subroutine bissection(f,a,b,tolerance,racine)
    real, intent(inout):: a,b
    real, intent(out):: racine
    double precision, intent(in) :: tolerance
    real:: m
    logical:: go

    go=.true.

    if go then
        m=(a+b)/2
        if (f(a)*f(m) .lt. 0):
            b=m
        else 
            a=m
        end if 

        go= (abs(b-a) .lt. tolerance .or. f(m) .lt. tolerance)


    end if 

end subroutine bissection