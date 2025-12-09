program bissection_main
    implicit none
    real :: a, b, racine
    integer :: statut
    
    a = 0.0
    b = 2.0
    
    call bissection_method(f, a, b, 1e-6, 100, racine, statut)
    
    if (statut == 0) then
        print '(A, F12.8)', "Racine trouvée: ", racine
    else
        print *, "Erreur dans la recherche de la racine"
    end if

contains

    real function f(x)
        real, intent(in) :: x
        f = x**2 - 2.0
    end function f

end program bissection_main

subroutine bissection_method(func, a, b, tolerance, max_iter, racine, statut)
    implicit none
    real, intent(inout) :: a, b
    real, intent(in) :: tolerance
    integer, intent(in) :: max_iter
    real, intent(out) :: racine
    integer, intent(out) :: statut
    real, external :: func
    
    real :: fa, fb, fc, c
    integer :: iteration
    
    fa = func(a)
    fb = func(b)
    statut = -1
    
    if (fa * fb > 0.0) then
        print *, "Erreur: f(a) et f(b) doivent avoir des signes opposés"
        return
    end if
    
    do iteration = 1, max_iter
        c = (a + b) / 2.0
        fc = func(c)
        
        if (abs(fc) < tolerance .or. (b - a) / 2.0 < tolerance) then
            racine = c
            statut = 0
            return
        end if
        
        if (fa * fc < 0.0) then
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end if
    end do
    
    print *, "Maximum d'itérations atteint"
    racine = c
    statut = 1
    
end subroutine bissection_method