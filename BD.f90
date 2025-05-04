! Numerical resolution of the dimensionless BD field equations in vacuum 
! for a static, spherically symmetric spacetime.
! Joel Argudo Panes, spring 2025.
! TFG.
! Universitat de Barcelona.

PROGRAM BD

    IMPLICIT NONE
    
    INTEGER:: i, j, m, ndim = 3, N, N_var
    DOUBLE PRECISION:: A0, B0_inv, y0,dy0, yini(3), xini
    DOUBLE PRECISION:: derivatives, x0, xf, k, w, h, error, x_var
    DOUBLE PRECISION:: variation,variations(3), y_previous(3), y_out(3), y_in(3)
    COMMON / PARAMETERS / w, k, x0
    EXTERNAL derivatives 
    
    ! End point.
    xf = 0.5d0
    ! Point where we study the stability of the solution.
    x_var = 1.1d0
    ! Integration step.
    h = -0.1d0
    ! Maximum desired variation of the solutions at x_var with respect to the 
    ! previous iteration.
    error = 1d-5
    
    DO i = 1, 300 ! w (BD parameter) loop.
    
        j = 0
    
        IF ((i.GT.10).AND.(mod(i,10).NE.0)) THEN 
            CYCLE
        ELSE
            w = dble(i)
        ENDIF 
    
        y_previous = [0.d0, 0.d0, 0.d0]  ! Solutions in the previous iteration.
        variation = 1.d0  ! Maximum variation of the solutions at x_var
    
        ! x0 (initial point) loop.
        DO WHILE (variation.GE.error) 
    
            x0 = 10.1d0 + dble(j)
            j = j + 1
    
            ! Number of steps.
            N = int((xf-x0)/h)
            N_var = int((x_var-x0)/h)
    
       ! Initial conditions (weak field solution at x0).
        A0 = -1.d0/x0+1.d0 
        B0_inv = 1.d0/(1.d0+1.d0/x0*(1.d0+w)/(2.d0+w))
        y0 = dlog((2.d0*w+4.d0)/(2.d0*w+3)+1.d0/((2.d0*w+3.d0)*x0))
        dy0 = (-1.d0/((2.d0*w+3.d0)*x0**2))/((2.d0*w+4.d0)/(2.d0*w+3)+1.d0/((2.d0*w+3.d0)*x0))

        k = x0**2*dsqrt(A0*B0_inv)*dy0*dexp(y0)

        yini = [y0,A0,B0_inv]
    
        ! We solve the equations with RK4 and compare the solutions at x_var with 
        ! the previous iteration.
        y_in = yini
        xini = x0
    
        DO m = 1, N_var
    
            CALL RK4_step(ndim,xini,y_in,h,y_out,derivatives)
            xini = xini + h
            y_in = y_out
    
        ENDDO 
    
        variations = [dabs(y_out(1)-y_previous(1)), dabs(y_out(2)-y_previous(2)),&
                     dabs(y_out(3)-y_previous(3))]
    
        variation = max(variations(1), variations(2),variations(3))
        y_previous = y_out
    
        ENDDO
    
        ! We solve the equations with RK4 with the chosen x0 and write the solutions
        ! in a file.
        CALL RK4(ndim,x0,xf,yini,h,derivatives)
    
    ENDDO
    
    
    END PROGRAM 
    
    ! Subroutine that takes in the values of x and the functions y (A(x), 1/B(x), y(x))
    ! and returns their derivatives f.
    
    SUBROUTINE derivatives(x,y,f)
        IMPLICIT NONE
        DOUBLE PRECISION:: x, y(3), f(3)
        DOUBLE PRECISION:: w, k
        COMMON / PARAMETERS / w,k 
    
        f(1) = k*dexp(-y(1))/x**2/dsqrt(dabs(y(3)*y(2)))
        f(2) = (-(1.d0-1.d0/y(3))/x**2+w/2.d0*(f(1))**2-2.d0/x*f(1))/(1.d0/(y(2)*x)+f(1)/(2.d0*y(2)))
        f(3) = x*y(3)*(-(1.d0-1.d0/y(3))/x**2-w/2.d0*(f(1))**2+f(2)/(2.d0*y(2))*f(1))
    
        RETURN 
    END
    
    ! Subroutine that solves a set of ordinary differential equations with the RK4 algorithm.
    
    ! IN: ndim (number of equations), xini and yini (initial conditions), h (step),
    !     N (number of steps), derivatives (subroutine that computes the derivatives).
    ! OUT: the values of x, y, w and x0 are written in a file.
    
    SUBROUTINE RK4(ndim,xini,xfin,yini,step,derivatives)
        IMPLICIT NONE
        INTEGER:: ndim, N, i
        DOUBLE PRECISION:: xini, xfin, yini(ndim), h,step, step2
        DOUBLE PRECISION:: x, xlim, y(ndim)
        DOUBLE PRECISION:: k1(ndim), k2(ndim), k3(ndim), k4(ndim)
        DOUBLE PRECISION:: w, k, x0
        COMMON / PARAMETERS / w,k,x0
    
        x = xini
        y = yini

    
        step2 = step/100.d0 ! Integration step for x < xlim.

        xlim = 2.d0

        N = int(-(xini-xlim)/step) + int(-(xlim+step-xfin)/step2) + 1

        OPEN(11, file = "BD_solutions.dat", position = "append")

        WRITE(11,"(/,/)")
        WRITE(11,*) "# x / y / A / B^{-1} / ω / x0"
        WRITE(11,"(6f42.12)") xini, yini, w, x0
    
        DO i = 1, N

            IF (x.LE.xlim) THEN 
                h = step2
            ELSE
                h = step
            ENDIF 

            CALL derivatives(x,y,k1)
            
            CALL derivatives(x+h/2.d0,y+h/2.d0*k1,k2)
    
            CALL derivatives(x+h/2.d0,y+h/2.d0*k2,k3)
    
            CALL derivatives(x+h,y+h*k3,k4)
    
            x = x + h
    
            y = y + h/6.d0*(k1 + 2.d0*k2 + 2.d0*k3 + k4)
            
            WRITE(11,"(6f42.12)") x, y, w, x0
    
        ENDDO
    
        CLOSE(11)
    
        RETURN
    END
    
    ! Subroutine that computes one integration step with the RK4 algorithm.
    
    ! In: ndim (number of equations), xini and yini (initial conditions),
    !     h (step), derivatives (subroutine that computes the derivatives).
    ! OUT: y values after one step.
    
    SUBROUTINE RK4_step(ndim,xini,yini,h,y,derivatives)
        IMPLICIT NONE
        INTEGER:: ndim
        DOUBLE PRECISION:: xini, yini(ndim), h, x, y(ndim)
        DOUBLE PRECISION:: k1(ndim), k2(ndim), k3(ndim), k4(ndim)
    
        x = xini
        y = yini
    
            CALL derivatives(x,y,k1)
            
            CALL derivatives(x+h/2.d0,y+h/2.d0*k1,k2)
    
            CALL derivatives(x+h/2.d0,y+h/2.d0*k2,k3)
    
            CALL derivatives(x+h,y+h*k3,k4)
    
            y = y + h/6.d0*(k1 + 2.d0*k2 + 2.d0*k3 + k4)
    
        RETURN
    END