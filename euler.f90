!-------------------------------------------------------------------------!
! This is the code to calculate the radial wavefunctions of Hydrogen atom !
! by using Euler method. It doesn't find energy eigenvalues. It           !
! calculates the radial wavefunctions for given values of n and l. The    !
! value of n assigns the energy.                                          !
!                                                                         !
! In this code, u(r) is replaced by psi and R(r) is replaced by psi_plot  !
! v(r) is replaced by dpsi                                                !
!-------------------------------------------------------------------------!

program numerov_method
   implicit none

   !parameters
   !----------------------------------------------------------------------!
   real(8), parameter :: h_b = 1.054571800d-34  ! Planck constant
   real(8), parameter :: m_e = 9.10938356d-31   ! electron mass
   real(8), parameter :: e_c = 1.6021766208d-19 ! electronic charge
   real(8), parameter :: e_0 = 8.8541878128d-12 ! vacuum permitivitty
   real(8), parameter :: pi  = 4.d0 * atan(1.d0)
   real(8), parameter :: A   = ( 2.0d0 * m_e ) / ( h_b * h_b )
   real(8), parameter :: B   = e_c * e_c / ( 4.0d0 * pi * e_0 )
   real(8), parameter :: a_0 =  5.2917721067d-11 ! BOhr radius

   !----------------------------------------------------------------------!

   real(8) :: r_min, r_max, h, h2, h12, E, norm, fac, psi_b
   integer :: n, l, nr, i, match, j
   real(8), allocatable :: psi(:), r(:), f(:), dpsi(:), psi_plot(:)
   real(8), external :: Energy, psi_infinity, psi_zero, test_func1
   character (len = 50 ) :: filename
   !======================================================================!

   10 print*, "# Give the values of n and l"
   read(*,*) n, l

   if ( l < 0 .or. l > n - 1) then
      print*, ""
      print*, "n and l must satisfy n = 0, 1, 2..... and l = 0, 1, 2...n-1"
      print*, "Input the correct values"
      print*, " "
      goto 10
   end if        
   
   r_min = 1.0d-105
   r_max = 300.d0 * a_0
   
   print*,"Input the number of points to calculate wave function values at"
   read(*,*) nr

   h = ( r_max - r_min ) / nr

   allocate(psi(0:nr+1))
   allocate(psi_plot(0:nr+1))
   allocate(r(0:nr+1))
   allocate(f(0:nr+1))
   allocate(dpsi(0:nr))

   print*, "Calculating energy"
   print*, "------------------"
   E = Energy(n)
   print*, "Energy for n =", n, "is", E
        
   !----------------------------------------------------------------------!
   !                                                                      !
   ! The trick is to calculate psi at two very close points in            !
   ! r ~ infinity and r ~ 0 region. Then use the values of psi to calcu   !
   ! -late the derivative of psi at r(0) and r(nr) using the basic defin  !
   ! -ition of derivative. Then use the Euler method to solve second order!
   ! ordinary differential equation.                                      !
   !                                                                      !
   !----------------------------------------------------------------------!

   do i = 0, nr + 1
      r(i) = r_min + dble(i) * h
      f(i) = A * ( l * ( l + 1 ) / ( A * r(i) * r(i)) - B / r(i) - E )
   end do
   
   !-----------------------------------------!
   ! To calculate the psi and derivative dpsi!
   !-----------------------------------------!

   psi(nr + 1) = psi_infinity(E,r(nr+1))
   
   psi(nr) = psi_infinity( E, r(nr) )
   dpsi(nr) = ( psi(nr+1) - psi(nr) ) / h

   psi(nr-1) = psi(nr) - h * dpsi(nr)
   dpsi(nr-1) = dpsi(nr) - h * f(nr) * psi(nr-1)

   infinity_to_zero : do
      do i = nr - 2, 0, -1
         psi(i)  = psi(i+1) - h * dpsi(i+1)
         dpsi(i) = dpsi(i+1) - h * f(i+1) * psi(i)
         
         if ( abs(psi(i)) <= abs(psi(i+1)) ) exit infinity_to_zero
      end do
   end do infinity_to_zero

   match = i            ! r(i) = r(match) is the point psi_infinity fails

   psi_b = psi(match)      ! value of psi at boundary in infinity region
 
   !------------------------------------------!
   ! To calculate psi and dpsi at zero region !
   !------------------------------------------!

   psi(1) = psi_zero( E, r(1), l )
   
   psi(0) = psi_zero( E, r(0), l )
   dpsi(0) = ( psi(1) - psi(0) ) / h

   dpsi(1) = dpsi(0) + h * f(0) * psi(1)
   
   do i = 2, match
      psi(i) = psi(i-1) + h * dpsi(i-1)
      dpsi(i) = dpsi(i-1) + h * f(i-1) * psi(i)
   end do

   !----------------------------------------------------------------------!   
   !                                                                      !
   !   Now we scale the calculated wave functions at zero region to make  !
   !   the wavefunction continuous                                        !
   !----------------------------------------------------------------------!

   fac = psi_b / psi(match)
   print*, "fac =", fac
   psi(0:match) = psi(0:match) * fac

   !-----------------!
   !  Normalization  !
   !-----------------!

   norm = 1.0d0 / sqrt ( h * sum(psi*psi) )
   psi = psi * norm

   !----------------------------------------------------------------------!
   ! Now actually the psi we calculated are u in Griffith's notation and  !
   ! the relation is:                                                     !
   !                                                                      !
   !                            u(r) = r * psi(r)                         !
   !                                                                      !
   ! So divide our calculated psi(i)'s by r(i)'s to find actual solutions !
   ! to radial part of Schrodinger equation for H-atom. Then we multiply  !
   ! them by a_0**1.5 to produce plots similar to in Griffiths' QM.       !
   !                                                                      !
   !----------------------------------------------------------------------!

   do i = 0, nr
      psi_plot(i) = psi(i) * a_0**(1.5d0)/r(i)
   end do

   if ( abs( maxval(psi_plot)) .lt. abs( minval(psi_plot)) ) then 
      psi_plot = psi_plot * (-1.0d0)
   else
      psi_plot = psi_plot
   end if   

   print*, "Give the name of file to save the data in."
   print*, "If you input 'euler.dat', plot will be generated automatically"
   read(*,*) filename   
   !----------------------------------------------------------------------!
   open( 3, file = filename )
   write(3,'( i0,2x,f0.10,1x,f0.10)') (i,r(i)/a_0,psi_plot(i),i=0,nr) 
   !----------------------------------------------------------------------!

   deallocate(psi)
   deallocate(dpsi)
   deallocate(r)
   deallocate(f)
   deallocate(psi_plot)
   
   !----------------------------------------------------------------------!

   call system('gnuplot -p euler.plt')
   
end program numerov_method      

!=========================================================================!

real(8) function Energy(n)
   implicit none
   !parameters
   !----------------------------------------------------------------------!
   real(8), parameter :: h_b = 1.054571800d-34
   real(8), parameter :: m_e = 9.10938356d-31
   real(8), parameter :: e_c = 1.6021766208d-19
   real(8), parameter :: e_0 = 8.8541878128d-12
   real(8), parameter :: pi  = 4.d0 * atan(1.d0)
   !----------------------------------------------------------------------!
   integer :: n

   Energy = - m_e * e_c**4 / ( 8.d0 * ( pi * e_0 * h_b * 2.d0 * n )**2)

end function Energy        

!=========================================================================!

real(8) function psi_infinity(E, r)
   implicit none
   real(8), parameter :: h_b = 1.054571800d-34
   real(8), parameter :: m_e = 9.10938356d-31
   real(8) :: E, r

   psi_infinity = exp( - sqrt(- 2.0d0 * m_e * E) * r / h_b )

end function psi_infinity   

!=========================================================================!

real(8) function test_func1(r)
   implicit none
   real(8) :: r
   real(8), parameter :: a =  5.2917721067d-11
   
   test_func1 = 2.0d0 * r * a**(-1.5d0) / exp( r / a )
end function test_func1

!=========================================================================!

real(8) function psi_zero(E,r,l)
   implicit none
   real(8), parameter :: h_b = 1.054571800d-34
   real(8), parameter :: m_e = 9.10938356d-31
   real(8) :: E, r
   integer :: l

   psi_zero = ( sqrt( - 2.0d0 * m_e * E ) * r / h_b )**(l+1)
end function psi_zero
