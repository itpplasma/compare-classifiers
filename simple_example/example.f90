
program example
    use simple, only: Tracer, Field, CutDetector

    implicit none

    ! Constants
    real(8), parameter :: pi = 3.14159265358979d0
    real(8), parameter :: c = 2.9979d10
    real(8), parameter :: e_charge = 4.8032d-10
    real(8), parameter :: e_mass = 9.1094d-28
    real(8), parameter :: p_mass = 1.6726d-24
    real(8), parameter :: ev = 1.6022d-12

    ! Initial parameters
    real(8) :: s0, th0, ph0, lam
    real(8) :: bmod_ref, v0, rlarm, tau
    real(8) :: Z_charge, m_mass, E_kin
    real(8) :: s, theta, varphi, theta_vmec, varphi_vmec

    integer :: ierr, kt, ntimstep

    type(Tracer) :: tracy

    ! Initialize parameters here

    ! Initialize Tracer
    call init_tracer(tracy)

    ! Initialize field and parameters
    call init_field("wout.nc", 5, 5, 3, 1)
    call init_params(Z_charge, m_mass, E_kin, 256, 1, 1e-13)

    ! Define initial coordinates here

    ! Convert coordinates
    call can_to_vmec(s, theta, varphi, theta_vmec, varphi_vmec)
    call vmec_to_cyl(s, theta_vmec, varphi_vmec, R, Z)

    ! Time integration
    integer :: ntimstep
    ! Define ntimstep and other variables

    ! Initialize integrator
    call init_integrator(z)

    do kt = 1, ntimstep
        ierr = timestep_sympl_z(f, z)
        ! Store or process z here
    end do

end program example
