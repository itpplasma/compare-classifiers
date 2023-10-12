program example
    use util, only: pi, c, e_charge, p_mass, ev
    use simple
    use cut_detector, only: CutDetector, init_cut_detector => init, trace_to_cut

    implicit none

    integer, parameter :: ntimstep = 10000
    integer, parameter :: ncut = 500

    ! Initial parameters
    real(8) :: s0, th0, ph0, lam
    real(8) :: bmod_ref
    real(8) :: E_kin
    integer :: Z_charge, m_mass

    ! To store orbit
    real(8) :: z(5)
    real(8) :: zs(5, ntimstep)

    ! To store cuts (first 5 entries are z, last one is normalized J_parallel)
    real(8) :: var_cut(6)
    real(8) :: var_tips(6, ncut)

    ! Converted VMEC and cylinder coordinates
    real(8) :: theta_vmec, varphi_vmec, Rcyl, Zcyl

    integer :: ierr, k, cut_type

    type(Tracer) :: tracy
    type(CutDetector) :: cutty

    ! Initialize parameters here
    s0 = 0.6
    th0 = 0.0d0
    ph0 = 0.25d0*pi
    lam = -0.4d0
    !lam = 0.8   ! passing regular

    bmod_ref = 5.0d4
    Z_charge = 2
    m_mass = 4
    E_kin = 3.5d6

    ! Initialize field and parameters
    call init_field(tracy, "wout.nc", 5, 5, 3, 1)
    call init_params(tracy, Z_charge, m_mass, E_kin, 256, 1, 1d-13)

    ! Example of coordinate conversion (there exists also vmec_to_can)
    call can_to_vmec(s0, th0, ph0, theta_vmec, varphi_vmec)
    call vmec_to_cyl(s0, theta_vmec, varphi_vmec, Rcyl, Zcyl)

    print *, "s, theta, varphi = ", s0, th0, ph0
    print *, "s, theta_vmec, varphi_vmec = ", s0, theta_vmec, varphi_vmec
    print *, "Rcyl, Zcyl = ", Rcyl, Zcyl

    ! -----------------------
    ! Example of orbit tracer
    ! -----------------------

    ! Initialize integrator
    z(1) = s0
    z(2) = th0
    z(3) = ph0
    z(4) = 1.0  ! Normalized momentum mv
    z(5) = lam  ! lambda = v_parallel/v

    call init_integrator(tracy, z)

    zs(:, 1) = z
    do k = 1, ntimstep-1
        call timestep_sympl_z(tracy%si, tracy%f, z, ierr)
        zs(:, k+1) = z
    end do

    ! Print out orbit
    open(101, file="orbit.dat")
    do k = 1, ntimstep
        write(101,*) zs(:, k)
    end do
    close(101)

    ! -----------------------
    ! Example of cut detector
    ! -----------------------

    ! Initialize integrator
    z(1) = s0
    z(2) = th0
    z(3) = ph0
    z(4) = 1.0  ! Normalized momentum mv
    z(5) = lam  ! lambda = v_parallel/v

    call init_integrator(tracy, z)

    call init_cut_detector(cutty, tracy%fper, z)
    ierr = 0
    var_cut = 0d0
    var_tips = 0d0
    do k = 1, ncut
        call trace_to_cut(cutty, tracy%si, tracy%f, z, var_cut, cut_type, ierr)
        var_tips(:, k) = var_cut
    end do

    ! Print out cut
    open(101, file="cut.dat")
    do k = 1, ncut
        write(101,*) var_tips(:, k)
    end do
    close(101)
end program example
