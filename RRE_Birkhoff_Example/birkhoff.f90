module birkhoff_procedures
    use iso_fortran_env, only:  int64, real64

    implicit none
contains 
    ! The standard weight for weighted Birkhoff averages
    function birkhoff_w(n, nmax) result(w)
        integer(int64), intent(in) :: n, nmax
        real(real64)                :: w
        real(real64)                :: x

        x = real(n, real64)/real(nmax + 1, real64)
        w = exp(-1/(x*(1-x)))
    end function

    ! Fills a sequence of the form e^sin(n), e^cos(n)
    subroutine fill_test_sequence(seq)
        real(real64), intent(inout)   :: seq(:, :)
        integer(int64)                :: n

        do n = 1, ubound(seq, 2)
            seq(1, n) = exp(sin(real(n, real64)))
            seq(2, n) = exp(cos(real(n, real64)))
        end do
    end subroutine

    
    ! Compute the birkhoff average of a vector-valued sequence
    ! Assumes that seq takes the matrix form [v1, v2, v3,...,vN]
    ! where each vi is length d
    !
    ! Outputs the length-d average
    subroutine birkhoff_average(seq, ave)
        real(real64), intent(in)  :: seq(:, :)
        real(real64), intent(out) :: ave(:)

        real(real64) :: total_weight = 0
        real(real64) :: w

        integer(int64) :: m, mmax
        
        ave = 0
        mmax = ubound(seq, 2)
        do m = 1, mmax
            w = birkhoff_w(m, mmax)
            total_weight = total_weight + w
            ave = ave + w * seq(:, m)
        end do

        ave = ave / total_weight
    end subroutine

    ! Solves RRE on a vector-valued sequence with filter length 2K+1. The
    ! other dimensions are determined by seq = d x N.
    !
    ! The size of the system solved is d T x K, where N = 2K + T + 1
    ! So, you need dT = d(N-2K-1) >= K for a full system. Typically, it's best
    ! to choose something like dT = 1.5K
    ! 
    ! Input: 
    ! - seq: The vector-valued sequence, stored in a matrix [v1, v2, v3,...,vN]
    !        where each vi is length d
    ! - K: Number of filter unkowns
    !
    ! Output:
    ! - c: The output filter of length 2K+1
    ! - ave: length d sequence average, as computed by applying c to the 
    !        first 2K+1 entries of seq
    ! - res: The residual of the least-squares problem
    subroutine birkhoff_RRE(seq, K, c, ave, res)
        real(real64), intent(in)   :: seq(:, :)
        integer(int64), intent(in) :: K

        real(real64), intent(out) :: c(:)
        real(real64), intent(out) :: ave(:)
        real(real64), intent(out) :: res

        real(real64), allocatable :: A(:, :), A2(:, :), b(:), b2(:), work(:), lsqrerr(:)
        real(real64), allocatable :: U(:, :)
        real(real64) :: wii, work_query(1), wsum


        integer(int64) :: d, N, T, ii, jj, kk, info, lwork, Nrhs
        d = ubound(seq, 1)
        N = ubound(seq, 2)
        T = N - 2*K - 1

        allocate(A(T*d, K))
        allocate(A2(T*d, K))
        allocate(b(T*d))
        allocate(b2(T*d))

        ! Find the difference sequence
        allocate(U(d, N-1))
        do ii = 1, N-1
            U(:, ii) = seq(:, ii+1) - seq(:, ii)
        end do

        ! Assemble the least squares system
        wsum = 0.0
        do ii = 1, T
            wii = sqrt(birkhoff_w(ii, T))
            wsum = wsum + birkhoff_w(ii, T)
            b(d*(ii-1)+1:d*ii) = - wii*U(:, ii+K)
            do jj = 1, K
                A(d*(ii-1)+1:d*ii, jj) = wii*(U(:, ii+jj-1)   - U(:, ii+jj) &
                                            - U(:, ii+2*K-jj) + U(:, ii+2*K+1-jj))
            end do
        end do
        b2(:) = b(:)
        A2(:, :) = A(:, :)


        ! Solve the least-squares system
        Nrhs = 1
        lwork = -1
        call dgels ('N', d*T, K, Nrhs, A, d*T, b, d*T, work_query, lwork, info)
        lwork = int(work_query(1), int64)
        allocate(work(lwork))
        call dgels ('N', d*T, K, Nrhs, A, d*T, b, d*T, work, lwork, info)
        
        ! Get the filter from the constrained least squares solution
        c(1) = b(1)
        c(2*K+1) = c(1)
        do ii = 2,K
            c(ii) = b(ii) - b(ii-1)
            c(2*K+2-ii) = c(ii)
        end do
        c(K+1) = 1 - 2*b(K)

        ! Get the Birkhoff average
        call dgemm('N', 'N', d, 1, 2*K+1, real(1.0, real64), seq, d, c, 2*K+1, real(0.0, real64), ave, d)


        ! Get the residual
        allocate(lsqrerr(T*d))
        call dgemm('N', 'N', T*d, 1, K, real(1.0, real64), A2, T*d, b(1:K), K, real(0.0, real64), lsqrerr, T*d)
        lsqrerr = lsqrerr - b2
        res = real(0.0, real64)
        do ii = 1,T*d
            res = res + lsqrerr(ii)**2
        end do
        res = sqrt(res / wsum)
    end subroutine
end module
