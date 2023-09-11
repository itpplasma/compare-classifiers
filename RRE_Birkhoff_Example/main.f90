program main
    use birkhoff_procedures

    implicit none
 
    integer(int64), parameter :: nmax=70
    integer(int64), parameter :: d=2
    integer(int64), parameter :: K=21

    
    real(real64), dimension(d, nmax) :: seq
    real(real64), dimension(d) :: ave = 0
    real(real64), dimension(d) :: ave2 = 0
    real(real64), dimension(2*K+1) :: c = 0
    real(real64) :: res


    call fill_test_sequence(seq)
    call birkhoff_average(seq, ave)
    call birkhoff_RRE(seq, K, c, ave2, res)
    
    print *, "Using a test sequence of length "
    print *, "N = ", nmax
    print *, " and number of filter unkowns " 
    print *, "K = ", K
    print *, "The filter length is 2K+1 = ", 2*K+1
    print *, "The least squares system is size d*(N-2K-1) x K = ", d*(nmax-(2*K)-1), " x ", K
    print *, " "

    print *, "The averages of e^cos(n) and e^sin(n) are" 
    print *, ave(:)
    print *, "The averages from RRE are"
    print *, ave2(:)
    print *, "for nmax = ", nmax
    print *, " "
    print *, "The difference of the two WBA averages is ", abs(ave(2)-ave(1))
    print *, "The difference of the two RRE averages is ", abs(ave2(2)-ave2(1))
    print *, "The norm of the difference of WBA from RRE is ", sqrt((ave(1)-ave2(1))**2 + (ave(2)-ave2(2))**2)
    print *, "The RRE residual is ", res

    
 end program