module princ_comp

    implicit none
    type, public :: princ_comp_type
        real(8), allocatable :: cov(:,:), contrib(:), cumul_contrib(:), evec(:,:), eval(:), x(:,:)
    contains
        procedure :: run => princ_comp_run
    end type

contains

    subroutine get_evect_sym(cov, evect, evalues)

        implicit none
        character :: jobz = 'V', uplo = 'U'
        real(8), allocatable ::  work(:), a(:,:)
        real(8), allocatable, intent(inout) ::  evalues(:), evect(:,:), cov(:,:)
        integer, allocatable :: order(:)
        integer :: i, lwork, lda, n, info
        logical, allocatable :: mask(:)

        interface 
          subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            integer, intent(out) :: info
            real(8), intent(inout) :: a(lda,*)
            real(8), intent(out) :: w(*), work(*)
          end subroutine
        end interface

        a = cov
        n = size(a,1)
        lwork = 26*n
        lda = n
        allocate(work(lwork), evect(n,n), evalues(n))

        ! Get eigenvectors and eigenvalues (LAPACK)
        ! matrix "a" is returned as the eigenvectors
        call dsyev(jobz, uplo, n, a, lda, evalues, work, lwork, info)
        if(info .ne. 0) then
            write(*,*) "ERROR: Eigenvector calculation failed."
        end if

        ! Sort the evalues in descending so we can use the top components for later calcs
        allocate(order(N), mask(n))
        mask = .true.
        do i = lbound(evalues,1), ubound(evalues,1)
            order(i) = maxloc(evalues,1,mask)
            mask(order(i)) = .false.
        end do

        evalues = evalues(order)
        evect = a(:,order)

    end subroutine get_evect_sym

    subroutine princ_comp_run(this, indata)

        integer :: ncomp, ndata, i, j, u
        real(8), allocatable, intent(inout) :: indata(:,:)
        class(princ_comp_type), intent(inout) :: this
character (len=1024) :: format_string
character (len=32) :: n_char

        ncomp = size(indata,1)
        ndata = size(indata,2)

        ! mean adjust
        do i = 1, ncomp
            indata(i,:) = indata(i,:) - sum(indata(i,:))/dble(ndata)
        end do

open(newunit=u, file="data.dat")
write(n_char,'(i0)') ncomp
format_string = "("//trim(n_char)//"f12.6)"
write(U, format_string) indata
close(u)


        allocate(this%cov(ncomp, ncomp))
        allocate(this%contrib(ncomp))
        allocate(this%cumul_contrib(ncomp))

        ! get covariance matrix
        do i = 1, ncomp
            do j = 1, ncomp
                this%cov(i,j) = dot_product(indata(j,:), indata(i,:))
                this%cov(j,i) = this%cov(i,j)
            end do
        end do

        this%cov = this%cov / dble(ndata - 1)

        call get_evect_sym(this%cov, this%evec, this%eval)

        this%contrib = this%eval / sum(this%eval)

        do i = 1, ncomp
            this%cumul_contrib(i) = sum(this%contrib(1:i))
        end do

        this%x = matmul(transpose(this%evec), indata)

    end subroutine princ_comp_run

end module princ_comp
