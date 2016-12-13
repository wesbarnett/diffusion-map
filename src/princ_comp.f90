module princ_comp

    implicit none
    type, public :: princ_comp_type
        real(8), allocatable :: cov(:,:), contrib(:), cumul_contrib(:), evec(:,:), eval(:), data(:,:)
    contains
        procedure :: run => princ_comp_run
    end type

contains

    ! TODO clean up
    subroutine get_evect_sym(cov,evect,evalues)

        implicit none
        integer :: N, OK
        double precision, dimension(:), allocatable ::  work
        double precision, dimension(:), allocatable, intent(inout) ::  evalues
        double precision, dimension(:,:), allocatable, intent(inout) :: evect
        double precision, dimension(:,:), allocatable, intent(inout) :: cov
        double precision, dimension(:,:), allocatable :: A
        integer, allocatable, dimension(:) :: ORDER
        logical, allocatable, dimension(:) :: mask
        integer :: I
        integer :: LWORK, LDA
        character :: jobz = 'V', uplo = 'U'

        interface 
          subroutine dsyev(jobz,uplo,N,A,LDA,w,work,LWORK,INFO)
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: N, LDA, LWORK
            integer, intent(out) :: INFO
            double precision, intent(inout) :: A(LDA,*)
            double precision, intent(out) :: w(*), work(*)
          end subroutine
        end interface

        A = cov
        N = size(A,1)
        LWORK = 26*N
        LDA = N
        allocate(work(LWORK))
        allocate(evect(N, N))
        allocate(evalues(N))

!       Get eigenvectors and eigenvalues (LAPACK)
        call dsyev(jobz,uplo,N,A,LDA,evalues,work,LWORK,OK)

!       Sort the evalues in descending so we can use the top components for
!       later calcs
        allocate(ORDER(N))
        allocate(mask(N))
        mask = .true.
        do I = lbound(evalues,1), ubound(evalues,1)
            ORDER(I) = maxloc(evalues,1,mask)
            mask(ORDER(I)) = .false.
        end do

        evalues = evalues(ORDER)
        evect = A(:,ORDER)

    end subroutine get_evect_sym

    subroutine princ_comp_run(this, indata)

        integer :: ncomp, ndata, i, j
        real(8), allocatable, intent(inout) :: indata(:,:)
        class(princ_comp_type), intent(inout) :: this

        ndata = size(indata,1)
        ncomp = size(indata,2)

        ! mean adjust
        do i = 1, ndata
            indata(i,:) = indata(i,:) - sum(indata(i,:))/dble(ndata)
        end do

        allocate(this%cov(ncomp, ncomp))
        allocate(this%contrib(ncomp))
        allocate(this%cumul_contrib(ncomp))
        allocate(this%data(ndata,ncomp))

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

        this%data = matmul(transpose(this%evec), indata)

    end subroutine princ_comp_run

end module princ_comp
