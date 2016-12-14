module diffusion_map

    implicit none
    type, public :: diffusion_map_type
        real(8), dimension(:,:), allocatable :: evec, map
        real(8), dimension(:), allocatable :: eval
    contains
        procedure :: run => diffusion_map_run
    end type

contains

    subroutine get_evect(p, evect, evalues)

        ! SLOW
        implicit none
        integer :: n, info, i, lwork, lda, ldvl, ldvr
        real(8), dimension(:), allocatable ::  wr, wi, work
        real(8), dimension(:,:), allocatable :: a, vl, vr
        real(8), dimension(:), allocatable, intent(inout) ::  evalues(:), evect(:,:), p(:,:)
        integer, allocatable :: order(:)
        logical, allocatable :: mask(:)
        character :: jobvl = 'N', jobvr = 'V'

        interface 
          subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
            character, intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            integer, intent(out) :: info
            real(8), intent(inout) :: a(lda,*)
            real(8), intent(out) :: wr(*), wi(*), vl(ldvl,*), vr(ldvr,*), work(*)
          end subroutine
        end interface

        a = p
        n = size(a,1)
        lwork = 26*n
        lda = n
        ldvl = n
        ldvr = n
        allocate(work(lwork), vl(ldvl, n), vr(ldvr, n), wr(n), wi(n))

!       Get eigenvectors and eigenvalues (LAPACK)
        call dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
        if(info .ne. 0) then
            write(*,*) "ERROR: Eigenvector calculation failed."
        end if

!       Sort the evalues in descending so we can use the top components for
!       later calcs
        allocate(order(n), mask(n))
        mask = .true.
        do i = lbound(wr,1), ubound(wr,1)
            order(i) = maxloc(wr,1,mask)
            mask(order(i)) = .false.
        end do

        allocate(evect(ldvr, n), evalues(n))
        evalues = wr(order)
        evect = vr(:,order)

    end subroutine get_evect

    function gauss_similarity(distance, bandwidth)

        implicit none
        real(8) :: bandwidth
        real(8), dimension(:,:), allocatable, intent(in) :: distance
        real(8), dimension(:,:), allocatable :: gauss_similarity

        allocate(gauss_similarity(size(distance,1),size(distance,2)))
        gauss_similarity = exp( - ( (distance**2) / (2*bandwidth) ) )

    end function gauss_similarity

    subroutine diffusion_map_run(this, distance, bandwidth, time, dimensions)

        implicit none
        integer :: i, j, n
        integer, intent(in), optional :: dimensions
        real(8), dimension(:,:), allocatable :: similarity, markov_transition, evec
        real(8), dimension(:), allocatable :: eval
        real(8), dimension(:,:), allocatable, intent(in) :: distance
        real(8), intent(in) :: time, bandwidth
        class(diffusion_map_type), intent(inout) :: this

        n = size(distance,1)

        ! Using a Gaussian kernel
        similarity = gauss_similarity(distance, bandwidth)

        ! Normalize the diffusion / similarity matrix
        allocate(markov_transition(n,n))
        do i = 1, n
            markov_transition(i,:) = similarity(i,:) / sum(similarity(i,:))
        end do

        call get_evect(markov_transition, evec, eval)

        allocate(this%eval(1:dimensions+1))
        this%eval = eval(1:dimensions+1)
        allocate(this%evec(size(evec,1),1:dimensions+1))
        this%evec = evec(:,1:dimensions+1)

        ! Note that we do not save the first eigenvector to the map since it is trivial (all 1's)
        ! We do keep it in "evect" just as a check
        allocate(this%map(n,dimensions))
        do j = 2, dimensions+1
            this%map(:,j-1) = this%evec(:,j)*this%eval(j)**time
        end do

    end subroutine diffusion_map_run

end module diffusion_map
