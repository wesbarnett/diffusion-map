
module subs
    
    implicit none
    integer, parameter :: X = 1, Y = 2, Z = 3, ndim = 3

contains

    subroutine get_evect_sym(matrix, evect, evalues)

        implicit none
        character :: jobz = 'V', uplo = 'U'
        real(8), allocatable ::  work(:), a(:,:)
        real(8), allocatable, intent(inout) ::  evalues(:), evect(:,:), matrix(:,:)
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

        a = matrix
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


end module subs

program main

    use gmxfort_trajectory
    use json_module
    use subs

    implicit none
    character (len=256), allocatable :: config_file
    character (len=32) :: arg
    real(8) :: rcm(3), atom(3), atomj(3), atomk(3)
    real(8), allocatable :: S(:,:), S_evect(:,:), S_eval(:)
    character (len=:), allocatable :: outfile, xtcfile
    type(json_file) :: config
    type(Trajectory) :: trj
    logical :: found
    integer :: i, j, k, l, natoms, nframes, u, m

    if (command_argument_count() .ne. 1) then
        write(0,*) "ERROR: First argument should be config file."
        stop
    end if

    call get_command_argument(1, arg)

    config_file = trim(arg)
    write(*,*) "Reading in configuration from "//trim(config_file)//"..."

    call config%initialize()
    call config%load_file(filename=config_file)
    call config%print_file()

    ! Simulation input
    call config%get('xtc',xtcfile,found)
    if (.not. found) then 
        xtcfile = "traj.xtc"
    end if

    call config%get('outfile',outfile,found)
    if (.not. found) then 
        xtcfile = "out.dat"
    end if

    call config%destroy()

    call trj%read(xtcfile)
    nframes = trj%nframes
    natoms = trj%natoms()

    allocate(S(3,3))

    open(newunit=u, file=outfile)
    do i = 1, nframes

        ! Get center of particles
        ! Gyration tensor does not use mass (inertia tensor does)
        rcm = 0.0d0
        do j = 1, natoms
            rcm = rcm + trj%x(i,j)        
        end do
        rcm = rcm / dble(natoms)

        ! Calculate gyration tensor
        S = 0.0d0
        do j = 1, natoms

            atom = trj%x(i,j) - rcm

            do k = 1, ndim-1

                do l = k+1, ndim

                    S(k,l) = S(k,l) + atom(k)*atom(l)
                    S(l,k) = S(k,l)

                end do

            end do

        end do
        S = S / dble(natoms)

        ! Diagonalize 
        call get_evect_sym(S, S_evect, S_eval)       
        write(u, "(3f12.6)") S_eval
        deallocate(S_evect, S_eval)

    end do

    close(u)

end program main
