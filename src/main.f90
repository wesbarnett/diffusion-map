! Diffusion map code example
! James W. Barnett

module subs

contains 

    subroutine check_file(myfile)

        implicit none
        character (len=*), intent(in) :: myfile
        logical ex

        inquire(file=trim(myfile), exist=ex)
        if (ex .eqv. .false.) then
            write(*,"(a)") "ERROR: "//trim(myfile)//" does not exist."
            stop
        end if

    end subroutine

    function get_rmsd(trj)

        use gmxfort_trajectory
        use omp_lib

        implicit none
        real(8), allocatable :: get_rmsd(:,:)
        type(Trajectory), intent(inout) :: trj
        integer :: nframe, natoms, nthreads

        nframe = trj%nframes
        natoms = trj%natoms()
        allocate(get_rmsd(nframe,nframe))

        !$omp parallel
        !$omp master
            nthreads = omp_get_num_threads()
            write(*,"(a,i0,a)") "Using ", nthreads, " OpenMP threads..."
        !$omp end master
        !$omp end parallel
     
        !$omp parallel 

        block

            integer :: i, j, k, l
            real(8) :: s, atom_i(3), atom_j(3)

            !$omp do 
            do i = 1, nframe-1

                write(*,*) i

                do j = i+1, nframe

                    s = 0.0d0
                    do k = 1, natoms

                        atom_i = trj%x(i,k)
                        atom_j = trj%x(j,k)
                        s = s + sum((atom_i - atom_j)**2)

                    end do

                    s = s / dble(natoms)

                    get_rmsd(i,j) = dsqrt(s)
                    get_rmsd(j,i) = get_rmsd(i,j)

                end do

            end do
            !$omp end do

        end block

        !$omp end parallel

    end function get_rmsd

end module subs

program main

    use gmxfort_trajectory
    use gmxfort_utils, only: dihedral_angle
    use diffusion_map
    use princ_comp
    use json_module
    use subs

    implicit none
    integer :: i, j, n, u, logbandwidth_l, logbandwidth_u, dimensions, nangles, ncomp
    real(8), dimension(:,:), allocatable :: distance, indata
    real(8), dimension(:), allocatable :: val, logsumsim
    real(8) :: bandwidth
    logical :: get_bandwidth, found, run_pca, run_dmap
    character (len=256), allocatable :: config_file
    character (len=32) :: arg, n_char
    character (len=:), allocatable :: bandwidth_file, evects_file, evalues_file, diffusionmap_file, pca_evects_file, &
        pca_evalues_file, pca_file, xtcfile, ndxfile, ndxgroup
    character (len=1024) :: format_string
    type(json_file) :: config
    real(8) :: time ! diffusion "time", not simulation time
    type(diffusion_map_type) :: dm
    type(princ_comp_type) :: pca
    type(Trajectory) :: trj

    if (command_argument_count() .ne. 1) then
        write(0,*) "ERROR: First argument should be config file."
        stop
    end if

    call get_command_argument(1, arg)

    config_file = trim(arg)
    write(*,*) "Reading in configuration from "//trim(config_file)//"..."
    call check_file(config_file)
    call config%initialize()
    call config%load_file(filename=config_file)
    call config%print_file()

    ! Bandwidth calc input
    call config%get('bandwidth_calc.get',get_bandwidth,found)
    if (.not. found) then 
        get_bandwidth=.false.
    end if
    call config%get("bandwidth_calc.log_l",logbandwidth_l,found)
    if (.not. found) then 
        logbandwidth_l = -10
    end if
    call config%get("bandwidth_calc.log_u",logbandwidth_u,found)
    if (.not. found) then 
        logbandwidth_l = 10
    end if
    call config%get("bandwidth_calc.file",bandwidth_file,found)
    if (.not. found) then 
        bandwidth_file = "eps.dat"
    end if

    ! DIffusion map input
    call config%get('dmap.run',run_dmap,found)
    if (.not. found) then 
        run_dmap = .false.
    end if
    call config%get("dmap.outfile",diffusionmap_file,found)
    if (.not. found) then 
        diffusionmap_file = "dmap.dat"
    end if
    call config%get("dmap.bandwidth",bandwidth,found)
    if (.not. found) then 
        bandwidth = 1.0
    end if
    call config%get("dmap.time",time,found)
    if (.not. found) then 
        time = 0.0
    end if
    call config%get("dmap.evects",evects_file,found)
    if (.not. found) then 
        evects_file = "evects.dat"
    end if
    call config%get("dmap.evalues",evalues_file,found)
    if (.not. found) then 
        evalues_file = "evalues.dat"
    end if

    ! PCA input
    call config%get('pca.run',run_pca,found)
    if (.not. found) then 
        run_pca = .false.
    end if
    call config%get("pca.evects",pca_evects_file,found)
    if (.not. found) then 
        pca_evects_file = "evects.dat"
    end if
    call config%get("pca.evalues",pca_evalues_file,found)
    if (.not. found) then 
        pca_evalues_file = "evalues.dat"
    end if
    call config%get("pca.outfile",pca_file,found)
    if (.not. found) then 
        pca_file = "pca.dat"
    end if

    ! Simulation input
    call config%get('sim.xtc',xtcfile,found)
    if (.not. found) then 
        xtcfile = "traj.xtc"
    end if
    call config%get('sim.ndx',ndxfile,found)
    if (.not. found) then 
        ndxfile = "index.ndx"
    end if
    call config%get('sim.group',ndxgroup,found)
    if (.not. found) then 
        ndxgroup = "C"
    end if

    call config%destroy()

! FIXME: use index file
    call trj%read(xtcfile)
!   call trj%read(xtcfile, ndxfile, ndxgroup)
    dimensions = trj%nframes
    nangles = trj%natoms() - 3
    ncomp = nangles*2

    if (run_dmap .or. get_bandwidth) then

        write(*,*) "Calculating RMSD..."
        distance = get_rmsd(trj)

        if (get_bandwidth) then

            write(*,*) "Performing diffusion map bandwidth analysis..."

            ! Cycle through different values of the bandwidth. See Fig. S1 in https://www.pnas.org/cgi/doi/10.1073/pnas.1003293107

            allocate(logsumsim(logbandwidth_l:logbandwidth_u))

            !$omp parallel 

            block

                real(8) :: bandwidth
                real(8), dimension(:,:), allocatable :: similarity
                integer :: i

                !$omp do
                do i = logbandwidth_l, logbandwidth_u
                    write(*,*) i
                    bandwidth = exp(dble(i))
                    similarity = gauss_similarity(distance, bandwidth)
                    logsumsim(i) = log(sum(similarity))
                end do
                !$omp end do

            end block

            !$omp end parallel 

            write(*,*) "Writing diffusion map bandwidth analysis to "//trim(bandwidth_file)//"..."
            open(newunit=u, file=trim(bandwidth_file))
            do i = logbandwidth_l, logbandwidth_u
                write(u,"(i0,f12.6)") i, logsumsim(i)
            end do
            close(u)

        else

            write(*,*) "Performing diffusion map calculation..."
            write(*,*) "bandwidth = ", bandwidth
            write(*,*) "time = ", time
            write(*,*) "dimensions = ", dimensions
            call dm%run(distance, bandwidth, time, trj%nframes)

            write(*,*) "Writing diffusion map eigenvalues to "//trim(evalues_file)//"..."
            open(newunit=u, file=trim(evalues_file))
            write(u,"(f12.6)") dm%eval
            close(u)

            write(*,*) "Writing diffusion map eigenectors to "//trim(evects_file)//"..."
            open(newunit=u, file=trim(evects_file))
            write(n_char,'(i0)') dimensions+1
            format_string = "("//trim(n_char)//"f12.6)"
            write(u,format_string) transpose(dm%evec)
            close(u)

            write(*,*) "Writing diffusion map to "//trim(diffusionmap_file)//"..."
            open(newunit=u, file=trim(diffusionmap_file))
            write(u,"(a)") "# First column is original location on swiss roll"
            write(u,"(a)") "# Next columns are points in diffusion space"
            write(u,"(a)") "# Note that first (trivial) eigenvector not used."
            write(u,"(a)") "# In gnuplot to plot the first dimensions try:"
            write(u,"(a)") "#   plot 'dmap.dat' u 2:3:1 w points palette"
            write(u,"(a,f12.6)") "# time = ", time
            write(u,"(a,f12.6)") "# bandwidth = ", bandwidth
            write(n_char,'(i0)') dimensions
            format_string = "("//trim(n_char)//"f12.6)"
            do i = 1, dimensions
                do j = 1, dimensions
                    write(u,"(f12.6)", advance="no") dm%map(i,j)
                end do
                write(u,*)
            end do
            close(u)

        end if

    end if

    if (run_pca) then

        write(*,*) "Performing PCA..."
        ! Calculate cosine and sine of dihedral angles for dihedral principal components analysis
        allocate(indata(ncomp, trj%nframes))
        !$omp parallel 

        block

            integer :: i, j, k
            real(8) :: a(3), b(3), c(3), d(3), box(3,3)
            real(8), allocatable :: ang(:)

            allocate(ang(nangles))

            !$omp do 
            do i = 1, trj%nframes

                write(*,*) i

                do j = 1, nangles
                    
                    a = trj%x(i, j)
                    b = trj%x(i, j+1)
                    c = trj%x(i, j+2)
                    d = trj%x(i, j+3)
                    box = trj%box(i)
                    ang(j) = dihedral_angle(a, b, c, d, box)

                end do

                k = 1
                do j = 1, ncomp, 2
                    indata(j,i) = dcos(ang(k));
                    indata(j+1,i) = dsin(ang(k));
                    K = K + 1
                end do

            end do
            !$omp end do

            deallocate(ang)

        end block

        !$omp end parallel

        call pca%run(indata)

        write(*,*) "Writing PCA eigenvalues to "//trim(pca_evalues_file)//"..."
        open(newunit=u, file=trim(pca_evalues_file))
            write(u,"(a)") "# Eigenvalues"
            write(u,"(f12.6)") pca%eval
            write(u,"(a)") "# Contribution"
            write(u,"(f12.6)") pca%contrib
            write(u,"(a)") "# Cumulative contribution"
            write(u,"(f12.6)") pca%cumul_contrib
        close(u)

        write(*,*) "Writing PCA eigenvectors to "//trim(pca_evects_file)//"..."
        open(newunit=u, file=trim(pca_evects_file))
        write(n_char,'(i0)') ncomp
        format_string = "("//trim(n_char)//"f12.6)"
        write(u,format_string) transpose(pca%evec)
        close(u)

        write(*,*) "Writing PCA data to "//trim(pca_file)//"..."
        open(newunit=u, file=trim(pca_file))
        do j = 1, trj%nframes
            do i = 1, ncomp
                write(u,"(f12.6)", advance="no") pca%x(i,j)
            end do
            write(u,*)
        end do
        close(u)
        
    end if

end program main
