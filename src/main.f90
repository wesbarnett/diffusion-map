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

    ! kept separate from diffusion_map module since distance metric is determined depending only the application
    function get_distance(indata)

        implicit none
        real(8), allocatable :: get_distance(:,:)
        real(8), allocatable, intent(in) :: indata(:,:)
        real(8) :: s
        integer :: ndata, ndim, i, j, k
    
        ndim = size(indata,1)
        ndata = size(indata,2)

        ! only appropriate for cartesian data. Real world use we would use a different metric here
        allocate(get_distance(ndata,ndata))

        do i = 1, ndata-1
            do j = i+1, ndata
                s = 0.0d0
                do k = 1, ndim
                    s = s + (indata(k,i)-indata(k,j))**2
                end do
                get_distance(i,j) = dsqrt(s)
                get_distance(j,i) = get_distance(i,j)
            end do
        end do

    end function get_distance

end module subs

program main

    use diffusion_map
    use princ_comp
    use json_module
    use subs

    implicit none
    integer :: i, j, n, u, logbandwidth_l, logbandwidth_u, max_output, dimensions
    real(8), dimension(:,:), allocatable :: distance, similarity, point
    real(8), dimension(:), allocatable :: val
    real(8) :: bandwidth
    logical :: get_bandwidth, found, run_pca, run_dmap
    character (len=256), allocatable :: config_file
    character (len=32) :: arg, n_char
    character (len=:), allocatable :: infile, bandwidth_file, evects_file, evalues_file, diffusionmap_file, pca_evects_file, &
        pca_evalues_file, pca_file
    character (len=1024) :: format_string
    type(json_file) :: config
    real(8) :: time ! diffusion "time", not simulation time
    type(diffusion_map_type) :: dm
    type(princ_comp_type) :: pca

    if (command_argument_count() .ne. 1) then
        write(0,*) "ERROR: First argument should be config file."
        stop
    end if

    call get_command_argument(1, arg)

    config_file = trim(arg)
    call check_file(config_file)
    call config%initialize()
    call config%load_file(filename=config_file)
    call config%print_file()
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
    call config%get("dimensions",dimensions,found)
    if (.not. found) then 
        dimensions = 3
    end if
    ! Default is 4 because we are using 3d data for this example (and first eigenvector is all 1's and is ignored)
    max_output = dimensions + 1
    call config%get("infile",infile,found)
    if (.not. found) then 
        infile = "infile.dat"
    end if
    call config%get("dmap.time",time,found)
    if (.not. found) then 
        time = 0.0
    end if
    call check_file(infile)
    call config%get("dmap.evects",evects_file,found)
    if (.not. found) then 
        evects_file = "evects.dat"
    end if
    call config%get("dmap.evalues",evalues_file,found)
    if (.not. found) then 
        evalues_file = "evalues.dat"
    end if
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
    call config%destroy()

    ! val is the original data's position on the swiss roll
    open(newunit=u, file=trim(infile), status="old")
    read(u,*) n
    allocate(point(dimensions,n))
    allocate(val(n))
    do i = 1, n
        read(u,*) point(:,i), val(i)
    end do
    close(u) 

    if (run_dmap .or. get_bandwidth) then

        distance = get_distance(point)

        if (get_bandwidth) then

            ! Cycle through different values of the bandwidth. See Fig. S1 in https://www.pnas.org/cgi/doi/10.1073/pnas.1003293107
            open(newunit=u, file=trim(bandwidth_file))
            do i = logbandwidth_l, logbandwidth_u
                bandwidth = exp(dble(i))
                similarity = gauss_similarity(distance, bandwidth)
                write(u,"(i0,f12.6)") i, log(sum(similarity))
            end do
            close(u)

        else

            call dm%run(distance, bandwidth, time, dimensions)

            open(newunit=u, file=trim(evalues_file))
            write(u,"(f12.6)") dm%eval
            close(u)

            open(newunit=u, file=trim(evects_file))
            write(n_char,'(i0)') max_output
            format_string = "("//trim(n_char)//"f12.6)"
            write(u,format_string) transpose(dm%evec)
            close(u)

            open(newunit=u, file=trim(diffusionmap_file))
            write(u,"(a)") "# First column is original location on swiss roll"
            write(u,"(a)") "# Next columns are points in diffusion space"
            write(u,"(a)") "# Note that first (trivial) eigenvector not used."
            write(u,"(a)") "# In gnuplot to plot the first dimensions try:"
            write(u,"(a)") "#   plot 'dmap.dat' u 2:3:1 w points palette"
            write(u,"(a,f12.6)") "# time = ", time
            write(u,"(a,f12.6)") "# bandwidth = ", bandwidth
            format_string = "("//trim(n_char)//"f12.6)"
            do i = 1, n
                write(u,"(f12.6)", advance="no") val(i)
                do j = 1, dimensions
                    write(u,"(f12.6)", advance="no") dm%map(i,j)
                end do
                write(u,*)
            end do
            close(u)

        end if

    end if

    if (run_pca) then

        call pca%run(point)

        open(newunit=u, file=trim(pca_evalues_file))
            write(u,"(a)") "# Eigenvalues"
            write(u,"(f12.6)") pca%eval
            write(u,"(a)") "# Contribution"
            write(u,"(f12.6)") pca%contrib
            write(u,"(a)") "# Cumulative contribution"
            write(u,"(f12.6)") pca%cumul_contrib
        close(u)

        open(newunit=u, file=trim(pca_evects_file))
        write(n_char,'(i0)') dimensions
        format_string = "("//trim(n_char)//"f12.6)"
        write(u,format_string) transpose(pca%evec)
        close(u)

        open(newunit=u, file=trim(pca_file))
        ! NOTE: dimensions for pca data are swapped
        do j = 1, n
            write(u,"(f12.6)", advance="no") val(j)
            do i = 1, dimensions
                write(u,"(f12.6)", advance="no") pca%x(i,j)
            end do
            write(u,*)
        end do
        close(u)
        
    end if

end program main
