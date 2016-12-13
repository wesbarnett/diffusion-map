module subs

    implicit none
    real(8), parameter :: pi = 2.0d0*acos(0.0d0)

contains

    subroutine init_seed()

        implicit none
        integer :: n, u
        integer, allocatable :: seed(:)

        call random_seed(size=n)

        allocate(seed(n))
        open(newunit=u, file="/dev/urandom", access="stream", form="unformatted", &
            action="read", status="old")
        read(u) seed(:)
        close(u)

        call random_seed(put=seed)

    end subroutine init_seed

    function gen_sphere_point(center, r)

        implicit none
        real(8) :: xi(2),  zeta(2),  zeta2, gen_sphere_point(3)
        real(8), intent(in) :: center(3), r

        zeta2 = 100.0d0

        do while (zeta2 .gt. 1.0d0)
            call random_number(xi)
            zeta = 1.0d0 - 2.0d0*xi
            zeta2 = sum(zeta**2) 
        end do

        gen_sphere_point = [ 2.0d0 * zeta(1) * sqrt(1.0d0 - zeta2), &
                             2.0d0 * zeta(2) * sqrt(1.0d0 - zeta2), &
                             1.0d0 - 2.0d0*zeta2 ]
        gen_sphere_point = gen_sphere_point * r + center

    end function gen_sphere_point

end module subs

program main

    use subs

    implicit none
    integer :: i, n, u
    real(8), allocatable :: z(:), data3d(:,:)
    real(8) :: Rx(3,3), Ry(3,3), Rz(3,3), angles(3), mindim, point(3), a, b, g, center(3), r, rotate(3,3), rand(3)
    logical :: do_rotate
    character (len=:), allocatable :: outfile

    ! TODO: Read in from config file
    do_rotate = .true.
    outfile = "puncsphere.dat"
    n = 2000
    center = [ 0.0, 0.0, 0.0 ]
    r = 0.5

    call init_seed()

    allocate(data3d(3,n))
    allocate(z(n))
    do i = 1, n
        do while (.true.)
            point = gen_sphere_point(center, r) 
            call random_number(rand)
            if ((rand(1)-0.5d0) .le. point(3) .and. point(3) .lt. (center(3)+0.75*r)) then 
                exit
            end if
        end do
        data3d(:,i) = point
        z(i) = point(3)
    end do

    ! Add some noise
    mindim = min(maxval(data3d(3,:)) - minval(data3d(3,:)), maxval(data3d(2,:)) - minval(data3d(2,:)), &
        maxval(data3d(1,:)) - minval(data3d(1,:)))
    do i = 1, n
        call random_number(rand)
        data3d(:,i) = data3d(:,i) + 0.02*rand*sqrt(mindim)
    end do

    ! Rotate all the points 
    if (do_rotate) then
        call random_number(angles)

        a = angles(1) * pi
        b = angles(2) * pi
        g = angles(3) * pi

        Rx = reshape((/ 1.0d0,  0.0d0,   0.0d0,  0.0d0, cos(a), -sin(a),   0.0d0, sin(a), cos(a) /), shape(Rx))
        Ry = reshape((/cos(b),  0.0d0,  sin(b),  0.0d0,  1.0d0,   0.0d0, -sin(b),  0.0d0, cos(b) /), shape(Ry))
        Rz = reshape((/cos(g), -sin(g),  0.0d0, sin(g), cos(g),   0.0d0,   0.0d0,  0.0d0,  1.0d0 /), shape(Rz))

        rotate = matmul(Rz,  matmul(Ry, Rx))
        data3d = matmul(rotate, data3d)
    end if

    open(newunit=u, file=trim(outfile))
    write(u, "(i0)") n
    do i = 1, n
        write(u, "(4f24.3)") data3d(:,i), z(i)
    end do
    close(u)

end program main
