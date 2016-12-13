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
    integer :: i, n, u, nclust, j, k
    real(8), allocatable :: z(:), data3d(:,:)
    real(8) :: Rx(3,3), Ry(3,3), Rz(3,3), angles(3), mindim, point(3), a, b, g, center(3), r, rotate(3,3), rand(3), rand2(3), &
        rand3(3)
    logical :: do_rotate
    character (len=:), allocatable :: outfile

    ! TODO: Read in from config file
    do_rotate = .true.
    outfile = "clusters.dat"
    n = 500
    nclust = 3

    call init_seed()

    allocate(data3d(3,n*nclust))
    allocate(z(n*nclust))
    k = 1
    do i = 1, nclust
        call random_number(rand)
        call random_number(rand3)
        rand = rand*10.0
        do j = 1, n
            call random_number(rand2)
            point =  gen_sphere_point(rand, rand2(1)*dble(i))
            data3d(:,k) = point
            z(k) = i
            k = k + 1
        end do
    end do

    ! Add some noise
    mindim = min(maxval(data3d(3,:)) - minval(data3d(3,:)), maxval(data3d(2,:)) - minval(data3d(2,:)), &
        maxval(data3d(1,:)) - minval(data3d(1,:)))
    do i = 1, n
        call random_number(rand)
        data3d(:,i) = data3d(:,i) + 0.10*rand*sqrt(mindim)
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
    write(u, "(i0)") n*nclust
    do i = 1, n*nclust
        write(u, "(4f24.3)") data3d(:,i), z(i)
    end do
    close(u)

end program main
