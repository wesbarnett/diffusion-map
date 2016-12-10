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

end module subs

program main

	use subs

    implicit none
    integer :: i, n, u
    real(8), allocatable :: data2d(:,:), data3d(:,:)
	real(8) :: rand(2), width, height, hole_x1, hole_x2, hole_y1, hole_y2, x, y, rotate(3,3), a, b, g, r
    real(8) :: Rx(3,3), Ry(3,3), Rz(3,3), angles(3), mindim

    call random_number(angles)

    a = angles(1) * pi
    b = angles(2) * pi
    g = angles(3) * pi

    Rx = reshape((/ 1.0d0,  0.0d0,   0.0d0,  0.0d0, cos(a), -sin(a),   0.0d0, sin(a), cos(a) /), shape(Rx))
    Ry = reshape((/cos(b),  0.0d0,  sin(b),  0.0d0,  1.0d0,   0.0d0, -sin(b),  0.0d0, cos(b) /), shape(Ry))
    Rz = reshape((/cos(g), -sin(g),  0.0d0, sin(g), cos(g),   0.0d0,   0.0d0,  0.0d0,  1.0d0 /), shape(Rz))

    rotate = matmul(Rz,  matmul(Ry, Rx))

	call init_seed()

    n = 1000
    height = 20.0

    ! Generate a plane of numbers randomly
    allocate(data2d(2,n))
    allocate(data3d(3,n))
    do i = 1, n
		call random_number(rand)
        !r = dble(i)-1.0d0
        !x = pi*(1.5d0+3.0d0*r)
        x = rand(1)*height
        y = rand(2)*height
		data2d(1,i) = x
		data2d(2,i) = y
		data3d(:,i) = [ x*dcos(x), x*dsin(x), y ]
    end do

    ! Add some noise
    mindim = min(maxval(data3d(2,:)) - minval(data3d(2,:)), maxval(data3d(1,:)) - minval(data3d(1,:)))
    do i = 1, n
        call random_number(rand)
        data3d(1:2,i) = data3d(1:2,i) + 0.02*rand*sqrt(mindim)
    end do

    data3d = matmul(rotate, data3d)

	open(newunit=u, file="swissroll.dat")
    write(U, "(i0)") n
    do i = 1, n
        write(U, "(4f24.3)") data3d(:,i), data2d(1,i)
    end do
	close(U)

end program main
