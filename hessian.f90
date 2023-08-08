! Name: density_histogram.f90
program density_histogram
    use kdtree2_module
    use kdtree2_precision_module
    use, intrinsic :: iso_fortran_env, only: iostat_end
    implicit none
    INCLUDE "omp_lib.h"
    type(kdtree2), pointer :: tree
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: i, j, k, l, nx, ny, nz, count, ios, hist_res, nn2
    real(dp), allocatable :: x(:,:), mass(:), wgt(:),gridxyz(:)
    real(dp) :: x_min, x_max, y_min, y_max, z_min, z_max,sigma,dlogsfil,dlogswall,dlogsc
    real(dp) :: dM1_max,dM2_max,dM3_max,thresh_c,thresh_f,thresh_w
    real(dp) :: bin_size, field(3,3,3), dM1, dM2, dM3,max_sfil,min_sfil,max_swall,min_swall,max_sc,min_sc
    real(dp), allocatable ::  density(:,:,:), smoothed_density(:,:,:)
    real(dp), allocatable :: eigenvalue1(:,:,:), eigenvalue2(:,:,:), eigenvalue3(:,:,:)
    real(dp), allocatable :: hessian(:,:), eigenvalues(:), sfil(:,:,:,:),swall(:,:,:,:),scluster(:,:,:,:)
    integer, allocatable :: structure(:,:,:)
    character(len=256) :: filename, out
    character(len=32) :: hist_res_str
    real(dp), dimension(3*3) :: work
    integer :: info, threads, nsig, nsmooth
    logical :: write_eig, write_density, calc_density
    type(kdtree2_result),allocatable :: resultsb(:)


!$OMP PARALLEL
!$ threads=OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

print*,'Code running with ',threads,' OMP threads'

    write_eig = .false.
    write_density = .true.
    nsmooth = 8
    calc_density = .false.

    ! Read filename, histogram resolution, calculate density option from command-line arguments
    if (command_argument_count() == 0) then
        print*, 'Usage: hessian <filename> <histogram resolution> <calculate density>'
        print*, '       <filename> is the name of the file containing the data'
        print*, '       <histogram resolution> is the size of the bins in the histogram'
        print*, '       <calculate density> is a boolean flag to calculate the density field'
        print*, 'e.g. hessian input_data.txt 1.1 .true.'
        stop
    else
        call get_command_argument(1, filename)
        call get_command_argument(2, hist_res_str)
        read(hist_res_str, *) bin_size
        call get_command_argument(3, hist_res_str)
        read(hist_res_str, *) calc_density
    end if


 ! Read in the data from a text file
    !filename = 'input_data.txt'
    !filename = 'tmptmptmp'
    !filename = 'test.gal.txt.den'
    !filename = 'tmp.den'

    print*, 'reading data from ', filename
    open(10, file=filename, status='old', action='read')
    count = 0
    do
        read(10, *, end=100, iostat=ios) x_min 
        if (ios == iostat_end) exit
        count = count + 1
    end do
100 continue
    close(10)

    print*, 'there are ', count, ' particles in the file'

    allocate(x(3, count), mass(count), wgt(count))

    call read_data()
  
    print*, 'data read in'



    ! Determine the bounding box of the input data
    x_min = minval(x(1,:)); x_max = maxval(x(1,:))
    y_min = minval(x(2,:)); y_max = maxval(x(2,:))
    z_min = minval(x(3,:)); z_max = maxval(x(3,:))

    print*, 'bounding box determined'
    print*, 'x_min = ', x_min, ' x_max = ', x_max
    print*, 'y_min = ', y_min, ' y_max = ', y_max
    print*, 'z_min = ', z_min, ' z_max = ', z_max

    ! determine the histogram resolution
    nx = floor((x_max - x_min) / bin_size) + 1
    ny = floor((y_max - y_min) / bin_size) + 1
    nz = floor((z_max - z_min) / bin_size) + 1

    print*, 'histogram resolution determined'
    print*, 'nx = ', nx, ' ny = ', ny, ' nz = ', nz

    ! Allocate memory for input data arrays
    allocate(density(nx, ny, nz))
    allocate(smoothed_density(nx, ny, nz))
    allocate(eigenvalue1(nx, ny, nz))
    allocate(eigenvalue2(nx, ny, nz))
    allocate(eigenvalue3(nx, ny, nz))
    allocate(hessian(3,3))
    allocate(eigenvalues(3))
    allocate(sfil(nx, ny, nz, nsmooth))
    allocate(swall(nx, ny, nz, nsmooth))
    allocate(scluster(nx, ny, nz,nsmooth))
    allocate(gridxyz(3))
    allocate(resultsb(count))
    allocate(structure(nx, ny, nz))

    print*, 'memory allocated'
    
    if (calc_density) then
        call calculate_density()
        print*, 'density calculated'
    else

    ! log the mass field (add 0.1 to avoid log(0))
    !mass = log10(mass+0.1)
    mass = log10(mass)

    ! create a kd-tree for the data
    tree => kdtree2_create(x,sort=.true.,rearrange=.true.)     ! this is how you create a tree.
    print*,'created tree'

    ! create a grid of points to interpolate the density field onto (x,y,z)
    density = 0.0_dp

    nn2=10 ! number of nearest neighbours used for density estimation
    !$OMP PARALLEL DO PRIVATE(i,j,k,gridxyz,wgt,resultsb) SHARED(density,tree)
    do i = 1,nx
        do j = 1,ny
            do k = 1,nz
                !l=l+1
                gridxyz(:) = [bin_size*(i-1)+x_min,bin_size*(j-1)+y_min,bin_size*(k-1)+z_min]       
                ! call kd-tree to find the nearest neighbours to each grid point
                !call kdtree2_r_nearest(tp=tree,qv=gridxyz(:,l),r2=(bin_size*5)**2,nfound=nn2,nalloc=(count),results=resultsb)  
                call kdtree2_n_nearest(tp=tree,qv=gridxyz,nn=nn2,results=resultsb(1:nn2))
                !if(nn2==0) then
                !    density(i,j,k) = 0.0_dp
                !    cycle
                !endif
                ! calculate the rbf weights
                call rbf_weight(3,nn2,x(:,resultsb(1:nn2)%idx),1.d0,phi5,mass(resultsb(1:nn2)%idx),wgt)
                ! interpolate the density field onto the grid using the rbf weights
                call rbf_interp_nd(3,nn2,x(:,resultsb(1:nn2)%idx),1.d0,phi5,wgt,1,gridxyz, density(i,j,k))

                !print*,i,j,k,density(i,j,k)
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    ! interpolate the density field onto the grid using the rbf weights
    !call rbf_interp_nd(3,count,x,1.d0,phi5,wgt,nx*ny*nz,gridxyz, tmp)

    print*, 'density grid created'
    print*, maxval(density)

endif

    ! normalise the log density field
    density = log10(10**density/sum(10**density))

    do l=1,nsmooth
        sigma=bin_size * sqrt(2.d0)**(l-1)
        ! Smooth the density field with a Gaussian function
        call gaussian_smooth(density, smoothed_density, nx, ny, nz, sigma)
        !call fft_gaussian_smooth(density, smoothed_density, nx, ny, nz, sigma)

        ! normalise the smoothed log density field
        smoothed_density = log10(10**smoothed_density/sum(10**smoothed_density))

        print*, 'density field smoothed with sigma = ', sigma
        print*, sum(10**density), sum(10**smoothed_density)

    ! calculate the hessian eigenvalues at each point in the smoothed log-density field
    ! calculate and save the filament and wall 'signatures' at each point and for each smoothing sigma
    ! parallelise over the grid points
!$OMP PARALLEL DO PRIVATE(i,j,k,hessian,eigenvalues,work,info, field) SHARED(smoothed_density,eigenvalue1,eigenvalue2,eigenvalue3)
    do i = 2, nx-1
        do j = 2, ny-1
            do k = 2, nz-1
                field = smoothed_density(i-1:i+1, j-1:j+1, k-1:k+1)
                call hessian_at_point(field, bin_size, hessian)
                hessian=hessian*sigma**2
                !write(20, '(3i5,1x,e23.15)') i, j, k , hessian(3,1)
                call dsyev('N', 'U', 3, hessian, 3, eigenvalues, work, 3*3, info)
                eigenvalue1(i, j, k) = eigenvalues(1)
                eigenvalue2(i, j, k) = eigenvalues(2)
                eigenvalue3(i, j, k) = eigenvalues(3)
                !print *, i, j, k, eigenvalues(1), eigenvalues(2), eigenvalues(3)

                scluster(i, j, k, l) = (eigenvalues(3)**2/abs(eigenvalues(1)))*heaviside(-eigenvalues(3))*heaviside(-eigenvalues(2))*heaviside(-eigenvalues(1))
                sfil(i,j,k,l) = (eigenvalues(2)**2/abs(eigenvalues(1)))*(1-abs(eigenvalues(3)/eigenvalues(1)))*heaviside(1-abs(eigenvalues(3)/eigenvalues(1)))*heaviside(-eigenvalues(2))*heaviside(-eigenvalues(1))
                swall(i,j,k,l) = abs(eigenvalues(1))*(1-abs(eigenvalues(3)/eigenvalues(1)))*(1-abs(eigenvalues(2)/eigenvalues(1)))*heaviside(1-abs(eigenvalues(3)/eigenvalues(1)))*heaviside(1-abs(eigenvalues(2)/eigenvalues(1)))*heaviside(-eigenvalues(1))

                !scluster(i, j, k, l) = eigenvalues(3)**2/abs(eigenvalues(1))
                !sfil(i,j,k,l) = (eigenvalues(2)**2/abs(eigenvalues(1)))*(1.d0-abs(eigenvalues(3)/eigenvalues(1)))
                !swall(i,j,k,l) = abs(eigenvalues(1))*(1.d0-abs(eigenvalues(3)/eigenvalues(1)))*(1.d0-abs(eigenvalues(2)/eigenvalues(1)))

                !if(eigenvalues(1) .gt. 0.0_dp) then
                !    swall(i,j,k,l) = 0.0_dp
                !    if(eigenvalues(2) .gt. 0.0_dp) then
                !        sfil(i,j,k,l) = 0.0_dp
                !        if(eigenvalues(3) .gt. 0.0_dp) then
                !            scluster(i,j,k,l) = 0.0_dp
                !        endif
                !    endif
                !endif

            end do
        end do
    end do
!$OMP END PARALLEL DO

    ! write the smoothed log-density field to file (if requested) smoothed_density_N.txt
    if(write_density) then
        write(out, '(a,i0,a)') 'smoothed_density_', l, '.txt'
        out = trim(adjustl(out)) 
        open(20, file=out, status='replace', action='write')
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    write(20, '(3i5,1x,e23.15)') i, j, k , smoothed_density(i, j, k)
                end do
            end do
        enddo
        close(20)
    endif

end do

print*, 'eigenvalues calculated for each smoothing scale'

    ! find the maximum of sfil and swall at each cell
    do i = 2, nx-1
        do j = 2, ny-1
            do k = 2, nz-1
                scluster(i, j, k,1) = maxval(scluster(i, j, k, :))
                sfil(i, j, k,1) = maxval(sfil(i, j, k, :))
                swall(i, j, k,1) = maxval(swall(i, j, k, :))
            end do
        end do
    end do

    !sfil=sfil+0.01
    !swall=swall+0.01
    max_sc = maxval(scluster(:, :, :,1))
    max_sfil = maxval(sfil(:, :, :,1))
    max_swall = maxval(swall(:, :, :,1))
    min_sc = 0.001 !minval(scluster(:, :, :,1))+0.01
    min_sfil = 0.001 !minval(sfil(:, :, :,1))+0.01
    min_swall = 0.001 !minval(swall(:, :, :,1))+0.01
    !min_sc    = minval(scluster(:, :, :,1))
    !min_sfil  = minval(sfil(:, :, :,1))
    !min_swall = minval(swall(:, :, :,1))

    print*, 'max_sc = ', max_sc
    print*, 'max_sfil = ', max_sfil
    print*, 'max_swall = ', max_swall
    print*, 'min_sc = ', min_sc
    print*, 'min_sfil = ', min_sfil
    print*, 'min_swall = ', min_swall

    nsig = 30
    dlogsc = (log10(max_sc) - log10(min_sc)) / nsig
    dlogsfil = (log10(max_sfil) - log10(min_sfil)) / nsig
    dlogswall = (log10(max_swall) - log10(min_swall)) / nsig

    do i = 1,nsig
        dM1 = sum(10.d0**smoothed_density, mask = (scluster(:,:,:,1) >= 10.0_dp**(log10(min_sc) + dlogsc * (i-1))) .and. scluster(:,:,:,1) < 10.0_dp**(log10(min_sc) + dlogsc * i))
        if(dM1>dM1_max) then
            thresh_c=10.0_dp**(log10(min_sc) + dlogsc * (i-1))
            dM1_max=dM1
        endif
        dM2 = sum(10.d0**smoothed_density, mask = (sfil(:,:,:,1) >= 10.0_dp**(log10(min_sfil) + dlogsfil * (i-1))) .and. sfil(:,:,:,1) < 10.0_dp**(log10(min_sfil) + dlogsfil * i))
        if(dM2>dM2_max) then
            thresh_f=10.0_dp**(log10(min_sfil) + dlogsfil * (i-1))
            dM2_max=dM2
        endif
        dM3 = sum(10.d0**smoothed_density, mask = (swall(:,:,:,1) >= 10.0_dp**(log10(min_swall) + dlogswall * (i-1))) .and. swall(:,:,:,1) < 10.0_dp**(log10(min_swall) + dlogswall * i))
        if(dM3>dM3_max) then
            thresh_w=10.0_dp**(log10(min_swall) + dlogswall * (i-1))
            dM3_max=dM3
        endif
        !print*, 10.0_dp**(log10(min_sfil) + dlogsfil * (i-1)), dM2, dM3
        print*, 10.0_dp**(log10(min_sc) + dlogsc * (i-1)), dM1,10.0_dp**(log10(min_sfil) + dlogsfil * (i-1)), dM2, 10.0_dp**(log10(min_swall) + dlogswall * (i-1)), dM3
    enddo

    do i = 1,nx
        do j = 1, ny
            do k = 1, nz
                if(scluster(i, j, k,1)>thresh_c) then
                    structure(i,j,k)=3
                elseif(sfil(i, j, k,1)>thresh_f) then
                    structure(i,j,k)=2
                elseif(swall(i, j, k,1)>thresh_w) then
                    structure(i,j,k)=1
                else
                    structure(i,j,k)=0
                endif
            enddo
        enddo
    enddo

    ! Output eigenvalues
    if(write_eig) then
    open(20, file='eigenvalues.txt', status='replace', action='write')
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                write(20, '(3i5,1x,3e23.15)') i, j, k, eigenvalue1(i, j, k), eigenvalue2(i, j, k), eigenvalue3(i, j, k)
            end do
        end do
    end do
    close(20)
    endif   

    ! Output the maximum of sfil
    open(20, file='sfil.txt', status='replace', action='write')
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                write(20, '(3i5,1x,10e23.15)') i, j, k, sfil(i, j, k,1)
            end do
        end do
    end do
    close(20)

    ! Output the maximum of swall
    open(20, file='swall.txt', status='replace', action='write')
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                write(20, '(3i5,1x,e23.15)') i, j, k, swall(i, j, k,1)
            end do
        end do
    end do
    close(20)

    ! Output the maximum of scluster
    open(20, file='scluster.txt', status='replace', action='write')
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                write(20, '(3i5,1x,e23.15)') i, j, k, scluster(i, j, k,1)
            end do
        end do
    end do

    ! Output the density field
    if(write_density) then
        open(20, file='density.txt', status='replace', action='write')
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    write(20, '(3i5,1x,e23.15)') i, j, k, density(i, j, k)
                end do
            end do
        end do
        close(20)
    endif

    ! Output final structure file
   !! if(write_density) then
        open(20, file='structure.txt', status='replace', action='write')
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    write(20, '(4i5)') i, j, k, structure(i, j, k)
                end do
            end do
        end do
        close(20)
    !endif


    ! output 1 component of the hessian
    !open(20, file='hessian.txt', status='replace', action='write')
    !do i = 1, nx
    !    do j = 1, ny
    !        do k = 1, nz
    !            write(20, '(3i5,1x,e23.15)') i, j, k, hessian(1,1)
    !        end do
    !    end do
    !end do
    !close(20)

contains


    subroutine calculate_density()
        ! bin the points into a grid (mass)

        integer :: i, l, m, n
        density=0.0_dp

        do i = 1, count
            ! Find the grid cell that the point is in
            l = int((x(1,i)-x_min)/bin_size) + 1
            m = int((x(2,i)-y_min)/bin_size) + 1
            n = int((x(3,i)-z_min)/bin_size) + 1

            ! Add the mass to the grid cell
            density(l,m,n) = density(l,m,n) + 1.0_dp

        end do
        density = density + 1e-1_dp
        density = density / bin_size**3
        end subroutine calculate_density

    subroutine read_data()

        ! Re-read the data and store it in the allocated arrays
        open(10, file=filename, status='old', action='read')
        if(calc_density) then
            do i = 1, count
                read(10, *) x(1,i), x(2,i), x(3,i)
            end do

         else

            do i = 1, count
                read(10, *) x(1,i), x(2,i), x(3,i), mass(i)
            end do
        end if
        close(10)

    end subroutine read_data

    subroutine gaussian_smooth(input, output, nx, ny, nz, sigma)
        real(dp), intent(in) :: input(nx,ny,nz), sigma
        real(dp), intent(out) :: output(nx,ny,nz)
        integer, intent(in) :: nx, ny, nz
        integer :: i, j, k, ii, jj, kk, sigma3
        real(dp) :: total_weight, gaussian_exp, dx, dy, dz, weight
        
        sigma3 = int(3 * sigma/bin_size)

        output = 0.0_dp
        ! smooth field
        ! add parallelization with dynamic schedule
        !$OMP PARALLEL DO PRIVATE( i,j,k, ii, jj, kk, weight, total_weight,gaussian_exp,dx,dy,dz) SHARED(input, output, sigma, sigma3, nx, ny, nz)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    total_weight = 0.0_dp
                    do ii = max(1, i-sigma3), min(nx, i+sigma3)
                        do jj = max(1, j-sigma3), min(ny, j+sigma3)
                            do kk = max(1, k-sigma3), min(nz, k+sigma3)
                                dx = real(ii - i, dp) * bin_size
                                dy = real(jj - j, dp) * bin_size
                                dz = real(kk - k, dp) * bin_size
                                gaussian_exp = - (dx*dx + dy*dy + dz*dz) / (2.0_dp * sigma * sigma)
                                weight = exp(gaussian_exp)
                                output(i, j, k) = output(i, j, k) + input(ii, jj, kk) * weight
                                total_weight = total_weight + weight
                            end do
                        end do
                    end do
                    output(i, j, k) = output(i, j, k) / total_weight
                end do
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine gaussian_smooth
    

    subroutine gaussian_smooth2(input, output, nx, ny, nz, sigma)
        real(dp), intent(in) :: input(nx,ny,nz), sigma
        real(dp), intent(out) :: output(nx,ny,nz)
        real(dp), allocatable :: tmp(:,:,:), kernal(:,:,:)
        integer, intent(in) :: nx, ny, nz
        integer :: i, j, k, ii, jj, kk, sigma3
        real(dp) :: total_weight, gaussian_exp, dx, dy, dz, weight
        
        sigma3 = int(3 * sigma/bin_size)

        allocate(tmp(nx+2*sigma3,ny+2*sigma3,nz+2*sigma3))
        allocate(kernal(2*sigma3+1,2*sigma3+1,2*sigma3+1))

        tmp = 0.0_dp
        tmp(sigma3+1:nx+sigma3,sigma3+1:ny+sigma3,sigma3+1:nz+sigma3) = input + 1d-10

        !precompute the gaussian kernal
        do i = 1, 2*sigma3+1
            do j = 1, 2*sigma3+1
                do k = 1, 2*sigma3+1
                    dx = real(i - sigma3 - 1, dp) * bin_size
                    dy = real(j - sigma3 - 1, dp) * bin_size
                    dz = real(k - sigma3 - 1, dp) * bin_size
                    gaussian_exp = - (dx*dx + dy*dy + dz*dz) / (2.0_dp * sigma * sigma)
                    kernal(i, j, k) = exp(gaussian_exp)
                end do
            end do
        end do
        !total_weight = sum(kernal)
        
        output = 0.0_dp
        ! smooth field with precomputed kernal
        ! use mask to avoid do loops
        ! add parallelization with dynamic schedule
        !$OMP PARALLEL DO PRIVATE( i,j,k) SHARED(tmp, output, sigma3)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    total_weight = sum(kernal, mask = tmp(i-sigma3:i+sigma3, j-sigma3:j+sigma3, k-sigma3:k+sigma3)>0.0_dp)
                    output(i, j, k) = sum(kernal*tmp(i-sigma3:i+sigma3, j-sigma3:j+sigma3, k-sigma3:k+sigma3))/total_weight
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        deallocate(kernal)
        deallocate(tmp)

    end subroutine gaussian_smooth2
    

function heaviside(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y

    if(x < 0.0_dp) then
        y = 0.0_dp
    else
        y = 1.0_dp
    endif
end function heaviside


subroutine hessian_at_point(density, bin_size, hessian)
    real(dp), intent(in) :: density(3, 3, 3), bin_size
    real(dp), intent(out) :: hessian(3, 3)

    real(dp) :: perturbation

    ! Define a small perturbation to ensure second derivatives aren't too close to zero
    perturbation = 0.0!1.0e-6_dp

    ! Calculate the second derivatives
    hessian(1,1) = (density(3,2,2) - 2.0_dp * density(2,2,2) + density(1,2,2)) / (bin_size * bin_size) + perturbation
    hessian(2,2) = (density(2,3,2) - 2.0_dp * density(2,2,2) + density(2,1,2)) / (bin_size * bin_size) + perturbation
    hessian(3,3) = (density(2,2,3) - 2.0_dp * density(2,2,2) + density(2,2,1)) / (bin_size * bin_size) + perturbation
    hessian(1,2) = (density(3,3,2) - density(3,1,2) - density(1,3,2) + density(1,1,2)) / (4.0_dp * bin_size * bin_size)
    hessian(1,3) = (density(3,2,3) - density(3,2,1) - density(1,2,3) + density(1,2,1)) / (4.0_dp * bin_size * bin_size)
    hessian(2,3) = (density(2,3,3) - density(2,3,1) - density(2,1,3) + density(2,1,1)) / (4.0_dp * bin_size * bin_size)

    ! Symmetrise the hessian
    hessian(2,1) = hessian(1,2)
    hessian(3,1) = hessian(1,3)
    hessian(3,2) = hessian(2,3)

end subroutine hessian_at_point

recursive subroutine flood_fill(volume, visited, x, y, z)
    logical, intent(in) :: volume(:,:,:)
    logical, intent(inout) :: visited(:,:,:)
    integer, intent(in) :: x, y, z
    integer :: dx, dy, dz

    if (.not. volume(x, y, z) .or. visited(x, y, z)) return

    visited(x, y, z) = .true.

    do dx = -1, 1
        do dy = -1, 1
            do dz = -1, 1
                if (x+dx > 0 .and. x+dx <= size(volume, dim=1) .and. &
                    y+dy > 0 .and. y+dy <= size(volume, dim=2) .and. &
                    z+dz > 0 .and. z+dz <= size(volume, dim=3)) then
                    call flood_fill(volume, visited, x+dx, y+dy, z+dz)
                end if
            end do
        end do
    end do
end subroutine flood_fill

function count_clusters(volume)
    implicit none
    logical, intent(in) :: volume(:,:,:)
    logical :: visited(size(volume, dim=1), size(volume, dim=2), size(volume, dim=3))
    integer :: count_clusters, x, y, z

    visited = .false.
    count_clusters = 0

    do x = 1, size(volume, dim=1)
        do y = 1, size(volume, dim=2)
            do z = 1, size(volume, dim=3)
                if (volume(x, y, z) .and. .not. visited(x, y, z)) then
                    count_clusters = count_clusters + 1
                    call flood_fill(volume, visited, x, y, z)
                end if
            end do
        end do
    end do
end function count_clusters

function count_voids(volume)
    implicit none
    logical, intent(in) :: volume(:,:,:)
    logical :: inverted_volume(size(volume, dim=1), size(volume, dim=2), size(volume, dim=3))
    integer :: count_voids

    inverted_volume = .not. volume
    count_voids = count_clusters(inverted_volume)
end function count_voids

function calculate_genus(volume)
    implicit none
    logical, intent(in) :: volume(:,:,:)
    integer :: calculate_genus

    calculate_genus = count_voids(volume) - count_clusters(volume)
end function calculate_genus

function field_to_volume(field, threshold)
    implicit none
    real(dp), intent(in) :: field(:,:,:)
    real(dp), intent(in) :: threshold
    logical :: field_to_volume(size(field, dim=1), size(field, dim=2), size(field, dim=3))
    integer :: i, j, k

    do i = 1, size(field, dim=1)
        do j = 1, size(field, dim=2)
            do k = 1, size(field, dim=3)
                field_to_volume(i, j, k) = (field(i, j, k) > threshold)
            end do
        end do
    end do
end function field_to_volume


subroutine daxpy ( n, da, dx, incx, dy, incy )
    !*****************************************************************************80
    !
    !! DAXPY computes constant times a vector plus a vector.
    !
    !  Discussion:
    !
    !    This routine uses double precision real arithmetic.
    !
    !    This routine uses unrolled loops for increments equal to one.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    16 May 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
    !    David Kincaid, Fred Krogh.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !    LINPACK User's Guide,
    !    SIAM, 1979,
    !    ISBN13: 978-0-898711-72-1,
    !    LC: QA214.L56.
    !
    !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    !    Algorithm 539, 
    !    Basic Linear Algebra Subprograms for Fortran Usage,
    !    ACM Transactions on Mathematical Software, 
    !    Volume 5, Number 3, September 1979, pages 308-323.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of elements in DX and DY.
    !
    !    Input, real ( kind = rk ) DA, the multiplier of DX.
    !
    !    Input, real ( kind = rk ) DX(*), the first vector.
    !
    !    Input, integer INCX, the increment between successive 
    !    entries of DX.
    !
    !    Input/output, real ( kind = rk ) DY(*), the second vector.
    !    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
    !
    !    Input, integer INCY, the increment between successive 
    !    entries of DY.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      real ( kind = rk ) da
      real ( kind = rk ) dx(*)
      real ( kind = rk ) dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n
    
      if ( n <= 0 ) then
        return
      end if
    
      if ( da == 0.0D+00 ) then
        return
      end if
    !
    !  Code for unequal increments or equal increments
    !  not equal to 1.
    !
      if ( incx /= 1 .or. incy /= 1 ) then
    
        if ( 0 <= incx ) then
          ix = 1
        else
          ix = ( - n + 1 ) * incx + 1
        end if
    
        if ( 0 <= incy ) then
          iy = 1
        else
          iy = ( - n + 1 ) * incy + 1
        end if
    
        do i = 1, n
          dy(iy) = dy(iy) + da * dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do
    !
    !  Code for both increments equal to 1.
    !
      else
    
        m = mod ( n, 4 )
    
        dy(1:m) = dy(1:m) + da * dx(1:m)
    
        do i = m+1, n, 4
          dy(i  ) = dy(i  ) + da * dx(i  )
          dy(i+1) = dy(i+1) + da * dx(i+1)
          dy(i+2) = dy(i+2) + da * dx(i+2)
          dy(i+3) = dy(i+3) + da * dx(i+3)
        end do
    
      end if
    
      return
    end subroutine daxpy

    function ddot ( n, dx, incx, dy, incy )
    !*****************************************************************************80
    !
    !! DDOT forms the dot product of two vectors.
    !
    !  Discussion:
    !
    !    This routine uses double precision real arithmetic.
    !
    !    This routine uses unrolled loops for increments equal to one.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    16 May 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
    !    David Kincaid, Fred Krogh.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !    LINPACK User's Guide,
    !    SIAM, 1979,
    !    ISBN13: 978-0-898711-72-1,
    !    LC: QA214.L56.
    !
    !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    !    Algorithm 539, 
    !    Basic Linear Algebra Subprograms for Fortran Usage,
    !    ACM Transactions on Mathematical Software, 
    !    Volume 5, Number 3, September 1979, pages 308-323.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of entries in the vectors.
    !
    !    Input, real ( kind = rk ) DX(*), the first vector.
    !
    !    Input, integer INCX, the increment between successive 
    !    entries in DX.
    !
    !    Input, real ( kind = rk ) DY(*), the second vector.
    !
    !    Input, integer INCY, the increment between successive 
    !    entries in DY.
    !
    !    Output, real ( kind = rk ) DDOT, the sum of the product of the 
    !    corresponding entries of DX and DY.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      real ( kind = rk ) ddot
      real ( kind = rk ) dtemp
      real ( kind = rk ) dx(*)
      real ( kind = rk ) dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n
    
      ddot = 0.0D+00
      dtemp = 0.0D+00
    
      if ( n <= 0 ) then
        return
      end if
    !
    !  Code for unequal increments or equal increments
    !  not equal to 1.
    !
      if ( incx /= 1 .or. incy /= 1 ) then
    
        if ( 0 <= incx ) then
          ix = 1
        else
          ix = ( - n + 1 ) * incx + 1
        end if
    
        if ( 0 <= incy ) then
          iy = 1
        else
          iy = ( - n + 1 ) * incy + 1
        end if
    
        do i = 1, n
          dtemp = dtemp + dx(ix) * dy(iy)
          ix = ix + incx
          iy = iy + incy
        end do
    !
    !  Code for both increments equal to 1.
    !
      else
    
        m = mod ( n, 5 )
    
        do i = 1, m
          dtemp = dtemp + dx(i) * dy(i)
        end do
    
        do i = m+1, n, 5
    
          dtemp = dtemp + dx(i  ) * dy(i  ) &
                        + dx(i+1) * dy(i+1) &
                        + dx(i+2) * dy(i+2) &
                        + dx(i+3) * dy(i+3) &
                        + dx(i+4) * dy(i+4)
        end do
    
      end if
    
      ddot = dtemp
    
      return
    end function ddot

    function dnrm2 ( n, x, incx )
    !*****************************************************************************80
    !
    !! DNRM2 returns the euclidean norm of a vector.
    !
    !  Discussion:
    !
    !    This routine uses double precision real arithmetic.
    !
    !     DNRM2 ( X ) = sqrt ( X' * X )
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    16 May 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
    !    David Kincaid, Fred Krogh.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !    LINPACK User's Guide,
    !    SIAM, 1979,
    !    ISBN13: 978-0-898711-72-1,
    !    LC: QA214.L56.
    !
    !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    !    Algorithm 539, 
    !    Basic Linear Algebra Subprograms for Fortran Usage,
    !    ACM Transactions on Mathematical Software,
    !    Volume 5, Number 3, September 1979, pages 308-323.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of entries in the vector.
    !
    !    Input, real ( kind = rk ) X(*), the vector whose norm is to be computed.
    !
    !    Input, integer INCX, the increment between successive 
    !    entries of X.
    !
    !    Output, real ( kind = rk ) DNRM2, the Euclidean norm of X.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      real ( kind = rk ) absxi
      real ( kind = rk ) dnrm2
      integer incx
      integer ix
      integer n
      real ( kind = rk ) norm
      real ( kind = rk ) scale
      real ( kind = rk ) ssq
      real ( kind = rk ) x(*)
    
      if ( n < 1 .or. incx < 1 ) then
    
        norm  = 0.0D+00
    
      else if ( n == 1 ) then
    
        norm  = abs ( x(1) )
    
      else
    
        scale = 0.0D+00
        ssq = 1.0D+00
    
        do ix = 1, 1 + ( n - 1 )*incx, incx
          if ( x(ix) /= 0.0D+00 ) then
            absxi = abs ( x(ix) )
            if ( scale < absxi ) then
              ssq = 1.0D+00 + ssq * ( scale / absxi )**2
              scale = absxi
            else
              ssq = ssq + ( absxi / scale )**2
            end if
          end if
        end do
        norm  = scale * sqrt ( ssq )
      end if
    
      dnrm2 = norm
    
      return
    end function dnrm2

    subroutine drot ( n, x, incx, y, incy, c, s )
    !*****************************************************************************80
    !
    !! DROT applies a plane rotation.
    !
    !  Discussion:
    !
    !    This routine uses double precision real arithmetic.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    08 April 1999
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
    !    David Kincaid, Fred Krogh.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !    LINPACK User's Guide,
    !    SIAM, 1979,
    !    ISBN13: 978-0-898711-72-1,
    !    LC: QA214.L56.
    !
    !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    !    Algorithm 539, 
    !    Basic Linear Algebra Subprograms for Fortran Usage,
    !    ACM Transactions on Mathematical Software,
    !    Volume 5, Number 3, September 1979, pages 308-323.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of entries in the vectors.
    !
    !    Input/output, real ( kind = rk ) X(*), one of the vectors to be rotated.
    !
    !    Input, integer INCX, the increment between successive 
    !    entries of X.
    !
    !    Input/output, real ( kind = rk ) Y(*), one of the vectors to be rotated.
    !
    !    Input, integer INCY, the increment between successive
    !    elements of Y.
    !
    !    Input, real ( kind = rk ) C, S, parameters (presumably the cosine and
    !    sine of some angle) that define a plane rotation.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      real ( kind = rk ) c
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer n
      real ( kind = rk ) s
      real ( kind = rk ) stemp
      real ( kind = rk ) x(*)
      real ( kind = rk ) y(*)
    
      if ( n <= 0 ) then
    
      else if ( incx == 1 .and. incy == 1 ) then
    
        do i = 1, n
          stemp = c * x(i) + s * y(i)
          y(i) = c * y(i) - s * x(i)
          x(i) = stemp
        end do
    
      else
    
        if ( 0 <= incx ) then
          ix = 1
        else
          ix = ( - n + 1 ) * incx + 1
        end if
    
        if ( 0 <= incy ) then
          iy = 1
        else
          iy = ( - n + 1 ) * incy + 1
        end if
    
        do i = 1, n
          stemp = c * x(ix) + s * y(iy)
          y(iy) = c * y(iy) - s * x(ix)
          x(ix) = stemp
          ix = ix + incx
          iy = iy + incy
        end do
    
      end if
    
      return
    end subroutine drot

    subroutine drotg ( sa, sb, c, s )
    !*****************************************************************************80
    !
    !! DROTG constructs a Givens plane rotation.
    !
    !  Discussion:
    !
    !    This routine uses double precision real arithmetic.
    !
    !    Given values A and B, this routine computes
    !
    !    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
    !          = sign ( B ) if abs ( A ) <= abs ( B );
    !
    !    R     = SIGMA * ( A * A + B * B );
    !
    !    C = A / R if R is not 0
    !      = 1     if R is 0;
    !
    !    S = B / R if R is not 0,
    !        0     if R is 0.
    !
    !    The computed numbers then satisfy the equation
    !
    !    (  C  S ) ( A ) = ( R )
    !    ( -S  C ) ( B ) = ( 0 )
    !
    !    The routine also computes
    !
    !    Z = S     if abs ( A ) > abs ( B ),
    !      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
    !      = 1     if C is 0.
    !
    !    The single value Z encodes C and S, and hence the rotation:
    !
    !    If Z = 1, set C = 0 and S = 1;
    !    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
    !    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    15 May 2006
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
    !    David Kincaid, Fred Krogh.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !    LINPACK User's Guide,
    !    SIAM, 1979,
    !    ISBN13: 978-0-898711-72-1,
    !    LC: QA214.L56.
    !
    !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    !    Algorithm 539, 
    !    Basic Linear Algebra Subprograms for Fortran Usage,
    !    ACM Transactions on Mathematical Software,
    !    Volume 5, Number 3, September 1979, pages 308-323.
    !
    !  Parameters:
    !
    !    Input/output, real ( kind = rk ) SA, SB.  On input, SA and SB are the values
    !    A and B.  On output, SA is overwritten with R, and SB is
    !    overwritten with Z.
    !
    !    Output, real ( kind = rk ) C, S, the cosine and sine of the
    !    Givens rotation.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      real ( kind = rk ) c
      real ( kind = rk ) r
      real ( kind = rk ) roe
      real ( kind = rk ) s
      real ( kind = rk ) sa
      real ( kind = rk ) sb
      real ( kind = rk ) scale
      real ( kind = rk ) z
    
      if ( abs ( sb ) < abs ( sa ) ) then
        roe = sa
      else
        roe = sb
      end if
    
      scale = abs ( sa ) + abs ( sb )
    
      if ( scale == 0.0D+00 ) then
        c = 1.0D+00
        s = 0.0D+00
        r = 0.0D+00
      else
        r = scale * sqrt ( ( sa / scale )**2 + ( sb / scale )**2 )
        r = sign ( 1.0D+00, roe ) * r
        c = sa / r
        s = sb / r
      end if
    
      if ( 0.0D+00 < abs ( c ) .and. abs ( c ) <= s ) then
        z = 1.0D+00 / c
      else
        z = s
      end if
    
      sa = r
      sb = z
    
      return
    end subroutine drotg

    subroutine dscal ( n, sa, x, incx )
    !*****************************************************************************80
    !
    !! DSCAL scales a vector by a constant.
    !
    !  Discussion:
    !
    !    This routine uses double precision real arithmetic.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    08 April 1999
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
    !    David Kincaid, Fred Krogh.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !    LINPACK User's Guide,
    !    SIAM, 1979,
    !    ISBN13: 978-0-898711-72-1,
    !    LC: QA214.L56.
    !
    !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    !    Algorithm 539, 
    !    Basic Linear Algebra Subprograms for Fortran Usage,
    !    ACM Transactions on Mathematical Software,
    !    Volume 5, Number 3, September 1979, pages 308-323.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of entries in the vector.
    !
    !    Input, real ( kind = rk ) SA, the multiplier.
    !
    !    Input/output, real ( kind = rk ) X(*), the vector to be scaled.
    !
    !    Input, integer INCX, the increment between successive 
    !    entries of X.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer i
      integer incx
      integer ix
      integer m
      integer n
      real ( kind = rk ) sa
      real ( kind = rk ) x(*)
    
      if ( n <= 0 ) then
    
      else if ( incx == 1 ) then
    
        m = mod ( n, 5 )
    
        x(1:m) = sa * x(1:m)
    
        do i = m+1, n, 5
          x(i)   = sa * x(i)
          x(i+1) = sa * x(i+1)
          x(i+2) = sa * x(i+2)
          x(i+3) = sa * x(i+3)
          x(i+4) = sa * x(i+4)
        end do
    
      else
    
        if ( 0 <= incx ) then
          ix = 1
        else
          ix = ( - n + 1 ) * incx + 1
        end if
    
        do i = 1, n
          x(ix) = sa * x(ix)
          ix = ix + incx
        end do
    
      end if
    
      return
    end subroutine dscal

    subroutine dsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )
    !*****************************************************************************80
    !
    !! DSVDC computes the singular value decomposition of a real rectangular matrix.
    !
    !  Discussion:
    !
    !    This routine reduces an M by N matrix A to diagonal form by orthogonal
    !    transformations U and V.  The diagonal elements S(I) are the singular
    !    values of A.  The columns of U are the corresponding left singular
    !    vectors, and the columns of V the right singular vectors.
    !
    !    The form of the singular value decomposition is then
    !
    !      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    16 September 2006
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Jack Dongarra, Jim Bunch, Cleve Moler, 
    !    Pete Stewart.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !    LINPACK User's Guide,
    !    SIAM, 1979,
    !    ISBN13: 978-0-898711-72-1,
    !    LC: QA214.L56.
    !
    !  Parameters:
    !
    !    Input/output, real ( kind = rk ) A(LDA,N).  On input, the M by N
    !    matrix whose singular value decomposition is to be computed.
    !    On output, the matrix has been destroyed.  Depending on the user's
    !    requests, the matrix may contain other useful information.
    !
    !    Input, integer LDA, the leading dimension of the array A.
    !    LDA must be at least N.
    !
    !    Input, integer M, the number of rows of the matrix.
    !
    !    Input, integer N, the number of columns of the matrix A.
    !
    !    Output, real ( kind = rk ) S(MM), where MM = max(M+1,N).  The first
    !    min(M,N) entries of S contain the singular values of A arranged in
    !    descending order of magnitude.
    !
    !    Output, real ( kind = rk ) E(MM), where MM = max(M+1,N).  Ordinarily
    !    contains zeros.  However see the discussion of INFO for exceptions.
    !
    !    Output, real ( kind = rk ) U(LDU,K).  If JOBA = 1 then K = M;
    !    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of
    !    left singular vectors.  U is not referenced if JOBA = 0.  If M <= N
    !    or if JOBA = 2, then U may be identified with A in the subroutine call.
    !
    !    Input, integer LDU, the leading dimension of the array U.
    !    LDU must be at least M.
    !
    !    Output, real ( kind = rk ) V(LDV,N), the N by N matrix of right singular
    !    vectors.  V is not referenced if JOB is 0.  If N <= M, then V may be
    !    identified with A in the subroutine call.
    !
    !    Input, integer LDV, the leading dimension of the array V.
    !    LDV must be at least N.
    !
    !    Workspace, real ( kind = rk ) WORK(M).
    !
    !    Input, integer JOB, controls the computation of the singular
    !    vectors.  It has the decimal expansion AB with the following meaning:
    !      A =  0, do not compute the left singular vectors.
    !      A =  1, return the M left singular vectors in U.
    !      A >= 2, return the first min(M,N) singular vectors in U.
    !      B =  0, do not compute the right singular vectors.
    !      B =  1, return the right singular vectors in V.
    !
    !    Output, integer INFO, status indicator.
    !    The singular values (and their corresponding singular vectors)
    !    S(INFO+1), S(INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
    !    Thus if INFO is 0, all the singular values and their vectors are
    !    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
    !    matrix with the elements of S on its diagonal and the elements of E on
    !    its superdiagonal.  Thus the singular values of A and B are the same.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer lda
      integer ldu
      integer ldv
      integer m
      integer n
    
      real ( kind = rk ) a(lda,n)
      real ( kind = rk ) b
      real ( kind = rk ) c
      real ( kind = rk ) cs
      real ( kind = rk ) e(*)
      real ( kind = rk ) el
      real ( kind = rk ) emm1
      real ( kind = rk ) f
      real ( kind = rk ) g
      integer info
      integer iter
      integer j
      integer job
      integer jobu
      integer k
      integer kase
      integer kk
      integer l
      integer ll
      integer lls
      integer ls
      integer lu
      integer, parameter :: maxit = 30
      integer mm
      integer mm1
      integer mn
      integer nct
      integer nctp1
      integer ncu
      integer nrt
      integer nrtp1
      real ( kind = rk ) s(*)
      real ( kind = rk ) scale
      !real ( kind = rk ) ddot
      real ( kind = rk ) shift
      real ( kind = rk ) sl
      real ( kind = rk ) sm
      real ( kind = rk ) smm1
      real ( kind = rk ) sn
      !real ( kind = rk ) dnrm2
      real ( kind = rk ) t
      real ( kind = rk ) t1
      real ( kind = rk ) test
      real ( kind = rk ) u(ldu,m)
      real ( kind = rk ) v(ldv,n)
      logical wantu
      logical wantv
      real ( kind = rk ) work(m)
      real ( kind = rk ) ztest
    !
    !  Determine what is to be computed.
    !
      wantu = .false.
      wantv = .false.
      jobu = mod ( job, 100 ) / 10
    
      if ( 1 < jobu ) then
        ncu = min ( m, n )
      else
        ncu = m
      end if
    
      if ( jobu /= 0 ) then
        wantu = .true.
      end if
    
      if ( mod ( job, 10 ) /= 0 ) then
        wantv = .true.
      end if
    !
    !  Reduce A to bidiagonal form, storing the diagonal elements
    !  in S and the super-diagonal elements in E.
    !
      info = 0
      nct = min ( m-1, n )
      nrt = max ( 0, min ( m, n-2 ) )
      lu = max ( nct, nrt )
    
      do l = 1, lu
    !
    !  Compute the transformation for the L-th column and
    !  place the L-th diagonal in S(L).
    !
        if ( l <= nct ) then
    
          s(l) = dnrm2 ( m-l+1, a(l,l), 1 )
    
          if ( s(l) /= 0.0D+00 ) then
            if ( a(l,l) /= 0.0D+00 ) then
              s(l) = sign ( s(l), a(l,l) )
            end if
            call dscal ( m-l+1, 1.0D+00 / s(l), a(l,l), 1 )
            a(l,l) = 1.0D+00 + a(l,l)
          end if
    
          s(l) = -s(l)
    
        end if
    
        do j = l+1, n
    !
    !  Apply the transformation.
    !
          if ( l <= nct .and. s(l) /= 0.0D+00 ) then
            t = -ddot ( m-l+1, a(l,l), 1, a(l,j), 1 ) / a(l,l)
            call daxpy ( m-l+1, t, a(l,l), 1, a(l,j), 1 )
          end if
    !
    !  Place the L-th row of A into E for the
    !  subsequent calculation of the row transformation.
    !
          e(j) = a(l,j)
    
        end do
    !
    !  Place the transformation in U for subsequent back multiplication.
    !
        if ( wantu .and. l <= nct ) then
          u(l:m,l) = a(l:m,l)
        end if
    !
    !  Compute the L-th row transformation and place the
    !  L-th superdiagonal in E(L).
    !
        if ( l <= nrt ) then
    
          e(l) = dnrm2 ( n-l, e(l+1), 1 )
    
          if ( e(l) /= 0.0D+00 ) then
            if ( e(l+1) /= 0.0D+00 ) then
              e(l) = sign ( e(l), e(l+1) )
            end if
            call dscal ( n-l, 1.0D+00 / e(l), e(l+1), 1 )
            e(l+1) = 1.0D+00 + e(l+1)
          end if
    
          e(l) = -e(l)
    !
    !  Apply the transformation.
    !
          if ( l + 1 <= m .and. e(l) /= 0.0D+00 ) then
    
            work(l+1:m) = 0.0D+00
    
            do j = l+1, n
              call daxpy ( m-l, e(j), a(l+1,j), 1, work(l+1), 1 )
            end do
    
            do j = l+1, n
              call daxpy ( m-l, -e(j)/e(l+1), work(l+1), 1, a(l+1,j), 1 )
            end do
    
          end if
    !
    !  Place the transformation in V for subsequent back multiplication.
    !
          if ( wantv ) then
            v(l+1:n,l) = e(l+1:n)
          end if
    
        end if
    
      end do
    !
    !  Set up the final bidiagonal matrix of order MN.
    !
      mn = min ( m + 1, n )
      nctp1 = nct + 1
      nrtp1 = nrt + 1
    
      if ( nct < n ) then
        s(nctp1) = a(nctp1,nctp1)
      end if
    
      if ( m < mn ) then
        s(mn) = 0.0D+00
      end if
    
      if ( nrtp1 < mn ) then
        e(nrtp1) = a(nrtp1,mn)
      end if
    
      e(mn) = 0.0D+00
    !
    !  If required, generate U.
    !
      if ( wantu ) then
    
        u(1:m,nctp1:ncu) = 0.0D+00
    
        do j = nctp1, ncu
          u(j,j) = 1.0D+00
        end do
    
        do ll = 1, nct
    
          l = nct - ll + 1
    
          if ( s(l) /= 0.0D+00 ) then
    
            do j = l+1, ncu
              t = -ddot ( m-l+1, u(l,l), 1, u(l,j), 1 ) / u(l,l)
              call daxpy ( m-l+1, t, u(l,l), 1, u(l,j), 1 )
            end do
    
            u(l:m,l) = -u(l:m,l)
            u(l,l) = 1.0D+00 + u(l,l)
            u(1:l-1,l) = 0.0D+00
    
          else
    
            u(1:m,l) = 0.0D+00
            u(l,l) = 1.0D+00
    
          end if
    
        end do
    
      end if
    !
    !  If it is required, generate V.
    !
      if ( wantv ) then
    
        do ll = 1, n
    
          l = n - ll + 1
    
          if ( l <= nrt .and. e(l) /= 0.0D+00 ) then
    
            do j = l + 1, n
              t = -ddot ( n-l, v(l+1,l), 1, v(l+1,j), 1 ) / v(l+1,l)
              call daxpy ( n-l, t, v(l+1,l), 1, v(l+1,j), 1 )
            end do
    
          end if
    
          v(1:n,l) = 0.0D+00
          v(l,l) = 1.0D+00
    
        end do
    
      end if
    !
    !  Main iteration loop for the singular values.
    !
      mm = mn
      iter = 0
    
      do while ( 0 < mn )
    !
    !  If too many iterations have been performed, set flag and return.
    !
        if ( maxit <= iter ) then
          info = mn
          return
        end if
    !
    !  This section of the program inspects for
    !  negligible elements in the S and E arrays.
    !
    !  On completion the variables KASE and L are set as follows:
    !
    !  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
    !  KASE = 2     if S(L) is negligible and L < MN
    !  KASE = 3     if E(L-1) is negligible, L < MN, and
    !               S(L), ..., S(MN) are not negligible (QR step).
    !  KASE = 4     if E(MN-1) is negligible (convergence).
    !
        do ll = 1, mn
    
          l = mn - ll
    
          if ( l == 0 ) then
            exit
          end if
    
          test = abs ( s(l) ) + abs ( s(l+1) )
          ztest = test + abs ( e(l) )
    
          if ( ztest == test ) then
            e(l) = 0.0D+00
            exit
          end if
    
        end do
    
        if ( l == mn - 1 ) then
    
          kase = 4
    
        else
    
          do lls = l + 1, mn + 1
    
            ls = mn - lls + l + 1
    
            if ( ls == l ) then
              exit
            end if
    
            test = 0.0D+00
            if ( ls /= mn ) then
              test = test + abs ( e(ls) )
            end if
    
            if ( ls /= l + 1 ) then
              test = test + abs ( e(ls-1) )
            end if
    
            ztest = test + abs ( s(ls) )
    
            if ( ztest == test ) then
              s(ls) = 0.0D+00
              exit
            end if
    
          end do
    
          if ( ls == l ) then
            kase = 3
          else if ( ls == mn ) then
            kase = 1
          else
            kase = 2
            l = ls
          end if
    
        end if
    
        l = l + 1
    !
    !  Deflate negligible S(MN).
    !
        if ( kase == 1 ) then
    
          mm1 = mn - 1
          f = e(mn-1)
          e(mn-1) = 0.0D+00
    
          do kk = l, mm1
    
            k = mm1 - kk + l
            t1 = s(k)
            call drotg ( t1, f, cs, sn )
            s(k) = t1
    
            if ( k /= l ) then
              f = -sn * e(k-1)
              e(k-1) = cs * e(k-1)
            end if
    
            if ( wantv ) then
              call drot ( n, v(1,k), 1, v(1,mn), 1, cs, sn )
            end if
    
          end do
    !
    !  Split at negligible S(L).
    !
        else if ( kase == 2 ) then
    
          f = e(l-1)
          e(l-1) = 0.0D+00
    
          do k = l, mn
    
            t1 = s(k)
            call drotg ( t1, f, cs, sn )
            s(k) = t1
            f = -sn * e(k)
            e(k) = cs * e(k)
            if ( wantu ) then
              call drot ( m, u(1,k), 1, u(1,l-1), 1, cs, sn )
            end if
    
          end do
    !
    !  Perform one QR step.
    !
        else if ( kase == 3 ) then
    !
    !  Calculate the shift.
    !
          scale = max ( abs ( s(mn) ), abs ( s(mn-1) ), abs ( e(mn-1) ), &
                        abs ( s(l) ), abs ( e(l) ) )
    
          sm = s(mn) / scale
          smm1 = s(mn-1) / scale
          emm1 = e(mn-1) / scale
          sl = s(l) / scale
          el = e(l) / scale
          b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0D+00
          c = sm  * sm * emm1 * emm1
          shift = 0.0D+00
    
          if ( b /= 0.0D+00 .or. c /= 0.0D+00 ) then
            shift = sqrt ( b * b + c )
            if ( b < 0.0D+00 ) then
              shift = -shift
            end if
            shift = c / ( b + shift )
          end if
    
          f = ( sl + sm ) * ( sl - sm ) + shift
          g = sl * el
    !
    !  Chase zeros.
    !
          mm1 = mn - 1
    
          do k = l, mm1
    
            call drotg ( f, g, cs, sn )
    
            if ( k /= l ) then
              e(k-1) = f
            end if
    
            f = cs * s(k) + sn * e(k)
            e(k) = cs * e(k) - sn * s(k)
            g = sn * s(k+1)
            s(k+1) = cs * s(k+1)
    
            if ( wantv ) then
              call drot ( n, v(1,k), 1, v(1,k+1), 1, cs, sn )
            end if
    
            call drotg ( f, g, cs, sn )
            s(k) = f
            f = cs * e(k) + sn * s(k+1)
            s(k+1) = -sn * e(k) + cs * s(k+1)
            g = sn * e(k+1)
            e(k+1) = cs * e(k+1)
    
            if ( wantu .and. k < m ) then
              call drot ( m, u(1,k), 1, u(1,k+1), 1, cs, sn )
            end if
    
          end do
    
          e(mn-1) = f
          iter = iter + 1
    !
    !  Convergence.
    !
        else if ( kase == 4 ) then
    !
    !  Make the singular value nonnegative.
    !
          if ( s(l) < 0.0D+00 ) then
            s(l) = -s(l)
            if ( wantv ) then
              v(1:n,l) = -v(1:n,l)
            end if
          end if
    !
    !  Order the singular value.
    !
          do
    
            if ( l == mm ) then
              exit
            end if
    
            if ( s(l+1) <= s(l) ) then
              exit
            end if
    
            t = s(l)
            s(l) = s(l+1)
            s(l+1) = t
    
            if ( wantv .and. l < n ) then
              call dswap ( n, v(1,l), 1, v(1,l+1), 1 )
            end if
    
            if ( wantu .and. l < m ) then
              call dswap ( m, u(1,l), 1, u(1,l+1), 1 )
            end if
    
            l = l + 1
    
          end do
    
          iter = 0
          mn = mn - 1
    
        end if
    
      end do
    
      return
    end subroutine dsvdc

    subroutine dswap ( n, x, incx, y, incy )
    !*****************************************************************************80
    !
    !! DSWAP interchanges two vectors.
    !
    !  Discussion:
    !
    !    This routine uses double precision real arithmetic.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    08 April 1999
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
    !    David Kincaid, Fred Krogh.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !    LINPACK User's Guide,
    !    SIAM, 1979,
    !    ISBN13: 978-0-898711-72-1,
    !    LC: QA214.L56.
    !
    !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    !    Algorithm 539, 
    !    Basic Linear Algebra Subprograms for Fortran Usage,
    !    ACM Transactions on Mathematical Software, 
    !    Volume 5, Number 3, September 1979, pages 308-323.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of entries in the vectors.
    !
    !    Input/output, real ( kind = rk ) X(*), one of the vectors to swap.
    !
    !    Input, integer INCX, the increment between successive 
    !    entries of X.
    !
    !    Input/output, real ( kind = rk ) Y(*), one of the vectors to swap.
    !
    !    Input, integer INCY, the increment between successive 
    !    elements of Y.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n
      real ( kind = rk ) temp
      real ( kind = rk ) x(*)
      real ( kind = rk ) y(*)
    
      if ( n <= 0 ) then
    
      else if ( incx == 1 .and. incy == 1 ) then
    
        m = mod ( n, 3 )
    
        do i = 1, m
          temp = x(i)
          x(i) = y(i)
          y(i) = temp
        end do
    
        do i = m+1, n, 3
    
          temp = x(i)
          x(i) = y(i)
          y(i) = temp
    
          temp = x(i+1)
          x(i+1) = y(i+1)
          y(i+1) = temp
    
          temp = x(i+2)
          x(i+2) = y(i+2)
          y(i+2) = temp
    
        end do
    
      else
    
        if ( 0 <= incx ) then
          ix = 1
        else
          ix = ( - n + 1 ) * incx + 1
        end if
    
        if ( 0 <= incy ) then
          iy = 1
        else
          iy = ( - n + 1 ) * incy + 1
        end if
    
        do i = 1, n
          temp = x(ix)
          x(ix) = y(iy)
          y(iy) = temp
          ix = ix + incx
          iy = iy + incy
        end do
    
      end if
    
      return
    end subroutine dswap

    subroutine phi1 ( n, r, r0, v )
    !*****************************************************************************80
    !
    !! PHI1 evaluates the multiquadric radial basis function.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 June 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    !    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    !    Third Edition,
    !    Cambridge University Press, 2007,
    !    ISBN13: 978-0-521-88068-8,
    !    LC: QA297.N866.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of points.
    !
    !    Input, real ( kind = rk ) R(N), the radial separation.
    !    0 < R.
    !
    !    Input, real ( kind = rk ) R0, a scale factor.
    !
    !    Output, real ( kind = rk ) V(N), the value of the radial basis function.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer n
    
      real ( kind = rk ) r(n)
      real ( kind = rk ) r0
      real ( kind = rk ) v(n)
    
      v(1:n) = sqrt ( r(1:n)**2 + r0**2 )
    
      return
    end subroutine phi1

    subroutine phi2 ( n, r, r0, v )
    !*****************************************************************************80
    !
    !! PHI2 evaluates the inverse multiquadric radial basis function.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 June 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    !    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    !    Third Edition,
    !    Cambridge University Press, 2007,
    !    ISBN13: 978-0-521-88068-8,
    !    LC: QA297.N866.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of points.
    !
    !    Input, real ( kind = rk ) R(N), the radial separation.
    !    0 < R.
    !
    !    Input, real ( kind = rk ) R0, a scale factor.
    !
    !    Output, real ( kind = rk ) V(N), the value of the radial basis function.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer n
    
      real ( kind = rk ) r(n)
      real ( kind = rk ) r0
      real ( kind = rk ) v(n)
    
      v = 1.0D+00 / sqrt ( r**2 + r0**2 )
    
      return
    end subroutine phi2

    subroutine phi3 ( n, r, r0, v )
    
    !*****************************************************************************80
    !
    !! PHI3 evaluates the thin-plate spline radial basis function.
    !
    !  Discussion:
    !
    !    Note that PHI3(R,R0) is negative if R < R0.  Thus, for this basis function,
    !    it may be desirable to choose a value of R0 smaller than any possible R.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 June 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    !    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    !    Third Edition,
    !    Cambridge University Press, 2007,
    !    ISBN13: 978-0-521-88068-8,
    !    LC: QA297.N866.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of points.
    !
    !    Input, real ( kind = rk ) R(N), the radial separation.
    !    0 < R.
    !
    !    Input, real ( kind = rk ) R0, a scale factor.
    !
    !    Output, real ( kind = rk ) V(N), the value of the radial basis function.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer n
    
      integer i
      real ( kind = rk ) r(n)
      real ( kind = rk ) r0
      real ( kind = rk ) v(n)
    
      do i = 1, n
        if ( r(i) .le. 0.0D+00 ) then
          v(i) = 0.0D+00
        else
          v(i) = r(i)**2 * log ( r(i) / r0 )
        end if
      end do
    
      return
    end 

    subroutine phi4 ( n, r, r0, v )
    !*****************************************************************************80
    !
    !! PHI4 evaluates the gaussian radial basis function.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 June 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    !    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    !    Third Edition,
    !    Cambridge University Press, 2007,
    !    ISBN13: 978-0-521-88068-8,
    !    LC: QA297.N866.
    !
    !  Parameters:
    !
    !    Input, integer N, the number of points.
    !
    !    Input, real ( kind = rk ) R(N), the radial separation.
    !    0 < R.
    !
    !    Input, real ( kind = rk ) R0, a scale factor.
    !
    !    Output, real ( kind = rk ) V(N), the value of the radial basis function.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer n
    
      real ( kind = rk ) r(n)
      real ( kind = rk ) r0
      real ( kind = rk ) v(n)
    
      v(1:n) = exp ( - 0.5D+00 * r(1:n)**2 / r0**2 )
    
      return
    end subroutine phi4

    subroutine phi5 ( n, r, r0, v )
        !*****************************************************************************80
        !
        !! PHI5 evaluates the linear basis function.
        !
        !  Discussion:
        !
        !    The linear basis function is also known as the "thin plate spline"
        !    basis function.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    30 June 2012
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Reference:
        !
        !    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
        !    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
        !    Third Edition,
        !    Cambridge University Press, 2007,
        !    ISBN13: 978-0-521-88068-8,
        !    LC: QA297.N866.
        !
        !  Parameters:
        !
        !    Input, integer N, the number of points.
        !
        !    Input, real ( kind = rk ) R(N), the radial separation.
        !    0 < R.
        !
        !    Output, real ( kind = rk ) V(N), the value of the radial basis function.
        !
        implicit none

        integer, parameter :: rk = kind ( 1.0D+00 )

        integer n

        integer i

        real ( kind = rk ) r(n)
        real ( kind = rk ) r0
        real ( kind = rk ) v(n)

        v(1:n) = r(1:n)

        return
    end subroutine phi5


    subroutine r8mat_solve_svd ( m, n, a, b, x )
    
    !*****************************************************************************80
    !
    !! R8MAT_SOLVE_SVD solves a linear system A*x=b using the SVD.
    !
    !  Discussion:
    !
    !    When the system is determined, the solution is the solution in the
    !    ordinary sense, and A*x = b.
    !
    !    When the system is overdetermined, the solution minimizes the
    !    L2 norm of the residual ||A*x-b||.
    !
    !    When the system is underdetermined, ||A*x-b|| should be zero, and
    !    the solution is the solution of minimum L2 norm, ||x||.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 June 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer M, N, the number of rows and columns
    !    in the matrix A.
    !
    !    Input, real ( kind = rk ) A(M,N), the matrix.
    !
    !    Input, real ( kind = rk ) B(M), the right hand side.
    !
    !    Output, real ( kind = rk ) X(N), the solution.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer m
      integer n
    
      real ( kind = rk ) a(m,n)
      real ( kind = rk ) a_copy(m,n)
      real ( kind = rk ) a_pseudo(n,m)
      real ( kind = rk ) b(m)
      real ( kind = rk ) e(max(m+1,n))
      integer i
      integer info
      integer lda
      integer ldu
      integer ldv
      integer job
      real ( kind = rk ) s(m,n)
      real ( kind = rk ) sp(n,m)
      real ( kind = rk ) sdiag(max(m+1,n))
      real ( kind = rk ) u(m,m)
      real ( kind = rk ) v(n,n)
      real ( kind = rk ) work(m)
      real ( kind = rk ) x(n)
    !
    !  Compute the SVD decomposition.
    !
      job = 11
      lda = m
      ldu = m
      ldv = n
    
      a_copy(1:m,1:n) = a(1:m,1:n)
    
      call dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job, info )
    
      if ( info /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SOLVE_SVD - Failure!'
        write ( *, '(a)' ) '  The SVD could not be calculated.'
        write ( *, '(a)' ) '  LINPACK routine DSVDC returned a nonzero'
        write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
        return
      end if
    
      s(1:m,1:n) = 0.0D+00
      do i = 1, min ( m, n )
        s(i,i) = sdiag(i)
      end do
    !
    !  Compute the pseudo inverse.
    !
      sp(1:n,1:m) = 0.0D+00
      do i = 1, min ( m, n )
        if ( s(i,i) /= 0.0D+00 ) then
          sp(i,i) = 1.0D+00 / s(i,i)
        end if
      end do
    
      a_pseudo(1:n,1:m) = matmul ( v(1:n,1:n), &
        matmul ( sp(1:n,1:m), transpose ( u(1:m,1:m) ) ) )
    !
    !  Compute x = A_pseudo * b
    !
      x(1:n) = matmul ( a_pseudo(1:n,1:m), b(1:m) )
    
      return
    end subroutine r8mat_solve_svd

    subroutine rbf_interp_nd ( m, nd, xd, r0, phi, w, ni, xi, fi )
    !*****************************************************************************80
    !
    !! RBF_INTERP_ND evaluates a radial basis function interpolant.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 June 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    !    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    !    Third Edition,
    !    Cambridge University Press, 2007,
    !    ISBN13: 978-0-521-88068-8,
    !    LC: QA297.N866.
    !
    !  Parameters:
    !
    !    Input, integer M, the spatial dimension.
    !
    !    Input, integer ND, the number of data points.
    !
    !    Input, real ( kind = rk ) XD(M,ND), the data points.
    !
    !    Input, real ( kind = rk ) R0, a scale factor.  R0 should be larger than 
    !    the typical separation between points, but smaller than the maximum 
    !    separation.  The value of R0 has a significant effect on the resulting 
    !    interpolant.
    !
    !    Input, subroutine PHI ( N, R, R0, V ), a subroutine to evaluate the radial
    !    basis functions.
    !
    !    Input, real ( kind = rk ) W(ND), the weights, as computed by RBF_WEIGHTS.
    !
    !    Input, integer NI, the number of interpolation points.
    !
    !    Input, real ( kind = rk ) XI(M,NI), the interpolation points.
    !
    !    Output, real ( kind = rk ) FI(NI), the interpolated values.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer m
      integer nd
      integer ni
    
      real ( kind = rk ) fi(ni)
      integer i
      integer j
      external phi
      real ( kind = rk ) r(nd)
      real ( kind = rk ) r0
      real ( kind = rk ) v(nd)
      real ( kind = rk ) w(nd)
      real ( kind = rk ) xd(m,nd)
      real ( kind = rk ) xi(m,ni)
    
      do i = 1, ni
    
        do j = 1, nd
          r(j) = sqrt ( sum ( ( xi(1:m,i) - xd(1:m,j) )**2 ) )
        end do
    
        call phi ( nd, r, r0, v )
    
        fi(i) = dot_product ( v, w )
    
      end do
    
      return
    end subroutine rbf_interp_nd

    subroutine rbf_weight ( m, nd, xd, r0, phi, fd, w )
    
    !*****************************************************************************80
    !
    !! RBF_WEIGHT computes weights for radial basis function interpolation.
    !
    !  Discussion:
    !
    !    We assume that there are N (nonsingular) equations in N unknowns.
    !
    !    However, it should be clear that, if we are willing to do some kind
    !    of least squares calculation, we could allow for singularity,
    !    inconsistency, or underdetermine systems.  This could be associated
    !    with data points that are very close or repeated, a smaller number
    !    of data points than function values, or some other ill-conditioning
    !    of the system arising from a peculiarity in the point spacing.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 June 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    !    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    !    Third Edition,
    !    Cambridge University Press, 2007,
    !    ISBN13: 978-0-521-88068-8,
    !    LC: QA297.N866.
    !
    !  Parameters:
    !
    !    Input, integer M, the spatial dimension.
    !
    !    Input, integer ND, the number of data points.
    !
    !    Input, real ( kind = rk ) XD(M,ND), the data points.
    !
    !    Input, real ( kind = rk ) R0, a scale factor.  R0 should be larger than 
    !    the typical separation between points, but smaller than the maximum 
    !    separation.  The value of R0 has a significant effect on the resulting 
    !    interpolant.
    !
    !    Input, subroutine PHI ( N, R, R0, V ), a subroutine to evaluate the radial
    !    basis functions.
    !
    !    Input, real ( kind = rk ) FD(ND), the function values at the data points.
    !
    !    Output, real ( kind = rk ) W(ND), the weights.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer m
      integer nd
    
      real ( kind = rk ) a(nd,nd)
      real ( kind = rk ) fd(nd)
      integer i
      integer j
      external phi
      real ( kind = rk ) r(nd)
      real ( kind = rk ) r0
      real ( kind = rk ) v(nd)
      real ( kind = rk ) w(nd)
      real ( kind = rk ) xd(m,nd)
    
    !$OMP PARALLEL DO  schedule(dynamic) PRIVATE ( i, j, r, v ), SHARED ( a, fd, nd, r0, xd )
      do i = 1, nd
        do j = 1, nd
          r(j) = sqrt ( sum ( ( xd(1:m,i) - xd(1:m,j) )**2 ) )
        end do
    
        call phi ( nd, r, r0, v )
    
        a(i,1:nd) = v(1:nd)
    
      end do
    !$OMP END PARALLEL DO
    !
    !  Solve for the weights.
    !
      call r8mat_solve_svd ( nd, nd, a, fd, w )
    
      return
    end subroutine rbf_weight

end program density_histogram

