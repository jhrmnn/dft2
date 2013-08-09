#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use types
    implicit none
    mwSize, intent(in) :: nlhs, nrhs
    mwPointer, intent(in), dimension(*) :: prhs
    mwPointer, intent(out), dimension(*) :: plhs
    mwSize :: mxGetN, mxGetM
    mwPointer :: mxGetField, mxGetPr, mxCreateDoubleMatrix
    real*8 :: mxGetScalar
    integer*4 :: mxGetString, stat
    integer :: i, l
    mwSize :: nbasis, ngrid, nshell, contr
    mwPointer :: P, grid, w, E, F
    character(len=50) :: func
    logical :: fderiv
    type(Basis_func_t), dimension(:), allocatable :: basis

    nbasis = mxGetM(prhs(1))
    P = mxGetPr(prhs(1))
    stat = mxGetString(prhs(2), func, mxGetN(prhs(2)))
    ngrid = mxGetM(mxGetField(prhs(3), 1, 'w'))
    grid = mxGetPr(mxGetField(prhs(3), 1, 'xyz'))
    w = mxGetPr(mxGetField(prhs(3), 1, 'w'))
    nshell = mxGetN(prhs(4))
    allocate (basis(nshell))
    do i = 1, nshell
        call mxCopyPtrToReal8(mxGetPr(mxGetField(prhs(4), i, 'R')), basis(i)%R, 3)
        l = mxGetScalar(mxGetField(prhs(4), i, 'l'))
        basis(i)%l = l
        contr = mxGetM(mxGetField(prhs(4), i, 'zeta'))
        basis(i)%contr = contr
        allocate (basis(i)%zeta(contr), basis(i)%d(contr))
        call mxCopyPtrToReal8(mxGetPr(mxGetField(prhs(4), i, 'zeta')), basis(i)%zeta, contr)
        call mxCopyPtrToReal8(mxGetPr(mxGetField(prhs(4), i, 'd')), basis(i)%d, contr)
    end do

    fderiv = (nlhs == 2)
    plhs(1) = mxCreateDoubleMatrix(1, 1, 0)
    E = mxGetPr(plhs(1))
    if (fderiv) then
        plhs(2) = mxCreateDoubleMatrix(nbasis, nbasis, 0)
        F = mxGetPr(plhs(2))
    end if

    call eval_xc(%val(E), %val(F), &
                 %val(P), %val(grid), %val(w), basis, &
                 nbasis, ngrid, nshell, &
                 func, fderiv)

    do i = 1, nshell
        deallocate (basis(i)%zeta, basis(i)%d)
    end do
    deallocate (basis)

    return
end

subroutine eval_xc(E, F, P, grid, w, basis, nbasis, ngrid, nshell, func, fderiv)
    use types
    implicit none
    integer, intent(in) :: nbasis, ngrid, nshell
    logical, intent(in) :: fderiv
    double precision, intent(in) :: P(nbasis, nbasis), grid(ngrid, 3),  w(ngrid)
    type(Basis_func_t), intent(in) :: basis(nshell)
    character(len=50), intent(in) :: func
    double precision, intent(out) :: E(1), F(nbasis, nbasis)
    integer, parameter :: batch = 1000
    double precision :: phi(batch, nbasis), dphi(batch, nbasis, 3)
    integer, allocatable :: idx(:)
    integer :: nbatch, i, npts

    nbatch = (ngrid-1)/batch+1
    allocate (idx(nbatch+1))
    idx = (/ (batch*(i-1)+1, i = 1, nbatch), ngrid+1 /)
    do i = 1, nbatch
        npts = idx(i+1)-idx(i)
        call eval_basis(phi, dphi, grid(idx(i):idx(i+1)-1, :), &
                        basis, npts, nbasis, nshell)
        call eval_func(E, F, P, phi, dphi, w(idx(i):idx(i+1)-1), &
                       nbasis, npts, &
                       func, fderiv)
    end do
    deallocate (idx)

    return
end

integer function l_to_n(l)
    implicit none
    integer, intent(in) :: l
    l_to_n = (l+1)*(l+2)/2
end function

subroutine eval_basis(phi, dphi, grid, basis, ngrid, nbasis, nshell)
    use types
    implicit none
    mwSize, intent(in) :: ngrid, nbasis, nshell
    type(Basis_func_t), intent(in) :: basis(nshell)
    double precision, intent(in) :: grid(ngrid, 3)
    double precision, intent(out) :: phi(ngrid, nbasis), dphi(ngrid, nbasis, 3)
    integer :: ibas(nshell), l_to_n, ncgto, i, j, k, ii, jj, kk
    type(Basis_func_t) :: b
    integer, allocatable :: n(:, :)
    double precision, allocatable :: t4(:)
    double precision :: R(3), R2, t1, t2(3), t3, t5(3), t6

    do i = 1, nshell
        ibas(i) = merge(0, ibas(i-1)+l_to_n(basis(i-1)%l), i == 1)
    end do
    do i = 1, nshell
        b = basis(i)
        ncgto = l_to_n(b%l)
        allocate (n(ncgto, 3), t4(b%contr))
        kk = 1
        do ii = 0, b%l
            do jj = 0, ii
                n(kk, :) = (/ b%l-ii, ii-jj, jj /)
                kk = kk+1
            end do
        end do
        do j = 1, ngrid
            R = grid(j, :)-b%R
            R2 = sum(R**2)
            t4 = b%d*exp(-b%zeta*R2)
            t1 = sum(t4)
            t3 = sum(b%zeta*t4)
            do k = 1, ncgto
                t2 = R**n(k, :)
                t6 = product(t2)
                phi(j, ibas(i)+k) = t6*t1
                t5 = -2*t6*R*t3
                if (n(k, 1) > 0) then
                    t5(1) = t5(1)+n(k, 1)*R(1)**(n(k, 1)-1)*R(2)*R(3)*t1
                end if
                if (n(k, 2) > 0) then
                    t5(2) = t5(2)+n(k, 2)*R(1)*R(2)**(n(k, 2)-1)*R(3)*t1
                end if
                if (n(k, 3) > 0) then
                    t5(3) = t5(3)+n(k, 3)*R(1)*R(2)*R(3)**(n(k, 3)-1)*t1
                end if
                dphi(j, ibas(i)+k, :) = t5
            end do
        end do
        deallocate (n, t4)
    end do

    return
end

subroutine eval_func(E, F, P, phi, dphi, w, nbasis, ngrid, func, fderiv)
    implicit none
    integer, intent(in) :: nbasis, ngrid
    double precision, intent(inout) :: E(1), F(nbasis, nbasis)
    double precision, intent(in) :: P(nbasis, nbasis), phi(ngrid, nbasis), &
                                    dphi(ngrid, nbasis, 3), w(ngrid)
    character(len=50), intent(in) :: func
    logical, intent(in) :: fderiv
    double precision :: phiP(nbasis), rhoc(ngrid), sumphiPdphi(ngrid, 3), &
                        sigmacc(ngrid), zk(ngrid), vrhoc(ngrid), &
                        vsigmacc(ngrid), drhoc, dsigmacc
    integer :: i, j, k

    do i = 1, ngrid
        phiP = matmul(phi(i, :), P)
        rhoc(i) = sum(phiP*phi(i, :))
        do j = 1, 3
            sumphiPdphi(i, j) = sum(phiP*dphi(i, :, j))
        end do
        sigmacc(i) = 4*sum(sumphiPdphi**2)
        zk(i) = 0
        vrhoc(i) = 0
        vsigmacc(i) = 0
    end do
    call xc_func(rhoc, sigmacc, zk, vrhoc, vsigmacc, ngrid, func, fderiv)
    E(1) = E(1)+sum(zk*w)
    if (fderiv) then
        do i = 1, ngrid
            do j = 1, nbasis
                do k = 1, nbasis
                    drhoc = phi(i, j)*phi(i, k)
                    dsigmacc = 8*sum(sumphiPdphi(i, :)*phi(i, j)*dphi(i, k, :))
                    F(j, k) = F(j, k)+(drhoc*vrhoc(i)+dsigmacc*vsigmacc(i))*w(i)
                end do
            end do
        end do
    end if

    return
end
