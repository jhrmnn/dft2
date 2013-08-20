#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use types
    implicit none
    mwSize, intent(in) :: nlhs, nrhs
    mwPointer, intent(in), dimension(*) :: prhs
    mwPointer, intent(out), dimension(*) :: plhs
    mwSize :: mxGetN, mxGetM
    mwPointer :: mxGetField, mxGetPr, mxCreateDoubleMatrix, &
                 mxGetFieldByNumber
    real*8 :: mxGetScalar
    integer*4 :: mxGetNumberOfFields, nfuncs
    integer :: i, l
    mwSize :: nbasis, ngrid, nshell, contr
    mwPointer :: P, grid, w, E, F
    character(len=50) :: mxGetFieldNameByNumber
    character(len=50), allocatable :: funcs(:)
    double precision, allocatable :: func_coeff(:)
    logical :: fderiv
    type(Basis_func_t), dimension(:), allocatable :: basis

    nbasis = mxGetM(prhs(1))
    P = mxGetPr(prhs(1))
    nfuncs = mxGetNumberOfFields(prhs(2))
    allocate (funcs(nfuncs), func_coeff(nfuncs))
    do i = 1, nfuncs
        funcs(i) = mxGetFieldNameByNumber(prhs(2), i)
        func_coeff(i) = mxGetScalar(mxGetFieldByNumber(prhs(2), 1, i))
    end do
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
    plhs(1) = mxCreateDoubleMatrix(1, nfuncs, 0)
    E = mxGetPr(plhs(1))
    if (fderiv) then
        plhs(2) = mxCreateDoubleMatrix(nbasis, nbasis, 0)
        F = mxGetPr(plhs(2))
    end if

    call eval_xc(%val(E), %val(F), &
                 %val(P), %val(grid), %val(w), basis, &
                 nbasis, ngrid, nshell, nfuncs, &
                 funcs, func_coeff, fderiv)

    deallocate (funcs, func_coeff)
    do i = 1, nshell
        deallocate (basis(i)%zeta, basis(i)%d)
    end do
    deallocate (basis)

    return
end

subroutine eval_xc(E, F, P, grid, w, basis, nbasis, ngrid, nshell, nfuncs, &
                  funcs, func_coeff, fderiv)
    use types
    implicit none
    integer, intent(in) :: nbasis, ngrid, nshell, nfuncs
    logical, intent(in) :: fderiv
    double precision, intent(in) :: P(nbasis, nbasis), grid(ngrid, 3), &
                                    w(ngrid), func_coeff(nfuncs)
    type(Basis_func_t), intent(in) :: basis(nshell)
    character(len=50), intent(in) :: funcs(nfuncs)
    double precision, intent(out) :: E(nfuncs), F(nbasis, nbasis)
    integer, parameter :: batch = 1000
    double precision :: phi(batch, nbasis), dphi(batch, nbasis, 3)
    integer, allocatable :: idx(:)
    integer :: nbatch, i, npts
    real :: t1, t2

    nbatch = (ngrid-1)/batch+1
    allocate (idx(nbatch+1))
    idx = (/ (batch*(i-1)+1, i = 1, nbatch), ngrid+1 /)
    do i = 1, nbatch
        npts = idx(i+1)-idx(i)
        call eval_basis(phi, dphi, grid(idx(i):idx(i+1)-1, :), &
                        basis, npts, nbasis, nshell)
        call eval_func(E, F, P, phi, dphi, w(idx(i):idx(i+1)-1), &
                       nbasis, npts, nfuncs, &
                       funcs, func_coeff, fderiv)
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
    double precision, allocatable :: prims(:)
    double precision :: R(3), R2, s1, s2, dphi0(3), ang

    do i = 1, nshell
        ibas(i) = merge(0, ibas(i-1)+l_to_n(basis(i-1)%l), i == 1)
    end do
    do i = 1, nshell
        b = basis(i)
        ncgto = l_to_n(b%l)
        allocate (n(ncgto, 3), prims(b%contr))
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
            prims = b%d*exp(-b%zeta*R2)
            s1 = sum(prims)
            s2 = sum(b%zeta*prims)
            do k = 1, ncgto
                ang = product(R**n(k, :))
                phi(j, ibas(i)+k) = ang*s1
                dphi0 = -2*ang*R*s2
                if (n(k, 1) > 0) then
                    dphi0(1) = dphi0(1)+n(k, 1)*R(1)**(n(k, 1)-1)*R(2)*R(3)*s1
                end if
                if (n(k, 2) > 0) then
                    dphi0(2) = dphi0(2)+n(k, 2)*R(1)*R(2)**(n(k, 2)-1)*R(3)*s1
                end if
                if (n(k, 3) > 0) then
                    dphi0(3) = dphi0(3)+n(k, 3)*R(1)*R(2)*R(3)**(n(k, 3)-1)*s1
                end if
                dphi(j, ibas(i)+k, :) = dphi0
            end do
        end do
        deallocate (n, prims)
    end do

    return
end

subroutine eval_func(E, F, P, phi, dphi, w, nbasis, ngrid, nfuncs, funcs, &
                     func_coeff, fderiv)
    implicit none
    integer, intent(in) :: nbasis, ngrid, nfuncs
    double precision, intent(inout) :: E(nfuncs), F(nbasis, nbasis)
    double precision, intent(in) :: P(nbasis, nbasis), phi(ngrid, nbasis), &
                                    dphi(ngrid, nbasis, 3), w(ngrid), &
                                    func_coeff(nfuncs)
    character(len=50), intent(in) :: funcs(nfuncs)
    logical, intent(in) :: fderiv
    double precision :: phiP(nbasis), rhoc(ngrid), sumphiPdphi(ngrid, 3), &
                        sigmacc(ngrid), zk(ngrid), vrhoc(ngrid), vrhoc0(ngrid), &
                        vsigmacc(ngrid), vsigmacc0(ngrid), drhoc, dsigmacc, &
                        rhoc0, sigmacc0, rfac(ngrid, nbasis)
    integer :: i, j, k
    real :: t1, t2

    do i = 1, ngrid
        phiP = matmul(phi(i, :), P)
        rhoc(i) = sum(phiP*phi(i, :))
        do j = 1, 3
            sumphiPdphi(i, j) = sum(phiP*dphi(i, :, j))
        end do
        sigmacc(i) = 4*sum(sumphiPdphi**2)
    end do
    vrhoc = 0
    vsigmacc = 0
    do i = 1, nfuncs
        vrhoc0 = 0
        vsigmacc0 = 0
        call xc_func(rhoc, sigmacc, zk, vrhoc0, vsigmacc0, ngrid, funcs(i), fderiv)
        E(i) = E(i)+func_coeff(i)*sum(zk*w)
        vrhoc = vrhoc+func_coeff(i)*vrhoc0
        vsigmacc = vsigmacc+func_coeff(i)*vsigmacc0
    end do
    if (fderiv) then
        do i = 1, ngrid
            do j = 1, nbasis
                rfac(i, j) = (phi(i, j)*vrhoc(i) &
                    +8*sum(sumphiPdphi(i, :)*dphi(i, j, :))*vsigmacc(i))*w(i)
            end do
        end do
        F = F+matmul(transpose(phi), rfac)
    end if

    return
end
