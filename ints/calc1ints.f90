#include "xkind.h"
#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use types
    implicit none
    mwSize, intent(in) :: nlhs, nrhs
    mwPointer, intent(in), dimension(*) :: prhs
    mwPointer, intent(out), dimension(*) :: plhs
    mwIndex :: i
    mwSize :: nshell, natoms, mxGetN, mxGetM, contr, nbasis
    mwPointer :: H, S, kin, coul, harm, mxGetField, mxGetPr, &
                 mxCreateDoubleMatrix, mxDuplicateArray
    real*8 :: mxGetScalar
    integer :: l
    type(Basis_func_t), dimension(:), allocatable :: basis
    type(Atom_t), dimension(:), allocatable :: atoms

    nshell = mxGetN(prhs(1))
    natoms = mxGetN(prhs(2))
    allocate (basis(nshell), atoms(natoms))
    nbasis = 0
    do i = 1, nshell
        call mxCopyPtrToReal8(mxGetPr(mxGetField(prhs(1), i, 'R')), basis(i)%R, 3)
        l = mxGetScalar(mxGetField(prhs(1), i, 'l'))
        basis(i)%l = l
        contr = mxGetM(mxGetField(prhs(1), i, 'zeta'))
        basis(i)%contr = contr
        allocate (basis(i)%zeta(contr), basis(i)%d(contr))
        call mxCopyPtrToReal8(mxGetPr(mxGetField(prhs(1), i, 'zeta')), basis(i)%zeta, contr)
        call mxCopyPtrToReal8(mxGetPr(mxGetField(prhs(1), i, 'd')), basis(i)%d, contr)
        nbasis = nbasis+(l+1)*(l+2)/2
    end do
    do i = 1, natoms
        atoms(i)%Z = mxGetScalar(mxGetField(prhs(2), i, 'Z'))
        call mxCopyPtrToReal8(mxGetPr(mxGetField(prhs(2), i, 'R')), atoms(i)%R, 3)
    end do
    plhs(1) = mxCreateDoubleMatrix(nbasis, nbasis, 0)
    plhs(2) = mxCreateDoubleMatrix(nbasis, nbasis, 0)
    plhs(3) = mxCreateDoubleMatrix(nbasis, nbasis, 0)
    plhs(4) = mxCreateDoubleMatrix(nbasis, nbasis, 0)
    plhs(5) = mxCreateDoubleMatrix(nbasis, nbasis, 0)
    H = mxGetPr(plhs(1))
    S = mxGetPr(plhs(2))
    kin = mxGetPr(plhs(3))
    coul = mxGetPr(plhs(4))
    harm = mxGetPr(plhs(5))
    call calc1ints(%val(H), %val(S), %val(kin), %val(coul), %val(harm), &
                   nbasis, basis, nshell, atoms, natoms)
    plhs(6) = mxDuplicateArray(prhs(1))
    do i = 1, nshell
        call mxCopyReal8ToPtr(basis(i)%d, mxGetPr(mxGetField(plhs(6), i, 'd')), &
                              basis(i)%contr)
        deallocate (basis(i)%zeta, basis(i)%d)
    end do
    deallocate (basis, atoms)

    return
end

subroutine calc1ints(S, H, kin, coul, harm, nbasis, basis, nshell, atoms, natoms)
    use types
    use gen1int
    implicit none
    mwSize, intent(in) :: nbasis, nshell, natoms
    type(Basis_func_t), dimension(nshell), intent(inout) :: basis
    type(Atom_t), dimension(natoms), intent(in) :: atoms
    real(REALK), dimension(nbasis, nbasis), intent(out) :: S, H, kin, coul, harm
    integer :: i, j, k
    mwIndex :: s0, s1
    integer :: n0, n1
    integer :: l0, l1
    type(one_prop_t) :: one_prop_S, one_prop_H, one_prop_kin, one_prop_coul, &
                        one_prop_harm
    integer :: num_prop_S, num_prop_H, num_prop_kin, num_prop_coul, &
               num_prop_harm
    real(REALK), dimension(:,:,:,:,:), allocatable :: &
        contr_ints_S, contr_ints_H, contr_ints_kin, contr_ints_coul, &
        contr_ints_harm
    integer :: ierr
    integer, dimension(:), allocatable :: idx
    real(REALK), dimension(3, natoms) :: coords
    real(REALK), dimension(natoms) :: charges
    integer, dimension(:,:), allocatable :: powers_bra, powers_ket
    allocate (idx(nshell))
    do i = 1, nshell
        if (i == 1) then
            idx(i) = 1
        else
            idx(i) = idx(i-1)+(basis(i-1)%l+1)*(basis(i-1)%l+2)/2
        end if
        call norm_contr_cgto(basis(i)%l, basis(i)%contr, &
                             basis(i)%zeta, 1, basis(i)%d)
    end do
    do i = 1, natoms
        coords(:, i) = atoms(i)%R
        charges(i) = -atoms(i)%Z
    end do
    call OnePropCreate(prop_name=INT_OVERLAP, one_prop=one_prop_S, info_prop=ierr)
    call OnePropGetNumProp(one_prop=one_prop_S, num_prop=num_prop_S)
    call OnePropCreate(prop_name=INT_ONE_HAMIL, one_prop=one_prop_H, &
                       coord_nuclei=coords, charge_nuclei=charges, info_prop=ierr)
    call OnePropGetNumProp(one_prop=one_prop_H, num_prop=num_prop_H)
    call OnePropCreate(prop_name=INT_KIN_ENERGY, one_prop=one_prop_kin, info_prop=ierr)
    call OnePropGetNumProp(one_prop=one_prop_kin, num_prop=num_prop_kin)
    call OnePropCreate(prop_name=INT_POT_ENERGY, one_prop=one_prop_coul, &
                       coord_nuclei=coords, charge_nuclei=charges, info_prop=ierr)
    call OnePropGetNumProp(one_prop=one_prop_H, num_prop=num_prop_H)
    do s0 = 1, nshell
    do s1 = 1, nshell
        l0 = basis(s0)%l
        l1 = basis(s1)%l
        if (l0 < l1) then
            cycle
        end if
        n0 = (l0+1)*(l0+2)/2
        n1 = (l1+1)*(l1+2)/2
        allocate (contr_ints_S(n0, 1, n1, 1, num_prop_S))
        allocate (contr_ints_H(n0, 1, n1, 1, num_prop_H))
        allocate (contr_ints_kin(n0, 1, n1, 1, num_prop_kin))
        allocate (contr_ints_coul(n0, 1, n1, 1, num_prop_coul))
        allocate (powers_bra(3, n0), powers_ket(3, n1))
        k = 1
        do i = 0, l0
            do j = 0, i
                powers_bra(:, k) = (/l0-i, i-j, j/)
                k = k+1
            end do
        end do
        k = 1
        do i = 0, l1
            do j = 0, i
                powers_ket(:, k) = (/l1-i, i-j, j/)
                k = k+1
            end do
        end do
        call OnePropGetIntegral(idx_bra=1, coord_bra=basis(s0)%R, &
                                angular_bra=basis(s0)%l, &
                                num_prim_bra=basis(s0)%contr, &
                                exponent_bra=basis(s0)%zeta, &
                                num_contr_bra=1, &
                                contr_coef_bra=basis(s0)%d, &
                                idx_ket=1, coord_ket=basis(s1)%R, &
                                angular_ket=basis(s1)%l, &
                                num_prim_ket=basis(s1)%contr, &
                                exponent_ket=basis(s1)%zeta, &
                                num_contr_ket=1, &
                                contr_coef_ket=basis(s1)%d, &
                                one_prop=one_prop_S, &
                                num_gto_bra=n0, num_gto_ket=n1, &
                                num_opt=num_prop_S, contr_ints=contr_ints_S, &
                                spher_gto=.false., &
                                powers_bra=powers_bra, powers_ket=powers_ket)
        call OnePropGetIntegral(idx_bra=1, coord_bra=basis(s0)%R, &
                                angular_bra=basis(s0)%l, &
                                num_prim_bra=basis(s0)%contr, &
                                exponent_bra=basis(s0)%zeta, &
                                num_contr_bra=1, &
                                contr_coef_bra=basis(s0)%d, &
                                idx_ket=1, coord_ket=basis(s1)%R, &
                                angular_ket=basis(s1)%l, &
                                num_prim_ket=basis(s1)%contr, &
                                exponent_ket=basis(s1)%zeta, &
                                num_contr_ket=1, &
                                contr_coef_ket=basis(s1)%d, &
                                one_prop=one_prop_H, &
                                num_gto_bra=n0, num_gto_ket=n1, &
                                num_opt=num_prop_H, contr_ints=contr_ints_H, &
                                spher_gto=.false., &
                                powers_bra=powers_bra, powers_ket=powers_ket)
        call OnePropGetIntegral(idx_bra=1, coord_bra=basis(s0)%R, &
                                angular_bra=basis(s0)%l, &
                                num_prim_bra=basis(s0)%contr, &
                                exponent_bra=basis(s0)%zeta, &
                                num_contr_bra=1, &
                                contr_coef_bra=basis(s0)%d, &
                                idx_ket=1, coord_ket=basis(s1)%R, &
                                angular_ket=basis(s1)%l, &
                                num_prim_ket=basis(s1)%contr, &
                                exponent_ket=basis(s1)%zeta, &
                                num_contr_ket=1, &
                                contr_coef_ket=basis(s1)%d, &
                                one_prop=one_prop_kin, &
                                num_gto_bra=n0, num_gto_ket=n1, &
                                num_opt=num_prop_kin, contr_ints=contr_ints_kin, &
                                spher_gto=.false., &
                                powers_bra=powers_bra, powers_ket=powers_ket)
        call OnePropGetIntegral(idx_bra=1, coord_bra=basis(s0)%R, &
                                angular_bra=basis(s0)%l, &
                                num_prim_bra=basis(s0)%contr, &
                                exponent_bra=basis(s0)%zeta, &
                                num_contr_bra=1, &
                                contr_coef_bra=basis(s0)%d, &
                                idx_ket=1, coord_ket=basis(s1)%R, &
                                angular_ket=basis(s1)%l, &
                                num_prim_ket=basis(s1)%contr, &
                                exponent_ket=basis(s1)%zeta, &
                                num_contr_ket=1, &
                                contr_coef_ket=basis(s1)%d, &
                                one_prop=one_prop_coul, &
                                num_gto_bra=n0, num_gto_ket=n1, &
                                num_opt=num_prop_coul, contr_ints=contr_ints_coul, &
                                spher_gto=.false., &
                                powers_bra=powers_bra, powers_ket=powers_ket)
        do i = 1, n0
        do j = 1, n1
            S(idx(s0)+i-1, idx(s1)+j-1) = contr_ints_S(i, 1, j, 1, 1)
            S(idx(s1)+j-1, idx(s0)+i-1) = contr_ints_S(i, 1, j, 1, 1)
            H(idx(s0)+i-1, idx(s1)+j-1) = contr_ints_H(i, 1, j, 1, 1)
            H(idx(s1)+j-1, idx(s0)+i-1) = contr_ints_H(i, 1, j, 1, 1)
            kin(idx(s0)+i-1, idx(s1)+j-1) = contr_ints_kin(i, 1, j, 1, 1)
            kin(idx(s1)+j-1, idx(s0)+i-1) = contr_ints_kin(i, 1, j, 1, 1)
            coul(idx(s0)+i-1, idx(s1)+j-1) = contr_ints_coul(i, 1, j, 1, 1)
            coul(idx(s1)+j-1, idx(s0)+i-1) = contr_ints_coul(i, 1, j, 1, 1)
        end do
        end do
        deallocate (powers_bra, powers_ket)
        deallocate (contr_ints_S)
        deallocate (contr_ints_H)
        deallocate (contr_ints_kin)
        deallocate (contr_ints_coul)
    end do
    end do
    call OnePropDestroy(one_prop=one_prop_S)
    call OnePropDestroy(one_prop=one_prop_H)
    call OnePropDestroy(one_prop=one_prop_kin)
    call OnePropDestroy(one_prop=one_prop_coul)
    return
end
