#include "xkind.h"
#include "fintrf.h"

module types
    implicit none
    type Basis_func_t
        real(REALK), dimension(3) :: R
        integer :: l
        mwSize :: contr
        real(REALK), dimension(:), allocatable :: zeta
        real(REALK), dimension(:), allocatable :: d
    end type
    type Atom_t
        integer :: Z
        real(REALK), dimension(3) :: R
    end type
end

