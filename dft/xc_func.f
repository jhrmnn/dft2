      subroutine xc_func(rhoc, sigmacc, zk, vrhoc, vsigmacc, npt, 
     >                   key, fderiv)

      implicit double precision (a-h, o-z)
      parameter (maxn=1000)
      dimension rhoc(maxn),rhoo(1)
      dimension sigmacc(maxn), sigmaco(1), sigmaoo(1)
      dimension tauc(1), tauo(1), upsilonc(1),upsilono(1)
      dimension zk(maxn), vrhoc(maxn), vrhoo(1)
      dimension vsigmacc(maxn), vsigmaco(1), vsigmaoo(1)
      dimension vtauc(1), vtauo(1), vupsilonc(1), vupsilono(1)
      logical fderiv, open
      integer igrad, npt
      character*32 name, key

      pi=acos(-1d0)
      open = .false.

      if (.false.) then
          write (*,*) 'This really should not happen'
c:dfauto
      end if

      end
