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
c:DENSITYcallstart
        else if (key.eq.'DENSITY') then
        call dftacg_density(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:DENSITYcallend
c:DIRACcallstart
        else if (key.eq.'DIRAC') then
        call dftacg_dirac(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:DIRACcallend
c:PW92Ccallstart
        else if (key.eq.'PW92C') then
        call dftacg_pw92c(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,tauc,tauo,
     >                   upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,vtauc,vtauo,
     >                   vupsilonc,vupsilono)
c:PW92Ccallend
c:dfauto
      end if

      end
c:DENSITYsubrstart

c    Generated: Fri  9 Aug 2013 16:53:09 CEST

      subroutine dftacg_density
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated DENSITY'
      igrad=0
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

             rho = rhoa + rhob
      zk(i) = rho
      vrhoc(i) = vrhoc(i) + 0.1D1
      vrhoo(i) = vrhoo(i)

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

             rho = rhoa + rhob
      zk(i) = rho

             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      zk(i) = rho
      vrhoc(i) = vrhoc(i) + 0.1D1

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      zk(i) = rho

             endif
           enddo
         endif
       endif

      return
      end

c:DENSITYsubrend
c:DIRACsubrstart

c    Generated: Fri  9 Aug 2013 16:53:09 CEST

      subroutine dftacg_dirac
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated DIRAC'
      igrad=0
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

             rho = rhoa + rhob
      t2 = (0.1D1 / pi) ** (0.1D1 / 0.3D1)
      t3 = rhoa ** (0.1D1 / 0.3D1)
      t7 = rhob ** (0.1D1 / 0.3D1)
      zk(i) = -0.136284044462410474416D1 * t2 * t3 * rhoa - 0.1362840444
     #62410474416D1 * t2 * t7 * rhob
      t13 = 0.908560296416069829442D0 * t2 * t3
      t15 = 0.908560296416069829442D0 * t2 * t7
      vrhoc(i) = vrhoc(i) - t13 - t15
      vrhoo(i) = vrhoo(i) - t13 + t15

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

             rho = rhoa + rhob
      t2 = (0.1D1 / pi) ** (0.1D1 / 0.3D1)
      t3 = rhoa ** (0.1D1 / 0.3D1)
      t7 = rhob ** (0.1D1 / 0.3D1)
      zk(i) = -0.136284044462410474416D1 * t2 * t3 * rhoa - 0.1362840444
     #62410474416D1 * t2 * t7 * rhob

             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      t4 = (0.1D1 / pi) ** (0.1D1 / 0.3D1)
      t5 = rhoa ** (0.1D1 / 0.3D1)
      t9 = rhob ** (0.1D1 / 0.3D1)
      zk(i) = -0.136284044462410474416D1 * t4 * t5 * rhoa - 0.1362840444
     #62410474416D1 * t4 * t9 * rhob
      vrhoc(i) = vrhoc(i) - 0.908560296416069829442D0 * t4 * t5 - 0.9085
     #60296416069829442D0 * t4 * t9

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      t4 = (0.1D1 / pi) ** (0.1D1 / 0.3D1)
      t5 = rhoa ** (0.1D1 / 0.3D1)
      t9 = rhob ** (0.1D1 / 0.3D1)
      zk(i) = -0.136284044462410474416D1 * t4 * t5 * rhoa - 0.1362840444
     #62410474416D1 * t4 * t9 * rhob

             endif
           enddo
         endif
       endif

      return
      end

c:DIRACsubrend
c:PW92Csubrstart

c    Generated: Fri  9 Aug 2013 16:53:10 CEST

      subroutine dftacg_pw92c
     > (name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   tauc,tauo,upsilonc,upsilono,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,
     >                   vtauc,vtauo,vupsilonc,vupsilono)
      implicit double precision (a-h,o-z)
      logical fderiv,open
      integer igrad,npt
      character*(*) name
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision tauc(*),tauo(*)
      double precision upsilonc(*),upsilono(*)
      double precision zk(*),vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)
      double precision vtauc(*),vtauo(*)
      double precision vupsilonc(*),vupsilono(*)
      include "common/cdft"
      include "common/tapes"
      parameter(tol=1d-12)
      pi=acos(-1d0)
      name='Automatically generated PW92C'
      igrad=0
       if(open) then
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

             rho = rhoa + rhob
      t1 = pi * rho
      t2 = t1 ** (0.1D1 / 0.3D1)
      t4 = 0.1D1 + 0.194159335344114122552D0 * t2
      t5 = t1 ** (0.1D1 / 0.6D1)
      t8 = sqrt(t1)
      t10 = t2 ** 2
      t12 = 0.724010193431683113327D1 * t5 + 0.325955091942229212011D1 *
     # t2 + 0.141872281647966739112D1 * t8 + 0.406913004517529319387D0 *
     # t10
      t15 = 0.1D1 + 0.160818243221511048214D2 / t12
      t16 = log(t15)
      t18 = 0.62182D-1 * t4 * t16
      t20 = 0.1D1 + 0.101077332976287768525D0 * t2
      t25 = 0.987212972256927209438D1 * t5 + 0.329180480994506259905D1 *
     # t2 + 0.762327521935289963194D0 * t8 + 0.410025070949612505036D0 *
     # t10
      t28 = 0.1D1 + 0.296085746432166755492D2 / t25
      t29 = log(t28)
      t30 = t20 * t29
      t32 = rhoa - 0.1D1 * rhob
      t33 = 0.1D1 / rho
      t34 = t32 * t33
      t35 = 0.1D1 + t34
      t36 = t35 ** (0.1D1 / 0.3D1)
      t39 = 0.1D1 - 0.1D1 * t34
      t40 = t39 ** (0.1D1 / 0.3D1)
      t42 = t36 * t35 + t40 * t39 - 0.2D1
      t43 = t32 ** 2
      t44 = t43 ** 2
      t45 = rho ** 2
      t46 = t45 ** 2
      t47 = 0.1D1 / t46
      t48 = t44 * t47
      t50 = 0.1D1 - 0.1D1 * t48
      t53 = 0.379957485370152817954D-1 * t30 * t42 * t50
      t55 = 0.1D1 + 0.186690969707574028554D0 * t2
      t60 = 0.134579137143944477912D2 * t5 + 0.563098414909787598194D1 *
     # t2 + 0.291521471421917737271D1 * t8 + 0.516066464547863440989D0 *
     # t10
      t63 = 0.1D1 + 0.321646831778706979736D2 / t60
      t64 = log(t63)
      t67 = -0.31090D-1 * t55 * t64 + t18
      t68 = t67 * t42
      t70 = 0.192366105093153631974D1 * t68 * t48
      zk(i) = rho * (-t18 + t53 + t70)
      t74 = 0.1D1 / t10 * pi
      t76 = 0.402440526345590145616D-2 * t74 * t16
      t77 = t12 ** 2
      t80 = t5 ** 2
      t81 = t80 ** 2
      t84 = 0.1D1 / t81 / t5 * pi
      t88 = 0.1D1 / t8 * pi
      t91 = 0.1D1 / t2 * pi
      t97 = 0.100000000000000000000D1 * t4 / t77 * (0.120668365571947185
     #555D1 * t84 + 0.108651697314076404004D1 * t74 + 0.7093614082398336
     #95563D0 * t88 + 0.271275336345019546258D0 * t91) / t15
      t101 = 0.128016964218639749325D-2 * t74 * t29 * t42 * t50
      t102 = t25 ** 2
      t115 = 0.112499995668310776915D1 * t20 / t102 * (0.164535495376154
     #534908D1 * t84 + 0.109726826998168753302D1 * t74 + 0.3811637609676
     #44981600D0 * t88 + 0.273350047299741670024D0 * t91) / t28 * t42 * 
     #t50
      t117 = t32 / t45
      t118 = 0.1D1 * t117
      t122 = 0.1D1 * t33
      t126 = 0.133333333333333333333D1 * t36 * (t33 - t118) + 0.13333333
     #3333333333333D1 * t40 * (-t122 + t117)
      t131 = t43 * t32 * t47
      t132 = 0.4D1 * t131
      t135 = t44 / t46 / rho
      t136 = 0.4D1 * t135
      t143 = t60 ** 2
      t158 = 0.192366105093153631974D1 * (-0.193474074940282551591D-2 * 
     #t74 * t64 + 0.999999999999999999999D0 * t55 / t143 * (0.2242985619
     #06574129855D1 * t84 + 0.187699471636595866065D1 * t74 + 0.14576073
     #5710958868637D1 * t88 + 0.344044309698575627326D0 * t91) / t63 + t
     #76 - t97) * t42 * t48
      t163 = 0.769464420372614527896D1 * t68 * t131
      t165 = 0.769464420372614527896D1 * t68 * t135
      t168 = 0.500000000000000000000D0 * rho * (-t76 + t97 + t101 - t115
     # + 0.379957485370152817954D-1 * t30 * t126 * t50 + 0.3799574853701
     #52817954D-1 * t30 * t42 * (-t132 + t136) + t158 + 0.19236610509315
     #3631974D1 * t67 * t126 * t48 + t163 - t165)
      t175 = 0.133333333333333333333D1 * t36 * (-t122 - t118) + 0.133333
     #333333333333333D1 * t40 * (t33 + t117)
      t188 = 0.500000000000000000000D0 * rho * (-t76 + t97 + t101 - t115
     # + 0.379957485370152817954D-1 * t30 * t175 * t50 + 0.3799574853701
     #52817954D-1 * t30 * t42 * (t132 + t136) + t158 + 0.192366105093153
     #631974D1 * t67 * t175 * t48 - t163 - t165)
      vrhoc(i) = vrhoc(i) + t168 + t188 - t18 + t53 + t70
      vrhoo(i) = vrhoo(i) + t168 - t188

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             t2 = 0.500000000000000000000D0 * rhoc(i)
      t4 = 0.500000000000000000000D0 * rhoo(i)
      rhoa = max(0.0D0, t2 + t4)
      rhob = max(0.0D0, t2 - t4)

             rho = rhoa + rhob
      t1 = pi * rho
      t2 = t1 ** (0.1D1 / 0.3D1)
      t5 = t1 ** (0.1D1 / 0.6D1)
      t8 = sqrt(t1)
      t10 = t2 ** 2
      t16 = log(0.1D1 + 0.160818243221511048214D2 / (0.72401019343168311
     #3327D1 * t5 + 0.325955091942229212011D1 * t2 + 0.14187228164796673
     #9112D1 * t8 + 0.406913004517529319387D0 * t10))
      t18 = 0.62182D-1 * (0.1D1 + 0.194159335344114122552D0 * t2) * t16
      t29 = log(0.1D1 + 0.296085746432166755492D2 / (0.98721297225692720
     #9438D1 * t5 + 0.329180480994506259905D1 * t2 + 0.76232752193528996
     #3194D0 * t8 + 0.410025070949612505036D0 * t10))
      t32 = rhoa - 0.1D1 * rhob
      t34 = t32 / rho
      t35 = 0.1D1 + t34
      t36 = t35 ** (0.1D1 / 0.3D1)
      t39 = 0.1D1 - 0.1D1 * t34
      t40 = t39 ** (0.1D1 / 0.3D1)
      t42 = t36 * t35 + t40 * t39 - 0.2D1
      t43 = t32 ** 2
      t44 = t43 ** 2
      t45 = rho ** 2
      t46 = t45 ** 2
      t48 = t44 / t46
      t64 = log(0.1D1 + 0.321646831778706979736D2 / (0.13457913714394447
     #7912D2 * t5 + 0.563098414909787598194D1 * t2 + 0.29152147142191773
     #7271D1 * t8 + 0.516066464547863440989D0 * t10))
      zk(i) = rho * (-t18 + 0.379957485370152817954D-1 * (0.1D1 + 0.1010
     #77332976287768525D0 * t2) * t29 * t42 * (0.1D1 - 0.1D1 * t48) + 0.
     #192366105093153631974D1 * (-0.31090D-1 * (0.1D1 + 0.18669096970757
     #4028554D0 * t2) * t64 + t18) * t42 * t48)

             endif
           enddo
         endif
       else
         if(fderiv) then
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      t3 = pi * rho
      t4 = t3 ** (0.1D1 / 0.3D1)
      t6 = 0.1D1 + 0.194159335344114122552D0 * t4
      t7 = t3 ** (0.1D1 / 0.6D1)
      t10 = sqrt(t3)
      t12 = t4 ** 2
      t14 = 0.724010193431683113327D1 * t7 + 0.325955091942229212011D1 *
     # t4 + 0.141872281647966739112D1 * t10 + 0.406913004517529319387D0 
     #* t12
      t17 = 0.1D1 + 0.160818243221511048214D2 / t14
      t18 = log(t17)
      t20 = 0.62182D-1 * t6 * t18
      t22 = 0.1D1 + 0.101077332976287768525D0 * t4
      t27 = 0.987212972256927209438D1 * t7 + 0.329180480994506259905D1 *
     # t4 + 0.762327521935289963194D0 * t10 + 0.410025070949612505036D0 
     #* t12
      t30 = 0.1D1 + 0.296085746432166755492D2 / t27
      t31 = log(t30)
      t32 = t22 * t31
      t34 = rhoa - 0.1D1 * rhob
      t35 = 0.1D1 / rho
      t36 = t34 * t35
      t37 = 0.1D1 + t36
      t38 = t37 ** (0.1D1 / 0.3D1)
      t41 = 0.1D1 - 0.1D1 * t36
      t42 = t41 ** (0.1D1 / 0.3D1)
      t44 = t38 * t37 + t42 * t41 - 0.2D1
      t45 = t34 ** 2
      t46 = t45 ** 2
      t47 = rho ** 2
      t48 = t47 ** 2
      t49 = 0.1D1 / t48
      t50 = t46 * t49
      t52 = 0.1D1 - 0.1D1 * t50
      t55 = 0.379957485370152817954D-1 * t32 * t44 * t52
      t57 = 0.1D1 + 0.186690969707574028554D0 * t4
      t62 = 0.134579137143944477912D2 * t7 + 0.563098414909787598194D1 *
     # t4 + 0.291521471421917737271D1 * t10 + 0.516066464547863440989D0 
     #* t12
      t65 = 0.1D1 + 0.321646831778706979736D2 / t62
      t66 = log(t65)
      t69 = -0.31090D-1 * t57 * t66 + t20
      t70 = t69 * t44
      t72 = 0.192366105093153631974D1 * t70 * t50
      zk(i) = rho * (-t20 + t55 + t72)
      t76 = 0.1D1 / t12 * pi
      t78 = 0.402440526345590145616D-2 * t76 * t18
      t79 = t14 ** 2
      t82 = t7 ** 2
      t83 = t82 ** 2
      t86 = 0.1D1 / t83 / t7 * pi
      t90 = 0.1D1 / t10 * pi
      t93 = 0.1D1 / t4 * pi
      t99 = 0.100000000000000000000D1 * t6 / t79 * (0.120668365571947185
     #555D1 * t86 + 0.108651697314076404004D1 * t76 + 0.7093614082398336
     #95563D0 * t90 + 0.271275336345019546258D0 * t93) / t17
      t103 = 0.128016964218639749325D-2 * t76 * t31 * t44 * t52
      t104 = t27 ** 2
      t117 = 0.112499995668310776915D1 * t22 / t104 * (0.164535495376154
     #534908D1 * t86 + 0.109726826998168753302D1 * t76 + 0.3811637609676
     #44981600D0 * t90 + 0.273350047299741670024D0 * t93) / t30 * t44 * 
     #t52
      t119 = t34 / t47
      t120 = 0.1D1 * t119
      t124 = 0.1D1 * t35
      t128 = 0.133333333333333333333D1 * t38 * (t35 - t120) + 0.13333333
     #3333333333333D1 * t42 * (-t124 + t119)
      t133 = t45 * t34 * t49
      t134 = 0.4D1 * t133
      t137 = t46 / t48 / rho
      t138 = 0.4D1 * t137
      t145 = t62 ** 2
      t160 = 0.192366105093153631974D1 * (-0.193474074940282551591D-2 * 
     #t76 * t66 + 0.999999999999999999999D0 * t57 / t145 * (0.2242985619
     #06574129855D1 * t86 + 0.187699471636595866065D1 * t76 + 0.14576073
     #5710958868637D1 * t90 + 0.344044309698575627326D0 * t93) / t65 + t
     #78 - t99) * t44 * t50
      t165 = 0.769464420372614527896D1 * t70 * t133
      t167 = 0.769464420372614527896D1 * t70 * t137
      t177 = 0.133333333333333333333D1 * t38 * (-t124 - t120) + 0.133333
     #333333333333333D1 * t42 * (t35 + t119)
      vrhoc(i) = vrhoc(i) + 0.5D0 * rho * (-t78 + t99 + t103 - t117 + 0.
     #379957485370152817954D-1 * t32 * t128 * t52 + 0.379957485370152817
     #954D-1 * t32 * t44 * (-t134 + t138) + t160 + 0.1923661050931536319
     #74D1 * t69 * t128 * t50 + t165 - t167) + 0.5D0 * rho * (-t78 + t99
     # + t103 - t117 + 0.379957485370152817954D-1 * t32 * t177 * t52 + 0
     #.379957485370152817954D-1 * t32 * t44 * (t134 + t138) + t160 + 0.1
     #92366105093153631974D1 * t69 * t177 * t50 - t165 - t167) - t20 + t
     #55 + t72

             endif
           enddo
         else
           do i=1,npt
             zk(i)=0.0d0
             if(rhoc(i).gt.tol) then
             rhoa = max(0.0D0, 0.500000000000000000000D0 * rhoc(i))
      rhob = rhoa
      rho = rhoa + rhob
      t3 = pi * rho
      t4 = t3 ** (0.1D1 / 0.3D1)
      t7 = t3 ** (0.1D1 / 0.6D1)
      t10 = sqrt(t3)
      t12 = t4 ** 2
      t18 = log(0.1D1 + 0.160818243221511048214D2 / (0.72401019343168311
     #3327D1 * t7 + 0.325955091942229212011D1 * t4 + 0.14187228164796673
     #9112D1 * t10 + 0.406913004517529319387D0 * t12))
      t20 = 0.62182D-1 * (0.1D1 + 0.194159335344114122552D0 * t4) * t18
      t31 = log(0.1D1 + 0.296085746432166755492D2 / (0.98721297225692720
     #9438D1 * t7 + 0.329180480994506259905D1 * t4 + 0.76232752193528996
     #3194D0 * t10 + 0.410025070949612505036D0 * t12))
      t34 = rhoa - 0.1D1 * rhob
      t36 = t34 / rho
      t37 = 0.1D1 + t36
      t38 = t37 ** (0.1D1 / 0.3D1)
      t41 = 0.1D1 - 0.1D1 * t36
      t42 = t41 ** (0.1D1 / 0.3D1)
      t44 = t38 * t37 + t42 * t41 - 0.2D1
      t45 = t34 ** 2
      t46 = t45 ** 2
      t47 = rho ** 2
      t48 = t47 ** 2
      t50 = t46 / t48
      t66 = log(0.1D1 + 0.321646831778706979736D2 / (0.13457913714394447
     #7912D2 * t7 + 0.563098414909787598194D1 * t4 + 0.29152147142191773
     #7271D1 * t10 + 0.516066464547863440989D0 * t12))
      zk(i) = rho * (-t20 + 0.379957485370152817954D-1 * (0.1D1 + 0.1010
     #77332976287768525D0 * t4) * t31 * t44 * (0.1D1 - 0.1D1 * t50) + 0.
     #192366105093153631974D1 * (-0.31090D-1 * (0.1D1 + 0.18669096970757
     #4028554D0 * t4) * t66 + t20) * t44 * t50)

             endif
           enddo
         endif
       endif

      return
      end

c:PW92Csubrend
