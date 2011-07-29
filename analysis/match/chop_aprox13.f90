
      subroutine burner(tstep,tin,din,ein,xin,tout,eout,xout)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'

! declare the pass
      double precision tstep,tin,din,ein,xin(1),tout,eout,xout(1)


! first we try a hydrostatic burn
      one_step             = .true.
      self_heat_const_den  = .false.
      call burner_driver(tstep,tin,din,ein,xin,tout,eout,xout)


! if the change in energy is large and negative,
! then do a self-heating burn instead

      if (eout .le. 0.7d0 * ein) then
       one_step             = .false.
       self_heat_const_den  = .true.
       call burner_driver(tstep,tin,din,ein,xin,tout,eout,xout)
      end if


      return
      end







      subroutine burner_driver(tstep,tin,din,ein,xin,tout,eout,xout)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'


! aprox13 in self-heating at constant density mode.
! given time step tstep, temperature tin, density din, internal energy ein and composition xin, 
! this routine returns the temperature tout, internal energy eout and burned composition xout.
! isotopes: he4, c12, o16, ne20, mg25, si28, ni56

! declare the pass
      double precision tstep,tin,din,ein,xin(1),tout,eout,xout(1)

! local variables
      integer          i,nok,nbad
      double precision abar,zbar,wbar,ye,xcess,snuda,snudz,asum,conv
      parameter        (conv = ev2erg*1.0d6*avo)

! for the integration driver
      external         aprox13,daprox13,stifbs_gift
      integer          kount
      double precision beg,stptry,stpmin,tend,ys2(abignet*nzmax),&
                       odescal,tol
      parameter        (tol = 1.0d-6, odescal = 1.0d-6)



! load the passed temperature and density in common block
      btemp = tin
      bden  = din 


! load the mass fractions
      xmass(ionbeg:ionend) = xin(ionbeg:ionend)


! get the molar fractions
      call azbar_burn(xmass(ionbeg),aion(ionbeg),zion(ionbeg),wion(ionbeg),ionmax, &
                 ymass(ionbeg),abar,zbar,wbar,ye,xcess)


! stuff the initial conditions into ys2
      ys2(ionbeg:ionend) = ymass(ionbeg:ionend)
      ys2(iener) = ein
      ys2(itemp) = tin


! set single step (tend=tstep) or hydrostatic ending times.
! the variable tstep has two meanings here. tstep in single step mode
! is the size of the time step to try. tstep in hydrostatic 
! mode is the ending integration time. the integration driver really
! gets some exercise if tstep is large in single step mode.

      beg  = 0.0d0
      tend = tstep
      if (one_step) then
       stptry = tstep
       stpmin = tstep * 1.0d-20
      else
       stptry = 1.0d-16
       stpmin = stptry * 1.0d-12
      end if



! integrate the aprox13 network
      call netint(beg,stptry,stpmin,tend,ys2, &
                  tol,neqs,nok,nbad,kount,odescal, &
                  aprox13,daprox13,stifbs_gift)

!      write(6,*) nbad,nok


! set the quantities to return
      xout(ionbeg:ionend) = ys2(ionbeg:ionend) * aion(ionbeg:ionend)
      tout = ys2(itemp)
      eout = ys2(iener)


! for hydro applications, renormalize the ouput composition
      asum = 1.0d0/sum(xout(ionbeg:ionend))
      xout(ionbeg:ionend) = max(xout(ionbeg:ionend) * asum, 1.0d-30)

      return
      end




!---------------------------------------------------------------------
! this file contains aprox13 network

! routine aprox13 sets up the odes
! routine rhs evaluates the right hand sides
! routine daprox13 sets up the dense aprox13 jacobian
! routine aprox13rat generates the reaction rates for routine aprox13
! routine screen_aprox13 applies screening corrections to the raw rates
! routine init_burner initializes the aprox13 network




      subroutine aprox13(tt,y,dydt)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'vector_eos.dek'

! this routine sets up the system of ode's for the aprox13 nuclear reactions.
! this is an alpha chain + heavy ion network with (a,p)(p,g) links
!
! isotopes: he4,  c12,  o16,  ne20, mg24, si28, s32,
!           ar36, ca40, ti44, cr48, fe52, ni56



! declare the pass
      double precision tt,y(1),dydt(1)

! local variables
      logical          deriva
      parameter        (deriva = .false.)
      integer          i
      double precision abar,zbar,snuda,snudz


! positive definite mass fractions
      do i=ionbeg,ionend
       y(i) = min(1.0d0,max(y(i),1.0d-30))
      enddo


! generate abar and zbar for this composition
      abar = 1.0d0/sum(y(ionbeg:ionend))
      zbar = sum(zion(ionbeg:ionend)*y(ionbeg:ionend)) * abar


! positive definite temperatures 
       y(itemp) = min(1.0d11,max(y(itemp),1.0d4))


! set the common block temperature and density
       btemp = y(itemp)


! get the reaction rates and screen them
      call aprox13rat
      call screen_aprox13(y)


! get the right hand side of the odes

      call rhs(y,ratdum,dydt,deriva)


! append the energy equation
      call ener_gener_rate(dydt,sdot)
      call sneut5(btemp,bden,abar,zbar, &
                  sneut,dsneutdt,dsneutdd,snuda,snudz)
      dydt(iener) = sdot - sneut



! hydrostatic cases
      if (one_step) then
       dydt(itemp) = 0.0d0


! append the temperature equation
      else if (self_heat_const_den) then
       temp_row(1) = y(itemp) ; den_row(1) = bden ; abar_row(1) = abar ; zbar_row(1) = zbar
       jlo_eos = 1 ; jhi_eos = 1
       call helmeos
       dydt(itemp) = dydt(iener)/cv_row(1)
      end if

      return
      end





      subroutine rhs(y,rate,dydt,deriva)
      include 'implno.dek'
      include 'network.dek'

! evaluates the right hand side of the aprox13 odes

! declare the pass
      logical          deriva
      double precision y(1),rate(1),dydt(1)


! local variables
      integer          i
      double precision sixth
      parameter        (sixth = 1.0d0/6.0d0)


! zero the abundance odes
      do i=ionbeg,ionend
       dydt(i) = 0.0d0
      enddo


! set up the system of odes:
! he4 reactions
! heavy ion reactions
      dydt(ihe4) = &
             0.5d0 * y(ic12) * y(ic12) * rate(ir1212) &
           + 0.5d0 * y(ic12) * y(io16) * rate(ir1216) &
           + 0.56d0* 0.5d0 * y(io16)*y(io16) * rate(ir1616)

! (a,g) and (g,a) reactions
      dydt(ihe4) =  dydt(ihe4) &
           - 0.5d0 * y(ihe4)*y(ihe4)*y(ihe4)*rate(ir3a) &
           + 3.0d0 * y(ic12) * rate(irg3a) &
           - y(ihe4)  * y(ic12) * rate(ircag) &
           + y(io16)  * rate(iroga) &
           - y(ihe4)  * y(io16) * rate(iroag) &
           + y(ine20) * rate(irnega) &
           - y(ihe4)  * y(ine20) * rate(irneag) &
           + y(img24) * rate(irmgga) &
           - y(ihe4)  * y(img24)* rate(irmgag) &
           + y(isi28) * rate(irsiga) &
           - y(ihe4)  * y(isi28)*rate(irsiag) &
           + y(is32)  * rate(irsga)

      dydt(ihe4) =  dydt(ihe4) &
           - y(ihe4)  * y(is32) * rate(irsag) &
           + y(iar36) * rate(irarga) &
           - y(ihe4)  * y(iar36)*rate(irarag) &
           + y(ica40) * rate(ircaga) &
           - y(ihe4)  * y(ica40)*rate(ircaag) &
           + y(iti44) * rate(irtiga) &
           - y(ihe4)  * y(iti44)*rate(irtiag) &
           + y(icr48) * rate(ircrga) &
           - y(ihe4)  * y(icr48)*rate(ircrag) &
           + y(ife52) * rate(irfega) &
           - y(ihe4)  * y(ife52) * rate(irfeag) &
           + y(ini56) * rate(irniga)


! (a,p)(p,g) and (g,p)(p,a) reactions

      if (.not.deriva) then
       dydt(ihe4) =  dydt(ihe4) &
           + 0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616) &
           - y(ihe4)  * y(img24) * rate(irmgap)*(1.0d0-rate(irr1)) &
           + y(isi28) * rate(irsigp) * rate(irr1) &
           - y(ihe4)  * y(isi28) * rate(irsiap)*(1.0d0-rate(irs1)) &
           + y(is32)  * rate(irsgp) * rate(irs1) &
           - y(ihe4)  * y(is32) * rate(irsap)*(1.0d0-rate(irt1)) &
           + y(iar36) * rate(irargp) * rate(irt1) &
           - y(ihe4)  * y(iar36) * rate(irarap)*(1.0d0-rate(iru1)) &
           + y(ica40) * rate(ircagp) * rate(iru1) &
           - y(ihe4)  * y(ica40) * rate(ircaap)*(1.0d0-rate(irv1)) &
           + y(iti44) * rate(irtigp) * rate(irv1)

       dydt(ihe4) =  dydt(ihe4) &
           - y(ihe4)  * y(iti44) * rate(irtiap)*(1.0d0-rate(irw1)) &
           + y(icr48) * rate(ircrgp) * rate(irw1) &
           - y(ihe4)  * y(icr48) * rate(ircrap)*(1.0d0-rate(irx1)) &
           + y(ife52) * rate(irfegp) * rate(irx1) &
           - y(ihe4)  * y(ife52) * rate(irfeap)*(1.0d0-rate(iry1)) &
           + y(ini56) * rate(irnigp) * rate(iry1)

      else
       dydt(ihe4) =  dydt(ihe4) &
           + 0.34d0*0.5d0*y(io16)*y(io16)* &
             (ratdum(irs1)*rate(ir1616) + rate(irs1)*ratdum(ir1616)) &
           - y(ihe4)*y(img24)*(rate(irmgap)*(1.0d0 - ratdum(irr1)) &
                                     - ratdum(irmgap)*rate(irr1)) &
           + y(isi28) * (ratdum(irsigp) * rate(irr1) + &
                            rate(irsigp) * ratdum(irr1)) &
           - y(ihe4)*y(isi28)*(rate(irsiap)*(1.0d0 - ratdum(irs1)) &
                                     - ratdum(irsiap)*rate(irs1)) &
           + y(is32)  * (ratdum(irsgp) * rate(irs1) + &
                            rate(irsgp) * ratdum(irs1))

       dydt(ihe4) =  dydt(ihe4) &
           - y(ihe4)*y(is32)*(rate(irsap)*(1.0d0 - ratdum(irt1)) &
                                    - ratdum(irsap)*rate(irt1)) &
           + y(iar36) * (ratdum(irargp) * rate(irt1) + &
                            rate(irargp) * ratdum(irt1)) &
           - y(ihe4)*y(iar36)*(rate(irarap)*(1.0d0 - ratdum(iru1)) &
                                     - ratdum(irarap)*rate(iru1)) &
           + y(ica40) * (ratdum(ircagp) * rate(iru1) + &
                            rate(ircagp) * ratdum(iru1)) &
           - y(ihe4)*y(ica40)*(rate(ircaap)*(1.0d0-ratdum (irv1)) &
                                     - ratdum(ircaap)*rate(irv1)) &
           + y(iti44) * (ratdum(irtigp) * rate(irv1) + &
                            rate(irtigp) * ratdum(irv1))

       dydt(ihe4) =  dydt(ihe4) &
           - y(ihe4)*y(iti44)*(rate(irtiap)*(1.0d0 - ratdum(irw1)) &
                                     - ratdum(irtiap)*rate(irw1)) &
           + y(icr48) * (ratdum(ircrgp) * rate(irw1) + &
                            rate(ircrgp) * ratdum(irw1)) &
           - y(ihe4)*y(icr48)*(rate(ircrap)*(1.0d0 - ratdum(irx1)) &
                                     - ratdum(ircrap)*rate(irx1)) &
           + y(ife52) * (ratdum(irfegp) * rate(irx1) + &
                            rate(irfegp) * ratdum(irx1)) &
           - y(ihe4)*y(ife52)*(rate(irfeap)*(1.0d0 - ratdum(iry1)) &
                                    - ratdum(irfeap)*rate(iry1)) &
           + y(ini56) * (ratdum(irnigp) * rate(iry1) + &
                            rate(irnigp) * ratdum(iry1))
      end if




! c12 reactions
      dydt(ic12) = -y(ic12) * y(ic12) * rate(ir1212) &
           - y(ic12) * y(io16) * rate(ir1216) &
           + sixth * y(ihe4)*y(ihe4)*y(ihe4)*rate(ir3a) &
           - y(ic12) * rate(irg3a) &
           - y(ic12) * y(ihe4) * rate(ircag) &
           + y(io16) * rate(iroga)

! o16 reactions
      dydt(io16) = -y(ic12) * y(io16) * rate(ir1216) &
           - y(io16) * y(io16) * rate(ir1616) &
           + y(ic12) * y(ihe4) * rate(ircag) &
           - y(io16) * y(ihe4) * rate(iroag) &
           - y(io16) * rate(iroga) &
           + y(ine20) * rate(irnega)

! ne20 reactions
      dydt(ine20) =  0.5d0 * y(ic12) * y(ic12) * rate(ir1212) &
           + y(io16) * y(ihe4) * rate(iroag) &
           - y(ine20) * y(ihe4) * rate(irneag) &
           - y(ine20) * rate(irnega) &
           + y(img24) * rate(irmgga)


! mg24 reactions
      dydt(img24)  = 0.5d0 * y(ic12) * y(io16) * rate(ir1216) &
           + y(ine20) * y(ihe4) * rate(irneag) &
           - y(img24) * y(ihe4) * rate(irmgag) &
           - y(img24) * rate(irmgga) &
           + y(isi28) * rate(irsiga)

      if (.not.deriva) then
       dydt(img24)  = dydt(img24) &
           - y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1)) &
           + y(isi28) * rate(irr1) * rate(irsigp)

      else
       dydt(img24)  = dydt(img24) &
           - y(img24)*y(ihe4)*(rate(irmgap)*(1.0d0 - ratdum(irr1)) &
                                     - ratdum(irmgap)*rate(irr1)) &
           + y(isi28) * (ratdum(irr1) * rate(irsigp) + &
                            rate(irr1) * ratdum(irsigp))
      end if



! si28 reactions
      dydt(isi28)  =  0.5d0 * y(ic12) * y(io16) * rate(ir1216) &
           + 0.56d0*0.5d0*y(io16)*y(io16)*rate(ir1616) &
           + y(img24) * y(ihe4) * rate(irmgag) &
           - y(isi28) * y(ihe4) * rate(irsiag) &
           - y(isi28) * rate(irsiga) &
           + y(is32)  * rate(irsga)

      if (.not.deriva) then
       dydt(isi28)  = dydt(isi28) &
           + 0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616) &
           + y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1)) &
           - y(isi28) * rate(irr1) * rate(irsigp) &
           - y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1)) &
           + y(is32)  * rate(irs1) * rate(irsgp)

      else
       dydt(isi28)  = dydt(isi28) &
           + 0.34d0*0.5d0*y(io16)*y(io16)* &
              (ratdum(irs1)*rate(ir1616) + rate(irs1)*ratdum(ir1616)) &
           + y(img24)*y(ihe4)*(rate(irmgap)*(1.0d0 - ratdum(irr1)) &
                                     - ratdum(irmgap)*rate(irr1)) &
           - y(isi28)*(ratdum(irr1) * rate(irsigp) + &
                          rate(irr1) * ratdum(irsigp)) &
           - y(isi28)*y(ihe4)*(rate(irsiap)*(1.0d0 - ratdum(irs1)) &
                                     - ratdum(irsiap)*rate(irs1)) &
           + y(is32)*(ratdum(irs1) * rate(irsgp) + &
                          rate(irs1) * ratdum(irsgp))
      end if



! s32 reactions
      dydt(is32) = 0.1d0 * 0.5d0 *y(io16)*y(io16)*rate(ir1616) &
            +  y(isi28) * y(ihe4) * rate(irsiag) &
            - y(is32) * y(ihe4) * rate(irsag) &
            - y(is32) * rate(irsga) &
            + y(iar36) * rate(irarga)

      if (.not.deriva) then
       dydt(is32)  = dydt(is32) &
            + 0.34d0*0.5d0*y(io16)**2*rate(ir1616)*(1.0d0-rate(irs1)) &
            + y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1)) &
            - y(is32) * rate(irs1) * rate(irsgp) &
            - y(is32) * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1)) &
            + y(iar36) * rate(irt1) * rate(irargp)
      else
       dydt(is32)  = dydt(is32) &
            + 0.34d0*0.5d0*y(io16)**2 * &
           (rate(ir1616)*(1.0d0-ratdum(irs1))-ratdum(ir1616)*rate(irs1)) &
            + y(isi28)*y(ihe4)*(rate(irsiap)*(1.0d0-ratdum(irs1)) &
                                     - ratdum(irsiap)*rate(irs1)) &
            - y(is32)*(ratdum(irs1) * rate(irsgp) + &
                          rate(irs1) * ratdum(irsgp)) &
            - y(is32)*y(ihe4)*(rate(irsap)*(1.0d0-ratdum(irt1)) &
                                     - ratdum(irsap)*rate(irt1)) &
            + y(iar36)*(ratdum(irt1) * rate(irargp) + &
                           rate(irt1) * ratdum(irargp))
      end if


! ar36 reactions
      dydt(iar36) =  y(is32)  * y(ihe4) * rate(irsag) &
            - y(iar36) * y(ihe4) * rate(irarag) &
            - y(iar36) * rate(irarga) &
            + y(ica40) * rate(ircaga)

      if (.not.deriva) then
       dydt(iar36)  = dydt(iar36) &
            + y(is32)  * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1)) &
            - y(iar36) * rate(irt1) * rate(irargp) &
            - y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1)) &
            + y(ica40) * rate(ircagp) * rate(iru1)

      else
       dydt(iar36)  = dydt(iar36) &
            + y(is32)*y(ihe4)*(rate(irsap)*(1.0d0 - ratdum(irt1)) &
                                     - ratdum(irsap)*rate(irt1)) &
            - y(iar36)*(ratdum(irt1) * rate(irargp) + &
                           rate(irt1) * ratdum(irargp)) &
            - y(iar36)*y(ihe4)*(rate(irarap)*(1.0d0-ratdum(iru1)) &
                                      - ratdum(irarap)*rate(iru1)) &
            + y(ica40)*(ratdum(ircagp) * rate(iru1) + &
                           rate(ircagp) * ratdum(iru1))
      end if


! ca40 reactions
      dydt(ica40) =  y(iar36) * y(ihe4) * rate(irarag) &
            - y(ica40) * y(ihe4) * rate(ircaag) &
            - y(ica40) * rate(ircaga) &
            + y(iti44) * rate(irtiga)


      if (.not.deriva) then
       dydt(ica40)  = dydt(ica40) &
            + y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1)) &
            - y(ica40) * rate(ircagp) * rate(iru1) &
            - y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1)) &
            + y(iti44) * rate(irtigp) * rate(irv1)

      else
       dydt(ica40)  = dydt(ica40) &
            + y(iar36)*y(ihe4)*(rate(irarap)*(1.0d0-ratdum(iru1)) &
                                      - ratdum(irarap)*rate(iru1)) &
            - y(ica40)*(ratdum(ircagp) * rate(iru1) + &
                           rate(ircagp) * ratdum(iru1)) &
            - y(ica40)*y(ihe4)*(rate(ircaap)*(1.0d0-ratdum(irv1)) &
                                      - ratdum(ircaap)*rate(irv1)) &
            + y(iti44)*(ratdum(irtigp) * rate(irv1) + &
                           rate(irtigp) * ratdum(irv1))
      end if



! ti44 reactions
      dydt(iti44) =  y(ica40) * y(ihe4) * rate(ircaag) &
            - y(iti44) * y(ihe4) * rate(irtiag) &
            - y(iti44) * rate(irtiga) &
            + y(icr48) * rate(ircrga)


      if (.not.deriva) then
       dydt(iti44)  = dydt(iti44) &
            + y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1)) &
            - y(iti44) * rate(irv1) * rate(irtigp) &
            - y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1)) &
            + y(icr48) * rate(irw1) * rate(ircrgp)
      else
       dydt(iti44)  = dydt(iti44) &
            + y(ica40)*y(ihe4)*(rate(ircaap)*(1.0d0-ratdum(irv1)) &
                                      - ratdum(ircaap)*rate(irv1)) &
            - y(iti44)*(ratdum(irv1) * rate(irtigp) + &
                           rate(irv1) * ratdum(irtigp)) &
            - y(iti44)*y(ihe4)*(rate(irtiap)*(1.0d0-ratdum(irw1)) &
                                      - ratdum(irtiap)*rate(irw1)) &
            + y(icr48)*(ratdum(irw1) * rate(ircrgp) + &
                           rate(irw1) * ratdum(ircrgp))
      end if


! cr48 reactions
      dydt(icr48) =  y(iti44) * y(ihe4) * rate(irtiag) &
            - y(icr48) * y(ihe4) * rate(ircrag) &
            - y(icr48) * rate(ircrga) &
            + y(ife52) * rate(irfega)

      if (.not.deriva) then
       dydt(icr48)  = dydt(icr48) &
            + y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1)) &
            - y(icr48) * rate(irw1) * rate(ircrgp) &
            - y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) &
            + y(ife52) * rate(irx1) * rate(irfegp)

      else
       dydt(icr48)  = dydt(icr48) &
            + y(iti44)*y(ihe4)*(rate(irtiap)*(1.0d0-ratdum(irw1)) &
                                      - ratdum(irtiap)*rate(irw1)) &
            - y(icr48)*(ratdum(irw1) * rate(ircrgp) + &
                           rate(irw1) * ratdum(ircrgp)) &
            - y(icr48)*y(ihe4)*(rate(ircrap)*(1.0d0-ratdum(irx1)) &
                                      - ratdum(ircrap)*rate(irx1)) &
            + y(ife52)*(ratdum(irx1) * rate(irfegp) + &
                           rate(irx1) * ratdum(irfegp))
      end if


! fe52 reactions
      dydt(ife52) =  y(icr48) * y(ihe4) * rate(ircrag) &
            - y(ife52) * y(ihe4) * rate(irfeag) &
            - y(ife52) * rate(irfega) &
            + y(ini56) * rate(irniga)

      if (.not.deriva) then
       dydt(ife52)  = dydt(ife52) &
            + y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) &
            - y(ife52) * rate(irx1) * rate(irfegp) &
            - y(ife52) * y(ihe4) * rate(irfeap)*(1.0d0-rate(iry1)) &
            + y(ini56) * rate(iry1) * rate(irnigp)

      else
       dydt(ife52)  = dydt(ife52) &
            + y(icr48)*y(ihe4)*(rate(ircrap)*(1.0d0-ratdum(irx1)) &
                                      - ratdum(ircrap)*rate(irx1)) &
            - y(ife52)*(ratdum(irx1) * rate(irfegp) + &
                           rate(irx1) * ratdum(irfegp)) &
            - y(ife52)*y(ihe4)*(rate(irfeap)*(1.0d0-ratdum(iry1)) &
                                      - ratdum(irfeap)*rate(iry1)) &
            + y(ini56)*(ratdum(iry1) * rate(irnigp) + &
                           rate(iry1) * ratdum(irnigp))
      end if


! ni56 reactions
      dydt(ini56) =  y(ife52) * y(ihe4) * rate(irfeag) &
            - y(ini56) * rate(irniga)

      if (.not.deriva) then
       dydt(ini56)  = dydt(ini56) &
            + y(ife52) * y(ihe4) * rate(irfeap)*(1.0d0-rate(iry1)) &
            - y(ini56) * rate(iry1) * rate(irnigp)

      else
       dydt(ini56)  = dydt(ini56) &
            + y(ife52)*y(ihe4)*(rate(irfeap)*(1.0d0-ratdum(iry1)) &
                                      - ratdum(irfeap)*rate(iry1)) &
            - y(ini56)*(ratdum(iry1) * rate(irnigp) + &
                           rate(iry1) * ratdum(irnigp))
      end if


      return
      end




      subroutine daprox13(tt,y,dfdy,nlog,nphys)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'vector_eos.dek'

! this routine sets up the dense aprox13 jacobian


! declare the pass
      integer          nlog,nphys
      double precision tt,y(1),dfdy(nphys,nphys)


! local variables
      logical          deriva
      parameter        (deriva = .true.)
      integer          i,j
      double precision abar,zbar,snuda,snudz,zz,enuc_conv2
      parameter        (enuc_conv2 = -avo*clight*clight)



! zero the jacobian
      dfdy(1:nlog,1:nlog) = 0.0d0


! positive definite mass fractions
      do i=ionbeg,ionend
       y(i) = min(1.0d0,max(y(i),1.0d-30))
      enddo


! generate abar and zbar for this composition
      abar = 1.0d0/sum(y(ionbeg:ionend))
      zbar = sum(zion(ionbeg:ionend)*y(ionbeg:ionend)) * abar


! positive definite temperatures 
       y(itemp) = min(1.0d11,max(y(itemp),1.0d4))


! set the common block temperature and density
       btemp = y(itemp)


! get the reaction rates and screen them
      call aprox13rat
      call screen_aprox13(y)



! he4 jacobian elements
      dfdy(ihe4,ihe4)  = -1.5d0*y(ihe4)*y(ihe4)*ratdum(ir3a) &
                  - y(ic12)  * ratdum(ircag) &
                  - y(io16)  * ratdum(iroag) &
                  - y(ine20) * ratdum(irneag) &
                  - y(img24) * ratdum(irmgag) &
                  - y(isi28) * ratdum(irsiag) &
                  - y(is32)  * ratdum(irsag) &
                  - y(iar36) * ratdum(irarag) &
                  - y(ica40) * ratdum(ircaag) &
                  - y(iti44) * ratdum(irtiag) &
                  - y(icr48) * ratdum(ircrag) &
                  - y(ife52) * ratdum(irfeag)

      dfdy(ihe4,ihe4)  = dfdy(ihe4,ihe4) &
                  - y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1)) &
                  - y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) &
                  - y(is32) * ratdum(irsap)   * (1.0d0-ratdum(irt1)) &
                  - y(iar36) * ratdum(irarap) * (1.0d0-ratdum(iru1)) &
                  - y(ica40) * ratdum(ircaap) * (1.0d0-ratdum(irv1)) &
                  - y(iti44) * ratdum(irtiap) * (1.0d0-ratdum(irw1)) &
                  - y(icr48) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) &
                  - y(ife52) * ratdum(irfeap) * (1.0d0-ratdum(iry1))


      dfdy(ihe4,ic12)  = y(ic12) * ratdum(ir1212) &
                  + 0.5d0 * y(io16) * ratdum(ir1216) &
                  + 3.0d0 * ratdum(irg3a) &
                  - y(ihe4) * ratdum(ircag)

      dfdy(ihe4,io16)  = 0.5d0 * y(ic12) * ratdum(ir1216) &
                  + 1.12d0 * 0.5d0 * y(io16) * ratdum(ir1616) &
                  + 0.68d0*ratdum(irs1)*0.5d0*y(io16)*ratdum(ir1616) &
                  + ratdum(iroga) &
                  - y(ihe4) * ratdum(iroag)

      dfdy(ihe4,ine20) =  ratdum(irnega) &
                  - y(ihe4) * ratdum(irneag)


      dfdy(ihe4,img24) =   ratdum(irmgga) &
                  - y(ihe4) * ratdum(irmgag) &
                  - y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

      dfdy(ihe4,isi28) =   ratdum(irsiga) &
                  - y(ihe4) * ratdum(irsiag) &
                  - y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) &
                  + ratdum(irr1) * ratdum(irsigp)


      dfdy(ihe4,is32)  =   ratdum(irsga) &
                  - y(ihe4) * ratdum(irsag) &
                  - y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1)) &
                  + ratdum(irs1) * ratdum(irsgp)

      dfdy(ihe4,iar36) =   ratdum(irarga) &
                  - y(ihe4) * ratdum(irarag) &
                  - y(ihe4) * ratdum(irarap) * (1.0d0-ratdum(iru1)) &
                  + ratdum(irt1) * ratdum(irargp)

      dfdy(ihe4,ica40) =   ratdum(ircaga) &
                  - y(ihe4) * ratdum(ircaag) &
                  - y(ihe4) * ratdum(ircaap) * (1.0d0-ratdum(irv1)) &
                  + ratdum(iru1) * ratdum(ircagp)

      dfdy(ihe4,iti44) =   ratdum(irtiga) &
                  - y(ihe4) * ratdum(irtiag) &
                  - y(ihe4) * ratdum(irtiap) * (1.0d0-ratdum(irw1)) &
                  + ratdum(irv1) * ratdum(irtigp)

      dfdy(ihe4,icr48) =   ratdum(ircrga) &
                  - y(ihe4) * ratdum(ircrag) &
                  - y(ihe4) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) &
                  + ratdum(irw1) * ratdum(ircrgp)

      dfdy(ihe4,ife52) =   ratdum(irfega) &
                  - y(ihe4) * ratdum(irfeag) &
                  - y(ihe4) * ratdum(irfeap) * (1.0d0-ratdum(iry1)) &
                  + ratdum(irx1) * ratdum(irfegp)

      dfdy(ihe4,ini56) =   ratdum(irniga) &
                  + ratdum(iry1) * ratdum(irnigp)


! c12 jacobian elements
      dfdy(ic12,ihe4) = 0.5d0*y(ihe4)*y(ihe4)*ratdum(ir3a) &
                - y(ic12) * ratdum(ircag)

      dfdy(ic12,ic12) = -2.0d0 * y(ic12) * ratdum(ir1212) &
                 - y(io16) * ratdum(ir1216) &
                 - ratdum(irg3a) &
                 - y(ihe4) * ratdum(ircag)

      dfdy(ic12,io16) = -y(ic12) * ratdum(ir1216) &
                + ratdum(iroga)



! o16 jacobian elements
      dfdy(io16,ihe4) = y(ic12)*ratdum(ircag) &
               - y(io16)*ratdum(iroag)

      dfdy(io16,ic12) = -y(io16)*ratdum(ir1216) &
                + y(ihe4)*ratdum(ircag)

      dfdy(io16,io16) = - y(ic12) * ratdum(ir1216) &
                 - 2.0d0 * y(io16) * ratdum(ir1616) &
                 - y(ihe4) * ratdum(iroag) &
                 - ratdum(iroga)

      dfdy(io16,ine20) = ratdum(irnega)



! ne20 jacobian elements
      dfdy(ine20,ihe4)  = y(io16) * ratdum(iroag) &
                 - y(ine20) * ratdum(irneag)

      dfdy(ine20,ic12)  = y(ic12) * ratdum(ir1212)

      dfdy(ine20,io16)  = y(ihe4) * ratdum(iroag)

      dfdy(ine20,ine20) = -y(ihe4) * ratdum(irneag) &
                  - ratdum(irnega)

      dfdy(ine20,img24) = ratdum(irmgga)



! mg24 jacobian elements
      dfdy(img24,ihe4)  = y(ine20) * ratdum(irneag) &
                  -y(img24) * ratdum(irmgag) &
                  -y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

      dfdy(img24,ic12)  = 0.5d0 * y(io16) * ratdum(ir1216)

      dfdy(img24,io16)  = 0.5d0 * y(ic12) * ratdum(ir1216)

      dfdy(img24,ine20) = y(ihe4) * ratdum(irneag)

      dfdy(img24,img24) = -y(ihe4) * ratdum(irmgag) &
                  - ratdum(irmgga) &
                  - y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

      dfdy(img24,isi28) = ratdum(irsiga) &
                 + ratdum(irr1) * ratdum(irsigp)



! si28 jacobian elements
      dfdy(isi28,ihe4)  = y(img24) * ratdum(irmgag) &
                 - y(isi28) * ratdum(irsiag) &
                 + y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1)) &
                 - y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1))

      dfdy(isi28,ic12)  = 0.5d0 * y(io16) * ratdum(ir1216)

      dfdy(isi28,io16)  =   0.5d0 * y(ic12) * ratdum(ir1216) &
                   + 1.12d0 * 0.5d0*y(io16) * ratdum(ir1616) &
                   + 0.68d0*0.5d0*y(io16)*ratdum(irs1)*ratdum(ir1616)

      dfdy(isi28,img24) = y(ihe4) * ratdum(irmgag) &
                 + y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

      dfdy(isi28,isi28) = -y(ihe4) * ratdum(irsiag) &
                  - ratdum(irsiga) &
                  - ratdum(irr1) * ratdum(irsigp) &
                  - y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))

      dfdy(isi28,is32)  = ratdum(irsga) &
                 + ratdum(irs1) * ratdum(irsgp)



! s32 jacobian elements
      dfdy(is32,ihe4)  = y(isi28) * ratdum(irsiag) &
                - y(is32) * ratdum(irsag) &
                + y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) &
                - y(is32) * ratdum(irsap) * (1.0d0-ratdum(irt1))

      dfdy(is32,io16)  = 0.68d0*0.5d0*y(io16) &
                        *ratdum(ir1616)*(1.0d0-ratdum(irs1)) &
                  + 0.2d0 * 0.5d0*y(io16) * ratdum(ir1616)

      dfdy(is32,isi28) = y(ihe4) * ratdum(irsiag) &
                 + y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))

      dfdy(is32,is32)  = -y(ihe4) * ratdum(irsag) &
                 - ratdum(irsga) &
                 - ratdum(irs1) * ratdum(irsgp) &
                 - y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))

      dfdy(is32,iar36) = ratdum(irarga) &
                + ratdum(irt1) * ratdum(irargp)



! ar36 jacobian elements
      dfdy(iar36,ihe4)  = y(is32)  * ratdum(irsag) &
                - y(iar36) * ratdum(irarag) &
                + y(is32)  * ratdum(irsap) * (1.0d0-ratdum(irt1)) &
                - y(iar36) * ratdum(irarap) * (1.0d0-ratdum(iru1))

      dfdy(iar36,is32)  = y(ihe4) * ratdum(irsag) &
                  + y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))

      dfdy(iar36,iar36) = -y(ihe4) * ratdum(irarag) &
                  - ratdum(irarga) &
                  - ratdum(irt1) * ratdum(irargp) &
                  - y(ihe4) * ratdum(irarap) * (1.0d0-ratdum(iru1))

      dfdy(iar36,ica40) = ratdum(ircaga) &
                        + ratdum(ircagp) * ratdum(iru1)



! ca40 jacobian elements
      dfdy(ica40,ihe4)   = y(iar36) * ratdum(irarag) &
                  - y(ica40) * ratdum(ircaag) &
                  + y(iar36) * ratdum(irarap)*(1.0d0-ratdum(iru1)) &
                  - y(ica40) * ratdum(ircaap)*(1.0d0-ratdum(irv1))

      dfdy(ica40,iar36)  = y(ihe4) * ratdum(irarag) &
                  + y(ihe4) * ratdum(irarap)*(1.0d0-ratdum(iru1))

      dfdy(ica40,ica40)  = -y(ihe4) * ratdum(ircaag) &
                   - ratdum(ircaga) &
                   - ratdum(ircagp) * ratdum(iru1) &
                   - y(ihe4) * ratdum(ircaap)*(1.0d0-ratdum(irv1))

      dfdy(ica40,iti44)  = ratdum(irtiga) &
                          + ratdum(irtigp) * ratdum(irv1)



! ti44 jacobian elements
      dfdy(iti44,ihe4)   = y(ica40) * ratdum(ircaag) &
                  - y(iti44) * ratdum(irtiag) &
                  + y(ica40) * ratdum(ircaap)*(1.0d0-ratdum(irv1)) &
                  - y(iti44) * ratdum(irtiap)*(1.0d0-ratdum(irw1))

      dfdy(iti44,ica40)  = y(ihe4) * ratdum(ircaag) &
                  + y(ihe4) * ratdum(ircaap)*(1.0d0-ratdum(irv1))

      dfdy(iti44,iti44)  = -y(ihe4) * ratdum(irtiag) &
                   - ratdum(irtiga) &
                   - ratdum(irv1) * ratdum(irtigp) &
                   - y(ihe4) * ratdum(irtiap)*(1.0d0-ratdum(irw1))

      dfdy(iti44,icr48)  = ratdum(ircrga) &
                  + ratdum(irw1) * ratdum(ircrgp)



! cr48 jacobian elements
      dfdy(icr48,ihe4)  = y(iti44) * ratdum(irtiag) &
                 - y(icr48) * ratdum(ircrag) &
                 + y(iti44) * ratdum(irtiap)*(1.0d0-ratdum(irw1)) &
                 - y(icr48) * ratdum(ircrap)*(1.0d0-ratdum(irx1))

      dfdy(icr48,iti44) = y(ihe4) * ratdum(irtiag) &
                 + y(ihe4) * ratdum(irtiap)*(1.0d0-ratdum(irw1))

      dfdy(icr48,icr48) = -y(ihe4) * ratdum(ircrag) &
                  - ratdum(ircrga) &
                  - ratdum(irw1) * ratdum(ircrgp) &
                  - y(ihe4) * ratdum(ircrap)*(1.0d0-ratdum(irx1))

      dfdy(icr48,ife52) = ratdum(irfega) &
                 + ratdum(irx1) * ratdum(irfegp)



! fe52 jacobian elements
      dfdy(ife52,ihe4)  = y(icr48) * ratdum(ircrag) &
                 - y(ife52) * ratdum(irfeag) &
                 + y(icr48) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) &
                 - y(ife52) * ratdum(irfeap) * (1.0d0-ratdum(iry1))

      dfdy(ife52,icr48) = y(ihe4) * ratdum(ircrag) &
                 + y(ihe4) * ratdum(ircrap) * (1.0d0-ratdum(irx1))

      dfdy(ife52,ife52) = - y(ihe4) * ratdum(irfeag) &
                   - ratdum(irfega) &
                   - ratdum(irx1) * ratdum(irfegp) &
                   - y(ihe4) * ratdum(irfeap) * (1.0d0-ratdum(iry1))

      dfdy(ife52,ini56) = ratdum(irniga) &
                 + ratdum(iry1) * ratdum(irnigp)



! ni56 jacobian elements
      dfdy(ini56,ihe4)  = y(ife52) * ratdum(irfeag) &
                 + y(ife52) * ratdum(irfeap) * (1.0d0-ratdum(iry1))

      dfdy(ini56,ife52) = y(ihe4) * ratdum(irfeag) &
                 + y(ihe4) * ratdum(irfeap) * (1.0d0-ratdum(iry1))

      dfdy(ini56,ini56) = -ratdum(irniga) &
                  - ratdum(iry1) * ratdum(irnigp)



! append the temperature derivatives of the rate equations

      call rhs(y,dratdumdt,zwork1,deriva)
      do i=ionbeg,ionend
       dfdy(i,itemp) = zwork1(i)
      enddo



! append the energy jacobian elements

      do j=ionbeg,ionend
       do i=ionbeg,ionend
        dfdy(iener,j) = dfdy(iener,j) + dfdy(i,j)*mion(i)
       enddo
       dfdy(iener,j)     = dfdy(iener,j)*enuc_conv2
       dfdy(iener,itemp) = dfdy(iener,itemp) &
                                 + dfdy(j,itemp)*mion(j)
      enddo
      dfdy(iener,itemp) = dfdy(iener,itemp) * enuc_conv2
      dsdotdt = dfdy(iener,itemp)


! account for the neutrino losses

      call sneut5(btemp,bden,abar,zbar, &
                  sneut,dsneutdt,dsneutdd,snuda,snudz)

      do j=ionbeg,ionend
       dfdy(iener,j) = dfdy(iener,j) &
                     - ((aion(j) - abar)*abar*snuda + (zion(j) - zbar)*abar*snudz)
      enddo
      dfdy(iener,itemp) = dfdy(iener,itemp) - dsneutdt


! append the temperature jacobian elements

      if (self_heat_const_den) then
       temp_row(1) = y(itemp) ; den_row(1) = bden ; abar_row(1) = abar ; zbar_row(1) = zbar
       jlo_eos = 1 ; jhi_eos = 1
       call helmeos

! d(itemp)/d(yi)
       zz = 1.0d0/cv_row(1)
       do j=ionbeg,ionend
        dfdy(itemp,j) = zz*dfdy(iener,j)
       enddo

! d(itemp)/d(temp)
       dfdy(itemp,itemp) = zz*dfdy(iener,itemp)
      end if



      return
      end



      subroutine aprox13rat
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine generates unscreened
! nuclear reaction rates for the aprox13 network.

! declare
      double precision rrate,drratedt,drratedd


! zero the rates
      ratraw(:) = 0.0d0
      dratrawdt(:) = 0.0d0
      dratrawdd(:) = 0.0d0


      if (btemp .lt. 1.0e6) return


! get the temperature factors
      call tfactors(btemp)



! c12(a,g)o16
      call rate_c12ag(btemp,bden, &
           ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
           ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))

! triple alpha to c12
      call rate_tripalf(btemp,bden, &
           ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a), &
           ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))

! c12 + c12
      call rate_c12c12(btemp,bden, &
           ratraw(ir1212),dratrawdt(ir1212),dratrawdd(ir1212), &
           rrate,drratedt,drratedd)

! c12 + o16
      call rate_c12o16(btemp,bden, &
           ratraw(ir1216),dratrawdt(ir1216),dratrawdd(ir1216), &
           rrate,drratedt,drratedd)

! o16 + o16
      call rate_o16o16(btemp,bden, &
           ratraw(ir1616),dratrawdt(ir1616),dratrawdd(ir1616), &
           rrate,drratedt,drratedd)

! o16(a,g)ne20
      call rate_o16ag(btemp,bden, &
           ratraw(iroag),dratrawdt(iroag),dratrawdd(iroag), &
           ratraw(irnega),dratrawdt(irnega),dratrawdd(irnega))

! ne20(a,g)mg24
      call rate_ne20ag(btemp,bden, &
           ratraw(irneag),dratrawdt(irneag),dratrawdd(irneag), &
           ratraw(irmgga),dratrawdt(irmgga),dratrawdd(irmgga))

! mg24(a,g)si28
      call rate_mg24ag(btemp,bden, &
           ratraw(irmgag),dratrawdt(irmgag),dratrawdd(irmgag), &
           ratraw(irsiga),dratrawdt(irsiga),dratrawdd(irsiga))

! mg24(a,p)al27
      call rate_mg24ap(btemp,bden, &
           ratraw(irmgap),dratrawdt(irmgap),dratrawdd(irmgap), &
           ratraw(iralpa),dratrawdt(iralpa),dratrawdd(iralpa))

! al27(p,g)si28
      call rate_al27pg(btemp,bden, &
           ratraw(iralpg),dratrawdt(iralpg),dratrawdd(iralpg), &
           ratraw(irsigp),dratrawdt(irsigp),dratrawdd(irsigp))

! si28(a,g)s32
      call rate_si28ag(btemp,bden, &
           ratraw(irsiag),dratrawdt(irsiag),dratrawdd(irsiag), &
           ratraw(irsga),dratrawdt(irsga),dratrawdd(irsga))

! si28(a,p)p31
      call rate_si28ap(btemp,bden, &
           ratraw(irsiap),dratrawdt(irsiap),dratrawdd(irsiap), &
           ratraw(irppa),dratrawdt(irppa),dratrawdd(irppa))

! p31(p,g)s32
      call rate_p31pg(btemp,bden, &
           ratraw(irppg),dratrawdt(irppg),dratrawdd(irppg), &
           ratraw(irsgp),dratrawdt(irsgp),dratrawdd(irsgp))

! s32(a,g)ar36
      call rate_s32ag(btemp,bden, &
           ratraw(irsag),dratrawdt(irsag),dratrawdd(irsag), &
           ratraw(irarga),dratrawdt(irarga),dratrawdd(irarga))

! s32(a,p)cl35
      call rate_s32ap(btemp,bden, &
           ratraw(irsap),dratrawdt(irsap),dratrawdd(irsap), &
           ratraw(irclpa),dratrawdt(irclpa),dratrawdd(irclpa))

! cl35(p,g)ar36
      call rate_cl35pg(btemp,bden, &
           ratraw(irclpg),dratrawdt(irclpg),dratrawdd(irclpg), &
           ratraw(irargp),dratrawdt(irargp),dratrawdd(irargp))

! ar36(a,g)ca40
      call rate_ar36ag(btemp,bden, &
           ratraw(irarag),dratrawdt(irarag),dratrawdd(irarag), &
           ratraw(ircaga),dratrawdt(ircaga),dratrawdd(ircaga))

! ar36(a,p)k39
      call rate_ar36ap(btemp,bden, &
           ratraw(irarap),dratrawdt(irarap),dratrawdd(irarap), &
           ratraw(irkpa),dratrawdt(irkpa),dratrawdd(irkpa))

! k39(p,g)ca40
      call rate_k39pg(btemp,bden, &
           ratraw(irkpg),dratrawdt(irkpg),dratrawdd(irkpg), &
           ratraw(ircagp),dratrawdt(ircagp),dratrawdd(ircagp))

! ca40(a,g)ti44
      call rate_ca40ag(btemp,bden, &
           ratraw(ircaag),dratrawdt(ircaag),dratrawdd(ircaag), &
           ratraw(irtiga),dratrawdt(irtiga),dratrawdd(irtiga))

! ca40(a,p)sc43
      call rate_ca40ap(btemp,bden, &
           ratraw(ircaap),dratrawdt(ircaap),dratrawdd(ircaap), &
           ratraw(irscpa),dratrawdt(irscpa),dratrawdd(irscpa))

! sc43(p,g)ti44
      call rate_sc43pg(btemp,bden, &
           ratraw(irscpg),dratrawdt(irscpg),dratrawdd(irscpg), &
           ratraw(irtigp),dratrawdt(irtigp),dratrawdd(irtigp))

! ti44(a,g)cr48
      call rate_ti44ag(btemp,bden, &
           ratraw(irtiag),dratrawdt(irtiag),dratrawdd(irtiag), &
           ratraw(ircrga),dratrawdt(ircrga),dratrawdd(ircrga))

! ti44(a,p)v47
      call rate_ti44ap(btemp,bden, &
           ratraw(irtiap),dratrawdt(irtiap),dratrawdd(irtiap), &
           ratraw(irvpa),dratrawdt(irvpa),dratrawdd(irvpa))

! v47(p,g)cr48
      call rate_v47pg(btemp,bden, &
           ratraw(irvpg),dratrawdt(irvpg),dratrawdd(irvpg), &
           ratraw(ircrgp),dratrawdt(ircrgp),dratrawdd(ircrgp))

! cr48(a,g)fe52
      call rate_cr48ag(btemp,bden, &
           ratraw(ircrag),dratrawdt(ircrag),dratrawdd(ircrag), &
           ratraw(irfega),dratrawdt(irfega),dratrawdd(irfega))

! cr48(a,p)mn51
      call rate_cr48ap(btemp,bden, &
           ratraw(ircrap),dratrawdt(ircrap),dratrawdd(ircrap), &
           ratraw(irmnpa),dratrawdt(irmnpa),dratrawdd(irmnpa))

! mn51(p,g)fe52
      call rate_mn51pg(btemp,bden, &
           ratraw(irmnpg),dratrawdt(irmnpg),dratrawdd(irmnpg), &
           ratraw(irfegp),dratrawdt(irfegp),dratrawdd(irfegp))

! fe52(a,g)ni56
      call rate_fe52ag(btemp,bden, &
           ratraw(irfeag),dratrawdt(irfeag),dratrawdd(irfeag), &
           ratraw(irniga),dratrawdt(irniga),dratrawdd(irniga))

! fe52(a,p)co55
      call rate_fe52ap(btemp,bden, &
           ratraw(irfeap),dratrawdt(irfeap),dratrawdd(irfeap), &
           ratraw(ircopa),dratrawdt(ircopa),dratrawdd(ircopa))

! co55(p,g)ni56
      call rate_co55pg(btemp,bden, &
           ratraw(ircopg),dratrawdt(ircopg),dratrawdd(ircopg), &
           ratraw(irnigp),dratrawdt(irnigp),dratrawdd(irnigp))

      return
      end




!------------------------------------------------------------------------------


      subroutine aprox13rattab
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! uses tables instead of analytical expressions to evaluate the
! raw reaction rates. a cubic polynomial is hardwired for speed.

      integer          i,j,imax,iat,mp,ifirst
      parameter        (mp = 4)
      double precision tlo,thi,tstp,zt, &
                       x,x1,x2,x3,x4,a,b,c,d,e,f,g,h,p,q, &
                       alfa,beta,gama,delt
      parameter        (thi = 10.0d0, tlo = 6.0d0, &
                        imax = int(thi-tlo)*250 + 1, &
                        tstp = (thi - tlo)/dfloat(imax-1))
      data             ifirst/0/


! read the table
      if (ifirst .eq. 0) then
       ifirst = 1
       if (imax .gt. nrattab) stop 'imax too small in aprox13tab'
       do i=1,imax
        zt = tlo + dfloat(i-1)*tstp
        zt = 10.0d0**(zt)
        ttab(i) = zt
       enddo
       open(unit=21,file='aprox13_rates.dat',status='old')
       read(21,129) ((rattab(j,i), j=1,nrat), i=1,imax)
129    format(1x,1p10e24.16)
       close(unit=21)
      end if


! normal execution starts here
! set the density dependence vector
      dtab(ircag)  = bden
      dtab(iroga)  = 1.0d0
      dtab(ir3a)   = bden*bden
      dtab(irg3a)  = 1.0d0
      dtab(ir1212) = bden
      dtab(ir1216) = bden
      dtab(ir1616) = bden
      dtab(iroag)  = bden
      dtab(irnega) = 1.0d0
      dtab(irneag) = bden
      dtab(irmgga) = 1.0d0
      dtab(irmgag) = bden
      dtab(irsiga) = 1.0d0
      dtab(irmgap) = bden
      dtab(iralpa) = bden
      dtab(iralpg) = bden
      dtab(irsigp) = 1.0d0
      dtab(irsiag) = bden
      dtab(irsga)  = 1.0d0
      dtab(irppa)  = bden
      dtab(irsiap) = bden
      dtab(irppg)  = bden
      dtab(irsgp)  = 1.0d0
      dtab(irsag)  = bden
      dtab(irarga) = 1.0d0
      dtab(irsap)  = bden
      dtab(irclpa) = bden
      dtab(irclpg) = bden
      dtab(irargp) = 1.0d0
      dtab(irarag) = bden
      dtab(ircaga) = 1.0d0
      dtab(irarap) = bden
      dtab(irkpa)  = bden
      dtab(irkpg)  = bden
      dtab(ircagp) = 1.0d0
      dtab(ircaag) = bden
      dtab(irtiga) = 1.0d0
      dtab(ircaap) = bden
      dtab(irscpa) = bden
      dtab(irscpg) = bden
      dtab(irtigp) = 1.0d0
      dtab(irtiag) = bden
      dtab(ircrga) = 1.0d0
      dtab(irtiap) = bden
      dtab(irvpa)  = bden
      dtab(irvpg)  = bden
      dtab(ircrgp) = 1.0d0
      dtab(ircrag) = bden
      dtab(irfega) = 1.0d0
      dtab(ircrap) = bden
      dtab(irmnpa) = bden
      dtab(irmnpg) = bden
      dtab(irfegp) = 1.0d0
      dtab(irfeag) = bden
      dtab(irniga) = 1.0d0
      dtab(irfeap) = bden
      dtab(ircopa) = bden
      dtab(ircopg) = bden
      dtab(irnigp) = 1.0d0


! hash locate the temperature
      iat = int((log10(btemp) - tlo)/tstp) + 1
      iat = max(1,min(iat - mp/2 + 1,imax - mp + 1))

! setup the lagrange interpolation coefficients for a cubic
      x  = btemp
      x1 = ttab(iat)
      x2 = ttab(iat+1)
      x3 = ttab(iat+2)
      x4 = ttab(iat+3)
      a  = x - x1
      b  = x - x2
      c  = x - x3
      d  = x - x4
      e  = x1 - x2
      f  = x1 - x3
      g  = x1 - x4
      h  = x2 - x3
      p  = x2 - x4
      q  = x3 - x4
      alfa =  b*c*d/(e*f*g)
      beta = -a*c*d/(e*h*p)
      gama =  a*b*d/(f*h*q)
      delt = -a*b*c/(g*p*q)

! crank off the raw reaction rates
      do j=1,nrat
       ratraw(j) = (alfa*rattab(j,iat) &
                  + beta*rattab(j,iat+1) &
                  + gama*rattab(j,iat+2) &
                  + delt*rattab(j,iat+3) &
                    ) * dtab(j)
      enddo
      return
      end


!------------------------------------------------------------------------------


      subroutine screen_aprox13(y)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine computes the screening factors
! and applies them to the raw reaction rates,
! producing the final reaction rates used by the
! right hand sides and jacobian matrix elements

! this routine assumes screen_on = 1 or = 0 has been set at a higher
! level presumably in the top level driver


! declare the pass
      double precision y(1)

! local variables
      integer          i,jscr,init
      double precision sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add, &
                       sc3a,sc3adt,sc3add

      double precision abar,zbar,z2bar,ytot1,zbarxx,z2barxx, &
                       denom,denomdt,denomdd,zz
      data             init/1/


! initialize
      ratdum(1:nrat)    = ratraw(1:nrat)
      dratdumdt(1:nrat) = dratrawdt(1:nrat)
      dratdumdd(1:nrat) = dratrawdd(1:nrat)
      scfac(1:nrat)     = 1.0d0
      dscfacdt(1:nrat)  = 0.0d0
      dscfacdt(1:nrat)  = 0.0d0


! with the passed composition, compute abar,zbar and other variables
      ytot1   = sum(y(ionbeg:ionend))
      zbarxx  = sum(zion(ionbeg:ionend)*y(ionbeg:ionend))
      z2barxx = sum(zion(ionbeg:ionend)*zion(ionbeg:ionend)*y(ionbeg:ionend))
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      z2bar  = z2barxx * abar



! first the always fun triple alpha and its inverse
      jscr = 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe4),aion(ihe4),4.0d0,8.0d0, &
                   jscr,init,sc2a,sc2adt,sc2add)

      sc3a   = sc1a * sc2a
      sc3adt = sc1adt*sc2a + sc1a*sc2adt
      sc3add = sc1add*sc2a + sc1a*sc2add

      ratdum(ir3a)    = ratraw(ir3a) * sc3a
      dratdumdt(ir3a) = dratrawdt(ir3a)*sc3a + ratraw(ir3a)*sc3adt
      dratdumdd(ir3a) = dratrawdd(ir3a)*sc3a + ratraw(ir3a)*sc3add

      scfac(ir3a)     = sc3a
      dscfacdt(ir3a)  = sc3adt
      dscfacdd(ir3a)  = sc3add


! c12 to o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircag)     = ratraw(ircag) * sc1a
      dratdumdt(ircag)  = dratrawdt(ircag)*sc1a + ratraw(ircag)*sc1adt
      dratdumdd(ircag)  = dratrawdd(ircag)*sc1a + ratraw(ircag)*sc1add

      scfac(ircag)      = sc1a
      dscfacdt(ircag)   = sc1adt
      dscfacdt(ircag)   = sc1add


! c12 + c12
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ic12),aion(ic12), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1212)    = ratraw(ir1212) * sc1a
      dratdumdt(ir1212) = dratrawdt(ir1212)*sc1a + ratraw(ir1212)*sc1adt
      dratdumdd(ir1212) = dratrawdd(ir1212)*sc1a + ratraw(ir1212)*sc1add

      scfac(ir1212)     = sc1a
      dscfacdt(ir1212)  = sc1adt
      dscfacdd(ir1212)  = sc1add



! c12 + o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(io16),aion(io16), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1216)    = ratraw(ir1216) * sc1a
      dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
      dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

      scfac(ir1216)     = sc1a
      dscfacdt(ir1216)  = sc1adt
      dscfacdd(ir1216)  = sc1add



! o16 + o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(io16),aion(io16),zion(io16),aion(io16), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1616)    = ratraw(ir1616) * sc1a
      dratdumdt(ir1616) = dratrawdt(ir1616)*sc1a + ratraw(ir1616)*sc1adt
      dratdumdd(ir1616) = dratrawdd(ir1616)*sc1a + ratraw(ir1616)*sc1add

      scfac(ir1616)     = sc1a
      dscfacdt(ir1616)  = sc1adt
      dscfacdd(ir1616)  = sc1add



! o16 to ne20
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(io16),aion(io16),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(iroag)    = ratraw(iroag) * sc1a
      dratdumdt(iroag) = dratrawdt(iroag)*sc1a + ratraw(iroag)*sc1adt
      dratdumdd(iroag) = dratrawdd(iroag)*sc1a + ratraw(iroag)*sc1add

      scfac(iroag)     = sc1a
      dscfacdt(iroag)  = sc1adt
      dscfacdd(iroag)  = sc1add



! ne20 to mg24
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ine20),aion(ine20),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irneag)    = ratraw(irneag) * sc1a
      dratdumdt(irneag) = dratrawdt(irneag)*sc1a + ratraw(irneag)*sc1adt
      dratdumdd(irneag) = dratrawdd(irneag)*sc1a + ratraw(irneag)*sc1add

      scfac(irneag)     = sc1a
      dscfacdt(irneag)  = sc1adt
      dscfacdd(irneag)  = sc1add


! mg24 to si28
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(img24),aion(img24),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irmgag)    = ratraw(irmgag) * sc1a
      dratdumdt(irmgag) = dratrawdt(irmgag)*sc1a + ratraw(irmgag)*sc1adt
      dratdumdd(irmgag) = dratrawdd(irmgag)*sc1a + ratraw(irmgag)*sc1add

      scfac(irmgag)     = sc1a
      dscfacdt(irmgag)  = sc1adt
      dscfacdd(irmgag)  = sc1add

      ratdum(irmgap)    = ratraw(irmgap) * sc1a
      dratdumdt(irmgap) = dratrawdt(irmgap)*sc1a + ratraw(irmgap)*sc1adt
      dratdumdd(irmgap) = dratrawdd(irmgap)*sc1a + ratraw(irmgap)*sc1add

      scfac(irmgap)     = sc1a
      dscfacdt(irmgap)  = sc1adt
      dscfacdd(irmgap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   13.0d0,27.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(iralpa)    = ratraw(iralpa) * sc1a
      dratdumdt(iralpa) = dratrawdt(iralpa)*sc1a + ratraw(iralpa)*sc1adt
      dratdumdd(iralpa) = dratrawdd(iralpa)*sc1a + ratraw(iralpa)*sc1add

      scfac(iralpa)     = sc1a
      dscfacdt(iralpa)  = sc1adt
      dscfacdd(iralpa)  = sc1add

      ratdum(iralpg)    = ratraw(iralpg) * sc1a
      dratdumdt(iralpg) = dratrawdt(iralpg)*sc1a + ratraw(iralpg)*sc1adt
      dratdumdd(iralpg) = dratrawdd(iralpg)*sc1a + ratraw(iralpg)*sc1add

      scfac(iralpg)     = sc1a
      dscfacdt(iralpg)  = sc1adt
      dscfacdd(iralpg)  = sc1add



! si28 to s32
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(isi28),aion(isi28),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irsiag)    = ratraw(irsiag) * sc1a
      dratdumdt(irsiag) = dratrawdt(irsiag)*sc1a + ratraw(irsiag)*sc1adt
      dratdumdd(irsiag) = dratrawdd(irsiag)*sc1a + ratraw(irsiag)*sc1add

      scfac(irsiag)     = sc1a
      dscfacdt(irsiag)  = sc1adt
      dscfacdd(irsiag)  = sc1add


      ratdum(irsiap)    = ratraw(irsiap) * sc1a
      dratdumdt(irsiap) = dratrawdt(irsiap)*sc1a + ratraw(irsiap)*sc1adt
      dratdumdd(irsiap) = dratrawdd(irsiap)*sc1a + ratraw(irsiap)*sc1add

      scfac(irsiap)     = sc1a
      dscfacdt(irsiap)  = sc1adt
      dscfacdd(irsiap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   15.0d0,31.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irppa)     = ratraw(irppa) * sc1a
      dratdumdt(irppa)  = dratrawdt(irppa)*sc1a  + ratraw(irppa)*sc1adt
      dratdumdd(irppa)  = dratrawdd(irppa)*sc1a  + ratraw(irppa)*sc1add

      scfac(irppa)      = sc1a
      dscfacdt(irppa)   = sc1adt
      dscfacdd(irppa)   = sc1add

      ratdum(irppg)     = ratraw(irppg) * sc1a
      dratdumdt(irppg)  = dratrawdt(irppg)*sc1a + ratraw(irppg)*sc1adt
      dratdumdd(irppg)  = dratrawdd(irppg)*sc1a + ratraw(irppg)*sc1add

      scfac(irppg)      = sc1a
      dscfacdt(irppg)   = sc1adt
      dscfacdd(irppg)   = sc1add



! s32 to ar36
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(is32),aion(is32),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irsag)     = ratraw(irsag) * sc1a
      dratdumdt(irsag)  = dratrawdt(irsag)*sc1a + ratraw(irsag)*sc1adt
      dratdumdd(irsag)  = dratrawdd(irsag)*sc1a + ratraw(irsag)*sc1add

      scfac(irsag)      = sc1a
      dscfacdt(irsag)   = sc1adt
      dscfacdd(irsag)   = sc1add

      ratdum(irsap)     = ratraw(irsap) * sc1a
      dratdumdt(irsap)  = dratrawdt(irsap)*sc1a + ratraw(irsap)*sc1adt
      dratdumdd(irsap)  = dratrawdd(irsap)*sc1a + ratraw(irsap)*sc1add

      scfac(irsap)      = sc1a
      dscfacdt(irsap)   = sc1adt
      dscfacdd(irsap)   = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   17.0d0,35.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irclpa)    = ratraw(irclpa) * sc1a
      dratdumdt(irclpa) = dratrawdt(irclpa)*sc1a + ratraw(irclpa)*sc1adt
      dratdumdd(irclpa) = dratrawdd(irclpa)*sc1a + ratraw(irclpa)*sc1add

      scfac(irclpa)     = sc1a
      dscfacdt(irclpa)  = sc1adt
      dscfacdt(irclpa)  = sc1add

      ratdum(irclpg)    = ratraw(irclpg) * sc1a
      dratdumdt(irclpg) = dratrawdt(irclpg)*sc1a + ratraw(irclpg)*sc1adt
      dratdumdd(irclpg) = dratrawdd(irclpg)*sc1a + ratraw(irclpg)*sc1add

      scfac(irclpg)     = sc1a
      dscfacdt(irclpg)  = sc1adt
      dscfacdd(irclpg)  = sc1add



! ar36 to ca40
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(iar36),aion(iar36),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irarag)    = ratraw(irarag) * sc1a
      dratdumdt(irarag) = dratrawdt(irarag)*sc1a + ratraw(irarag)*sc1adt
      dratdumdd(irarag) = dratrawdd(irarag)*sc1a + ratraw(irarag)*sc1add

      scfac(irarag)     = sc1a
      dscfacdt(irarag)  = sc1adt
      dscfacdd(irarag)  = sc1add

      ratdum(irarap)    = ratraw(irarap) * sc1a
      dratdumdt(irarap) = dratrawdt(irarap)*sc1a + ratraw(irarap)*sc1adt
      dratdumdd(irarap) = dratrawdd(irarap)*sc1a + ratraw(irarap)*sc1add

      scfac(irarap)     = sc1a
      dscfacdt(irarap)  = sc1adt
      dscfacdd(irarap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   19.0d0,39.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irkpa)     = ratraw(irkpa) * sc1a
      dratdumdt(irkpa)  = dratrawdt(irkpa)*sc1a  + ratraw(irkpa)*sc1adt
      dratdumdd(irkpa)  = dratrawdd(irkpa)*sc1a  + ratraw(irkpa)*sc1add

      scfac(irkpa)      = sc1a
      dscfacdt(irkpa)   = sc1adt
      dscfacdd(irkpa)   = sc1add

      ratdum(irkpg)     = ratraw(irkpg) * sc1a
      dratdumdt(irkpg)  = dratrawdt(irkpg)*sc1a  + ratraw(irkpg)*sc1adt
      dratdumdd(irkpg)  = dratrawdd(irkpg)*sc1a  + ratraw(irkpg)*sc1add

      scfac(irkpg)      = sc1a
      dscfacdt(irkpg)   = sc1adt
      dscfacdd(irkpg)   = sc1add



! ca40 to ti44
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ica40),aion(ica40),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircaag)    = ratraw(ircaag) * sc1a
      dratdumdt(ircaag) = dratrawdt(ircaag)*sc1a + ratraw(ircaag)*sc1adt
      dratdumdd(ircaag) = dratrawdd(ircaag)*sc1a + ratraw(ircaag)*sc1add

      scfac(ircaag)     = sc1a
      dscfacdt(ircaag)  = sc1adt
      dscfacdd(ircaag)  = sc1add

      ratdum(ircaap)    = ratraw(ircaap) * sc1a
      dratdumdt(ircaap) = dratrawdt(ircaap)*sc1a + ratraw(ircaap)*sc1adt
      dratdumdd(ircaap) = dratrawdd(ircaap)*sc1a + ratraw(ircaap)*sc1add

      scfac(ircaap)     = sc1a
      dscfacdt(ircaap)  = sc1adt
      dscfacdd(ircaap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   21.0d0,43.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irscpa)    = ratraw(irscpa) * sc1a
      dratdumdt(irscpa) = dratrawdt(irscpa)*sc1a + ratraw(irscpa)*sc1adt
      dratdumdd(irscpa) = dratrawdd(irscpa)*sc1a + ratraw(irscpa)*sc1add

      scfac(irscpa)     = sc1a
      dscfacdt(irscpa)  = sc1adt
      dscfacdd(irscpa)  = sc1add

      ratdum(irscpg)    = ratraw(irscpg) * sc1a
      dratdumdt(irscpg) = dratrawdt(irscpg)*sc1a + ratraw(irscpg)*sc1adt
      dratdumdd(irscpg) = dratrawdd(irscpg)*sc1a + ratraw(irscpg)*sc1add

      scfac(irscpg)     = sc1a
      dscfacdt(irscpg)  = sc1adt
      dscfacdd(irscpg)  = sc1add



! ti44 to cr48
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(iti44),aion(iti44),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irtiag)    = ratraw(irtiag) * sc1a
      dratdumdt(irtiag) = dratrawdt(irtiag)*sc1a + ratraw(irtiag)*sc1adt
      dratdumdd(irtiag) = dratrawdd(irtiag)*sc1a + ratraw(irtiag)*sc1add

      scfac(irtiag)     = sc1a
      dscfacdt(irtiag)  = sc1adt
      dscfacdd(irtiag)  = sc1add

      ratdum(irtiap)    = ratraw(irtiap) * sc1a
      dratdumdt(irtiap) = dratrawdt(irtiap)*sc1a + ratraw(irtiap)*sc1adt
      dratdumdd(irtiap) = dratrawdd(irtiap)*sc1a + ratraw(irtiap)*sc1add

      scfac(irtiap)  = sc1a
      dscfacdt(irtiap)  = sc1adt
      dscfacdd(irtiap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   23.0d0,47.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irvpa)     = ratraw(irvpa) * sc1a
      dratdumdt(irvpa)  = dratrawdt(irvpa)*sc1a  + ratraw(irvpa)*sc1adt
      dratdumdd(irvpa)  = dratrawdd(irvpa)*sc1a  + ratraw(irvpa)*sc1add

      scfac(irvpa)      = sc1a
      dscfacdt(irvpa)   = sc1adt
      dscfacdd(irvpa)   = sc1add

      ratdum(irvpg)     = ratraw(irvpg) * sc1a
      dratdumdt(irvpg)  = dratrawdt(irvpg)*sc1a  + ratraw(irvpg)*sc1adt
      dratdumdd(irvpg)  = dratrawdd(irvpg)*sc1a  + ratraw(irvpg)*sc1add

      scfac(irvpg)      = sc1a
      dscfacdt(irvpg)   = sc1adt
      dscfacdd(irvpg)   = sc1add



! cr48 to fe52
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(icr48),aion(icr48),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircrag)    = ratraw(ircrag) * sc1a
      dratdumdt(ircrag) = dratrawdt(ircrag)*sc1a + ratraw(ircrag)*sc1adt
      dratdumdd(ircrag) = dratrawdd(ircrag)*sc1a + ratraw(ircrag)*sc1add

      scfac(ircrag)     = sc1a
      dscfacdt(ircrag)  = sc1adt
      dscfacdd(ircrag)  = sc1add

      ratdum(ircrap)    = ratraw(ircrap) * sc1a
      dratdumdt(ircrap) = dratrawdt(ircrap)*sc1a + ratraw(ircrap)*sc1adt
      dratdumdd(ircrap) = dratrawdd(ircrap)*sc1a + ratraw(ircrap)*sc1add

      scfac(ircrap)     = sc1a
      dscfacdt(ircrap)  = sc1adt
      dscfacdd(ircrap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   25.0d0,51.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irmnpa)    = ratraw(irmnpa) * sc1a
      dratdumdt(irmnpa) = dratrawdt(irmnpa)*sc1a + ratraw(irmnpa)*sc1adt
      dratdumdd(irmnpa) = dratrawdd(irmnpa)*sc1a + ratraw(irmnpa)*sc1add

      scfac(irmnpa)     = sc1a
      dscfacdt(irmnpa)  = sc1adt
      dscfacdd(irmnpa)  = sc1add

      ratdum(irmnpg)    = ratraw(irmnpg) * sc1a
      dratdumdt(irmnpg) = dratrawdt(irmnpg)*sc1a + ratraw(irmnpg)*sc1adt
      dratdumdd(irmnpg) = dratrawdd(irmnpg)*sc1a + ratraw(irmnpg)*sc1add

      scfac(irmnpg)     = sc1a
      dscfacdt(irmnpg)  = sc1adt
      dscfacdd(irmnpg)  = sc1add


! fe52 to ni56
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ife52),aion(ife52),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irfeag)    = ratraw(irfeag) * sc1a
      dratdumdt(irfeag) = dratrawdt(irfeag)*sc1a + ratraw(irfeag)*sc1adt
      dratdumdd(irfeag) = dratrawdd(irfeag)*sc1a + ratraw(irfeag)*sc1add

      scfac(irfeag)     = sc1a
      dscfacdt(irfeag)  = sc1adt
      dscfacdd(irfeag)  = sc1add

      ratdum(irfeap) = ratraw(irfeap) * sc1a
      dratdumdt(irfeap) = dratrawdt(irfeap)*sc1a + ratraw(irfeap)*sc1adt
      dratdumdd(irfeap) = dratrawdd(irfeap)*sc1a + ratraw(irfeap)*sc1add

      scfac(irfeap)     = sc1a
      dscfacdt(irfeap)  = sc1adt
      dscfacdd(irfeap)  = sc1add

      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   27.0d0,55.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircopa)    = ratraw(ircopa) * sc1a
      dratdumdt(ircopa) = dratrawdt(ircopa)*sc1a + ratraw(ircopa)*sc1adt
      dratdumdd(ircopa) = dratrawdd(ircopa)*sc1a + ratraw(ircopa)*sc1add

      scfac(ircopa)     = sc1a
      dscfacdt(ircopa)  = sc1adt
      dscfacdd(ircopa)  = sc1add

      ratdum(ircopg)    = ratraw(ircopg) * sc1a
      dratdumdt(ircopg) = dratrawdt(ircopg)*sc1a + ratraw(ircopg)*sc1adt
      dratdumdd(ircopg) = dratrawdd(ircopg)*sc1a + ratraw(ircopg)*sc1add

      scfac(ircopg)     = sc1a
      dscfacdt(ircopg)  = sc1adt
      dscfacdd(ircopg)  = sc1add


! reset the screen initialization flag
      init = 0



! now form those lovely dummy proton link rates

      ratdum(irr1)     = 0.0d0
      dratdumdt(irr1)  = 0.0d0
      dratdumdd(irr1)  = 0.0d0
      denom    = ratdum(iralpa) + ratdum(iralpg)
      denomdt  = dratdumdt(iralpa) + dratdumdt(iralpg)
      denomdd  = dratdumdd(iralpa) + dratdumdd(iralpg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(irr1)    = ratdum(iralpa)*zz
       dratdumdt(irr1) = (dratdumdt(iralpa) - ratdum(irr1)*denomdt)*zz
       dratdumdd(irr1) = (dratdumdd(iralpa) - ratdum(irr1)*denomdd)*zz
      end if

      ratdum(irs1)     = 0.0d0
      dratdumdt(irs1)  = 0.0d0
      dratdumdd(irs1)  = 0.0d0
      denom    = ratdum(irppa) + ratdum(irppg)
      denomdt  = dratdumdt(irppa) + dratdumdt(irppg)
      denomdd  = dratdumdd(irppa) + dratdumdd(irppg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(irs1)    = ratdum(irppa)*zz
       dratdumdt(irs1) = (dratdumdt(irppa) - ratdum(irs1)*denomdt)*zz
       dratdumdd(irs1) = (dratdumdd(irppa) - ratdum(irs1)*denomdd)*zz
      end if

      ratdum(irt1)     = 0.0d0
      dratdumdt(irt1)  = 0.0d0
      dratdumdd(irt1)  = 0.0d0
      denom    = ratdum(irclpa) + ratdum(irclpg)
      denomdt  = dratdumdt(irclpa) + dratdumdt(irclpg)
      denomdd  = dratdumdd(irclpa) + dratdumdd(irclpg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(irt1)    = ratdum(irclpa)*zz
       dratdumdt(irt1) = (dratdumdt(irclpa) - ratdum(irt1)*denomdt)*zz
       dratdumdd(irt1) = (dratdumdd(irclpa) - ratdum(irt1)*denomdd)*zz
      end if

      ratdum(iru1)     = 0.0d0
      dratdumdt(iru1)  = 0.0d0
      dratdumdd(iru1)  = 0.0d0
      denom    = ratdum(irkpa) + ratdum(irkpg)
      denomdt  = dratdumdt(irkpa) + dratdumdt(irkpg)
      denomdd  = dratdumdd(irkpa) + dratdumdd(irkpg)
      if (denom .ne. 0.0) then
       zz   = 1.0d0/denom
       ratdum(iru1)   = ratdum(irkpa)*zz
       dratdumdt(iru1) = (dratdumdt(irkpa) - ratdum(iru1)*denomdt)*zz
       dratdumdd(iru1) = (dratdumdd(irkpa) - ratdum(iru1)*denomdd)*zz
      end if

      ratdum(irv1)     = 0.0d0
      dratdumdt(irv1)  = 0.0d0
      dratdumdd(irv1)  = 0.0d0
      denom    = ratdum(irscpa) + ratdum(irscpg)
      denomdt  = dratdumdt(irscpa) + dratdumdt(irscpg)
      denomdd  = dratdumdd(irscpa) + dratdumdd(irscpg)
      if (denom .ne. 0.0) then
       zz  = 1.0d0/denom
       ratdum(irv1)    = ratdum(irscpa)*zz
       dratdumdt(irv1) = (dratdumdt(irscpa) - ratdum(irv1)*denomdt)*zz
       dratdumdd(irv1) = (dratdumdd(irscpa) - ratdum(irv1)*denomdd)*zz
      end if

      ratdum(irw1)    = 0.0d0
      dratdumdt(irw1) = 0.0d0
      dratdumdd(irw1) = 0.0d0
      denom    = ratdum(irvpa) + ratdum(irvpg)
      denomdt  = dratdumdt(irvpa) + dratdumdt(irvpg)
      denomdd  = dratdumdd(irvpa) + dratdumdd(irvpg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(irw1)    = ratdum(irvpa)*zz
       dratdumdt(irw1) = (dratdumdt(irvpa) - ratdum(irw1)*denomdt)*zz
       dratdumdd(irw1) = (dratdumdd(irvpa) - ratdum(irw1)*denomdd)*zz
      end if

      ratdum(irx1)    = 0.0d0
      dratdumdt(irx1) = 0.0d0
      dratdumdd(irx1) = 0.0d0
      denom    = ratdum(irmnpa) + ratdum(irmnpg)
      denomdt  = dratdumdt(irmnpa) + dratdumdt(irmnpg)
      denomdd  = dratdumdd(irmnpa) + dratdumdd(irmnpg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(irx1)    = ratdum(irmnpa)*zz
       dratdumdt(irx1) = (dratdumdt(irmnpa) - ratdum(irx1)*denomdt)*zz
       dratdumdd(irx1) = (dratdumdd(irmnpa) - ratdum(irx1)*denomdd)*zz
      endif

      ratdum(iry1)    = 0.0d0
      dratdumdt(iry1) = 0.0d0
      dratdumdd(iry1) = 0.0d0
      denom    = ratdum(ircopa) + ratdum(ircopg)
      denomdt  = dratdumdt(ircopa) + dratdumdt(ircopg)
      denomdd  = dratdumdd(ircopa) + dratdumdd(ircopg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(iry1)    = ratdum(ircopa)*zz
       dratdumdt(iry1) = (dratdumdt(ircopa) - ratdum(iry1)*denomdt)*zz
       dratdumdd(iry1) = (dratdumdd(ircopa) - ratdum(iry1)*denomdd)*zz
      end if

      return
      end





      subroutine init_burner
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! this routine initializes stuff for the aprox13 network
!
! declare
      integer          i
      double precision  mev2erg,mev2gr
      parameter        (mev2erg = ev2erg*1.0d6, &
                        mev2gr  = mev2erg/clight**2)


! for easy zeroing of the isotope pointers
      integer          isotp(nisotp)
      equivalence      (isotp(1),ih1)

! for easy zeroing of the rate pointers
      integer          rts(numrates)
      equivalence      (rts(1),ir3a)


! zero all the isotope and rate pointers
      isotp(:) = 0
      rts(:)   = 0


! set the size of the network and the number of rates
      idnet   = idaprox13
      ionmax  = 13
      iener   = ionmax + 1
      itemp   = ionmax + 2
!      iener   = 1
!      itemp   = 2
      neqs    = 15
      nrat    = 67
      netname = 'aprox13'


! set the id numbers of the elements
      ionbeg = 1
      ionend = 13

      ihe4  = ionbeg
      ic12  = ionbeg + 1
      io16  = ionbeg + 2
      ine20 = ionbeg + 3
      img24 = ionbeg + 4
      isi28 = ionbeg + 5
      is32  = ionbeg + 6
      iar36 = ionbeg + 7
      ica40 = ionbeg + 8
      iti44 = ionbeg + 9
      icr48 = ionbeg + 10
      ife52 = ionbeg + 11
      ini56 = ionbeg + 12


! set the names of the elements
      ionam(ihe4)  = 'he4  '
      ionam(ic12)  = 'c12  '
      ionam(io16)  = 'o16  '
      ionam(ine20) = 'ne20 '
      ionam(img24) = 'mg24 '
      ionam(isi28) = 'si28 '
      ionam(is32)  = 's32  '
      ionam(iar36) = 'ar36 '
      ionam(ica40) = 'ca40 '
      ionam(iti44) = 'ti44 '
      ionam(icr48) = 'cr48 '
      ionam(ife52) = 'fe52 '
      ionam(ini56) = 'ni56 '

! set the number of nucleons in the element
      aion(ihe4)  = 4.0d0
      aion(ic12)  = 12.0d0
      aion(io16)  = 16.0d0
      aion(ine20) = 20.0d0
      aion(img24) = 24.0d0
      aion(isi28) = 28.0d0
      aion(is32)  = 32.0d0
      aion(iar36) = 36.0d0
      aion(ica40) = 40.0d0
      aion(iti44) = 44.0d0
      aion(icr48) = 48.0d0
      aion(ife52) = 52.0d0
      aion(ini56) = 56.0d0

! set the number of protons in the element
      zion(ihe4)  = 2.0d0
      zion(ic12)  = 6.0d0
      zion(io16)  = 8.0d0
      zion(ine20) = 10.0d0
      zion(img24) = 12.0d0
      zion(isi28) = 14.0d0
      zion(is32)  = 16.0d0
      zion(iar36) = 18.0d0
      zion(ica40) = 20.0d0
      zion(iti44) = 22.0d0
      zion(icr48) = 24.0d0
      zion(ife52) = 26.0d0
      zion(ini56) = 28.0d0


! set the binding energy of the element
      bion(ihe4)  =  28.29603d0
      bion(ic12)  =  92.16294d0
      bion(io16)  = 127.62093d0
      bion(ine20) = 160.64788d0
      bion(img24) = 198.25790d0
      bion(isi28) = 236.53790d0
      bion(is32)  = 271.78250d0
      bion(iar36) = 306.72020d0
      bion(ica40) = 342.05680d0
      bion(iti44) = 375.47720d0
      bion(icr48) = 411.46900d0
      bion(ife52) = 447.70800d0
      bion(ini56) = 484.00300d0


! set the number of neutrons and mass
      nion(ionbeg:ionend) = aion(ionbeg:ionend) - zion(ionbeg:ionend)

! mass of each isotope
      mion(ionbeg:ionend) = nion(ionbeg:ionend)*mn + zion(ionbeg:ionend)*mp - bion(ionbeg:ionend)*mev2gr

! molar mass
!      wion(ionbeg:ionend) = avo * mion(ionbeg:ionend)

! a common approximation
      wion(ionbeg:ionend) = aion(ionbeg:ionend)

! set the partition functions - statistical weights, ground-state only
      wpart(ionbeg:ionend) = 1.0d0


! set the id numbers of the reaction rates
      ir3a   = 1
      irg3a  = 2
      ircag  = 3
      ir1212 = 4
      ir1216 = 5
      ir1616 = 6
      iroga  = 7
      iroag  = 8
      irnega = 9
      irneag = 10
      irmgga = 11
      irmgag = 12
      irsiga = 13
      irmgap = 14
      iralpa = 15
      iralpg = 16
      irsigp = 17
      irsiag = 18
      irsga  = 19
      irsiap = 20
      irppa  = 21
      irppg  = 22
      irsgp  = 23
      irsag  = 24
      irarga = 25
      irsap  = 26
      irclpa = 27
      irclpg = 28
      irargp = 29
      irarag = 30
      ircaga = 31
      irarap = 32
      irkpa  = 33
      irkpg  = 34
      ircagp = 35
      ircaag = 36
      irtiga = 37
      ircaap = 38
      irscpa = 39
      irscpg = 40
      irtigp = 41
      irtiag = 42
      ircrga = 43
      irtiap = 44
      irvpa  = 45
      irvpg  = 46
      ircrgp = 47
      ircrag = 48
      irfega = 49
      ircrap = 50
      irmnpa = 51
      irmnpg = 52
      irfegp = 53
      irfeag = 54
      irniga = 55
      irfeap = 56
      ircopa = 57
      ircopg = 58
      irnigp = 59

      irr1   = 60
      irs1   = 61
      irt1   = 62
      iru1   = 63
      irv1   = 64
      irw1   = 65
      irx1   = 66
      iry1   = 67


! set the names of the reaction rates
      ratnam(ir3a)   = 'r3a  '
      ratnam(irg3a)  = 'rg3a '
      ratnam(ircag)  = 'rcag '
      ratnam(ir1212) = 'r1212'
      ratnam(ir1216) = 'r1216'
      ratnam(ir1616) = 'r1616'
      ratnam(iroga)  = 'roga '
      ratnam(iroag)  = 'roag '
      ratnam(irnega) = 'rnega'
      ratnam(irneag) = 'rneag'
      ratnam(irmgga) = 'rmgga'
      ratnam(irmgag) = 'rmgag'
      ratnam(irsiga) = 'rsiga'
      ratnam(irmgap) = 'rmgap'
      ratnam(iralpa) = 'ralpa'
      ratnam(iralpg) = 'ralpg'
      ratnam(irsigp) = 'rsigp'
      ratnam(irsiag) = 'rsiag'
      ratnam(irsga)  = 'rsga '
      ratnam(irsiap) = 'rsiap'
      ratnam(irppa)  = 'rppa '
      ratnam(irppg)  = 'rppg '
      ratnam(irsgp)  = 'rsgp '
      ratnam(irsag)  = 'rsag '
      ratnam(irarga) = 'rarga'
      ratnam(irsap)  = 'rsap '
      ratnam(irclpa) = 'rclpa'
      ratnam(irclpg) = 'rclpg'
      ratnam(irargp) = 'rargp'
      ratnam(irarag) = 'rarag'
      ratnam(ircaga) = 'rcaga'
      ratnam(irarap) = 'rarap'
      ratnam(irkpa)  = 'rkpa '
      ratnam(irkpg)  = 'rkpg '
      ratnam(ircagp) = 'rcagp'
      ratnam(ircaag) = 'rcaag'
      ratnam(irtiga) = 'rtiga'
      ratnam(ircaap) = 'rcaap'
      ratnam(irscpa) = 'rscpa'
      ratnam(irscpg) = 'rscpg'
      ratnam(irtigp) = 'rtigp'
      ratnam(irtiag) = 'rtiag'
      ratnam(ircrga) = 'rcrga'
      ratnam(irtiap) = 'rtiap'
      ratnam(irvpa)  = 'rvpa '
      ratnam(irvpg)  = 'rvpg '
      ratnam(ircrgp) = 'rcrgp'
      ratnam(ircrag) = 'rcrag'
      ratnam(irfega) = 'rfega'
      ratnam(ircrap) = 'rcrap'
      ratnam(irmnpa) = 'rmnpa'
      ratnam(irmnpg) = 'rmnpg'
      ratnam(irfegp) = 'rfegp'
      ratnam(irfeag) = 'rfeag'
      ratnam(irniga) = 'rniga'
      ratnam(irfeap) = 'rfeap'
      ratnam(ircopa) = 'rcopa'
      ratnam(ircopg) = 'rcopg'
      ratnam(irnigp) = 'rnigp'

      ratnam(irr1)   = 'r1   '
      ratnam(irs1)   = 's1   '
      ratnam(irt1)   = 't1   '
      ratnam(iru1)   = 'u1   '
      ratnam(irv1)   = 'v1   '
      ratnam(irw1)   = 'w1   '
      ratnam(irx1)   = 'x1   '
      ratnam(iry1)   = 'y1   '


! we'll initialze the network options in this routine as well
! most of these are meaningless in this stripped down version of aprox13

! general options
      screen_on      = 0
      use_tables     = 0
      weak_on        = 0
      ffn_on         = 0
      pure_network   = 0
      nse_analysis   = 0
      allow_nse_evol = 0


! printing information
      iprint_files  = 0
      iprint_screen = 0


! inititailize the burn type logicals
      one_step             = .false.
      hydrostatic          = .false.
      expansion            = .false.
      self_heat_const_den  = .true.
      self_heat_const_pres = .false.
      pt_hist              = .false.
      bbang                = .false.
      detonation           = .false.
      trho_hist            = .false.


! adiabatic expansion off
      psi       = 0.0d0
      temp_stop = 1.0d30


! isotope stopping criteria
       name_stop = 'he4 '
       xmass_stop = -1.0d30


! mass fractions above sthreshold are written to the summary file
      sthreshold = 1.0d30

      return
      end



!---------------------------------------------------------------------
      subroutine netint(start,stptry,stpmin,stopp,bc, &
                        eps,ylogi,nok,nbad,kount,odescal, &
                        derivs,jakob,steper)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'



! ode integrator for stiff odes
! tuned for nnuclear reacton networks

! input:
! start    = beginning integration point
! stptry   = suggested first step size
! stpmin   = minimum allowable step size
! stopp    = ending integration point
! bc       = initial conditions, array of physical dimension yphys
! eps      = desired fraction error during the integration
! odescal  = error scaling factor
! derivs   = name of the routine that contains the odes
! jakob    = name of the routine that contains the jacobian of the odes
! steper   = name of the routine that will take a single step

! output:
! nok      = number of succesful steps taken
! nbad     = number of bad steps taken, bad but retried and then succesful
! kount    = total number of steps taken



! declare the pass
      external         derivs,jakob,steper
      integer          ylogi,nok,nbad,kount
      double precision start,stptry,stpmin,stopp,bc(ylogi),eps,odescal


! local variables
      character*5      cdtname
      integer          nmax,stpmax,i,ii,nstp,idt,ioff
      parameter        (nmax = abignet*nzmax, stpmax=200000)
      double precision yscal(nmax),y(nmax),dydx(nmax),xdum(nmax), &
                       asum,cons,x,h,hdid,hnext,tiny
      parameter        (tiny=1.0d-15)


! for some more informative printouts
      double precision abar,zbar,wbar,ye,xcess


! here are the format statements for printouts as we integrate
100   format(1x,i6,' ',a,1pe11.4,a,a,a,1pe11.4, &
                3(a6,1pe10.3),5(a5,1pe9.2))
101   format(1x,1p12e10.2)



! initialize
      if (ylogi  .gt. nmax) stop 'ylogi > nmax in routine netint'
      x      = start
      h      = sign(stptry,stopp-start)
      nok    = 0
      nbad   = 0
      kount  = 0
      idt    = 0
      ioff   = ionbeg - 1 


! store the first step
      y(1:ylogi) = bc(1:ylogi)


! take at most stpmax steps
      do nstp=1,stpmax


! positive definite abundance fractions
       do i=ionbeg,ionend
        y(i) = min(1.0d0, max(y(i),1.0d-30))
       enddo



! get the right hand sides
        call derivs(x,y,dydx)


! scaling vector used to monitor accuracy
        do i=1,ylogi
         yscal(i) = max(odescal,abs(y(i)))
        enddo


! store intermediate results
       kount         = kount+1


! screen print
       if (iprint_screen .eq. 1) then


! form the mass fractions and nonconservation
        xdum(ionbeg:ionend) = y(ionbeg:ionend) * aion(ionbeg:ionend)
        cons = 1.0d0 - sum(xdum(ionbeg:ionend))

        if (idt .eq. 0) then
         cdtname = 'time '
        else
         cdtname = ionam(idt)
        end if

        call indexx(ionmax,xdum(ionbeg),izwork1(ionbeg))
        write(6,100) kount,' time',x, &
                     ' dt(',cdtname,')',hdid, &
                     (ionam(izwork1(ii)+ioff), xdum(izwork1(ii)+ioff),ii=ionend,ionend-2,-1), &
                     ' temp',btemp,' den ',bden,' enuc',sdot-sneut, &
                     ' ye  ',ye,' sum ',cons
       end if


! if the step can overshoot the stop point cut it
       if ((x+h-stopp)*(x+h-start) .gt. 0.0d0) h = stopp - x


! do an integration step
       call steper(y,dydx,ylogi,x,h,eps,yscal,hdid,hnext, &
                   derivs,jakob,idt)

       if (hdid.eq.h) then
        nok = nok+1
       else
        nbad = nbad+1
       end if



! this is the normal exit point, save the final step
       if ( (nstp .eq. stpmax)                 .or. &
            (x-stopp)*(stopp-start).ge. 0.0d0  ) then

        bc(1:ylogi) = y(1:ylogi)
        kount = kount+1


! screen print
        if (iprint_screen .eq. 1) then
         xdum(ionbeg:ionend) = y(ionbeg:ionend) * aion(ionbeg:ionend)
         cons = 1.0d0 - sum(xdum(ionbeg:ionend))
         call indexx(ionmax,xdum(ionbeg),izwork1(ionbeg))
         write(6,100) kount,' time',x, &
                     ' dt(',cdtname,')',hdid, &
                     (ionam(izwork1(ii)+ioff), xdum(izwork1(ii)+ioff),ii=ionend,ionend-2,-1), &
                     ' temp',btemp,' den ',bden,' enuc',sdot-sneut, &
                     ' ye  ',ye,' sum ',cons
        end if
        return
       end if


! set the step size for the next iteration; stay above stpmin

       h = hnext

       if (abs(h).lt.stpmin) then
        write(6,*) 'integration failed'
        write(6,*) 'nstp =',nstp
        write(6,*) 'nok nbad =',nok,nbad
        write(6,*) 'attempted time step =',stptry
        write(6,*) 'current time step =',h
        write(6,*) 'temperature btemp =',btemp
        write(6,*) 'density bden =',bden
        write(6,*) 'input composition:'
        write(6,*) (bc(i), i=ionbeg,ionend)
        write(6,*) 'current composition:'
        write(6,*) (y(i), i=ionbeg,ionend)
        stop 'h < stpmin in netint'
       end if


! back for another iteration or death
      enddo
      stop 'more than stpmax steps required in netint'
      end




!---------------------------------------------------------------------
!
! this routine contains auxillary network routine

! routine screen5 computes screening factors
! routine sneut5 computes neutrino loss rates
! routine ifermi12 does an inverse fermi integral of order 1/2
! routine zfermim12 does an inverse fermi integral of order -1/2


      subroutine screen5(temp,den,zbar,abar,z2bar, &
                         z1,a1,z2,a2,jscreen,init, &
                         scor,scordt,scordd)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! this subroutine calculates screening factors and their derivatives
! for nuclear reaction rates in the weak, intermediate and strong regimes.
! based on graboske, dewit, grossman and cooper apj 181 457 1973 for
! weak screening. based on alastuey and jancovici apj 226 1034 1978,
! with plasma parameters from itoh et al apj 234 1079 1979, for strong
! screening.

! input:
! temp    = temperature
! den     = density
! zbar    = mean charge per nucleus
! abar    = mean number of nucleons per nucleus
! z2bar   = mean square charge per nucleus
! z1 a1   = charge and number in the entrance channel
! z2 a2   = charge and number in the exit channel
! jscreen = counter of which reaction is being calculated
! init    = flag to compute the more expensive functions just once

! output:
! scor    = screening correction
! scordt  = derivative of screening correction with temperature
! scordd  = derivative of screening correction with density


! declare the pass
      integer          jscreen,init
      double precision temp,den,zbar,abar,z2bar,z1,a1,z2,a2, &
                       scor,scordt,scordd


! local variables
      double precision aa,daadt,daadd,bb,cc,dccdt,dccdd, &
                       pp,dppdt,dppdd,qq,dqqdt,dqqdd,rr,drrdt,drrdd, &
                       ss,dssdt,dssdd,tt,dttdt,dttdd,uu,duudt,duudd, &
                       vv,dvvdt,dvvdd,a3,da3,tempi,dtempi,deni, &
                       qlam0z,qlam0zdt,qlam0zdd, &
                       h12w,dh12wdt,dh12wdd,h12,dh12dt,dh12dd, &
                       h12x,dh12xdt,dh12xdd,alfa,beta, &
                       taufac,taufacdt,gamp,gampdt,gampdd, &
                       gamef,gamefdt,gamefdd, &
                       tau12,tau12dt,alph12,alph12dt,alph12dd, &
                       xlgfac,dxlgfacdt,dxlgfacdd, &
                       gamp14,gamp14dt,gamp14dd, &
                       xni,dxnidd,ytot, &
                       temp_old,den_old,zbar_old,abar_old


! screening variables
! zs13    = (z1+z2)**(1./3.)
! zhat    = combination of z1 and z2 raised to the 5/3 power
! zhat2   = combination of z1 and z2 raised to the 5/12 power
! lzav    = log of effective charge
! aznut   = combination of a1,z1,a2,z2 raised to 1/3 power


      integer          nscreen_max
      parameter        (nscreen_max = 2*abignet + 40)

      double precision zs13(nscreen_max),zhat(nscreen_max), &
                       zhat2(nscreen_max),lzav(nscreen_max), &
                       aznut(nscreen_max),zs13inv(nscreen_max)


! parameter fact is the cube root of 2
      double precision  x13,x14,x53,x532,x512,fact,co2,gamefx,gamefs, &
                        blend_frac
      parameter        (x13    = 1.0d0/3.0d0, &
                        x14    = 1.0d0/4.0d0, &
                        x53    = 5.0d0/3.0d0, &
                        x532   = 5.0d0/32.0d0, &
                        x512   = 5.0d0/12.0d0, &
                        fact   = 1.25992104989487d0, &
                        co2    = x13 * 4.248719d3, &
                        gamefx = 0.3d0, &
                        gamefs = 0.8d0, &
                        blend_frac = 0.05d0)


      data     temp_old/-1.0d0/, den_old/-1.0d0/, &
               zbar_old/-1.0d0/, abar_old/-1.0d0/




! compute and store the more expensive screening factors
      if (init .eq. 1) then
       if (jscreen .gt. nscreen_max) &
       stop 'jscreen > nscreen_max in screen5'
       zs13(jscreen)    = (z1 + z2)**x13
       zs13inv(jscreen) = 1.0d0/zs13(jscreen)
       zhat(jscreen)    = (z1 + z2)**x53  - z1**x53 - z2**x53
       zhat2(jscreen)   = (z1 + z2)**x512 - z1**x512 -z2**x512
       lzav(jscreen)    = x53 * log(z1*z2/(z1 + z2))
       aznut(jscreen)   = (z1**2 * z2**2 * a1*a2 / (a1 + a2))**x13
      endif


! calculate average plasma, if need be
      if (temp_old .ne. temp .or. &
          den_old  .ne. den  .or. &
          zbar_old  .ne. zbar  .or. &
          abar_old  .ne. abar ) then

       temp_old = temp
       den_old  = den
       zbar_old  = zbar
       abar_old  = abar

       ytot     = 1.0d0/abar
       rr       = den * ytot
       tempi   = 1.0d0/temp
       dtempi  = -tempi*tempi
       deni    = 1.0d0/den

       pp       = sqrt(rr*tempi*(z2bar + zbar))
       qq       = 0.5d0/pp *(z2bar + zbar)
       dppdt    = qq*rr*dtempi
       dppdd    = qq*ytot*tempi

       qlam0z   = 1.88d8 * tempi * pp
       qlam0zdt = 1.88d8 * (dtempi*pp + tempi*dppdt)
       qlam0zdd = 1.88d8 * tempi * dppdd

       taufac   = co2 * tempi**x13
       taufacdt = -x13*taufac*tempi

       qq      = rr*zbar
       xni     = qq**x13
       dxnidd  = x13 * xni * deni

       aa     = 2.27493d5 * tempi * xni
       daadt  = 2.27493d5 * dtempi * xni
       daadd  = 2.27493d5 * tempi * dxnidd
      end if


! calculate individual screening factors
      bb       = z1 * z2
      gamp     = aa
      gampdt   = daadt
      gampdd   = daadd

      qq       = fact * bb * zs13inv(jscreen)
      gamef    = qq * gamp
      gamefdt  = qq * gampdt
      gamefdd  = qq * gampdd

      tau12    = taufac * aznut(jscreen)
      tau12dt  = taufacdt * aznut(jscreen)

      qq       = 1.0d0/tau12
      alph12   = gamef * qq
      alph12dt = (gamefdt - alph12*tau12dt) * qq
      alph12dd = gamefdd * qq



! limit alph12 to 1.6 to prevent unphysical behavior.
! this should really be replaced by a pycnonuclear reaction rate formula
      if (alph12 .gt. 1.6) then
       alph12   = 1.6d0
       alph12dt = 0.0d0
       alph12dd = 0.0d0

       gamef    = 1.6d0 * tau12
       gamefdt  = 1.6d0 * tau12dt
       gamefdd  = 0.0d0

       qq       = zs13(jscreen)/(fact * bb)
       gamp     = gamef * qq
       gampdt   = gamefdt * qq
       gampdd   = 0.0d0
      end if



! weak screening regime
      h12w    = bb * qlam0z
      dh12wdt = bb * qlam0zdt
      dh12wdd = bb * qlam0zdd

      h12     = h12w
      dh12dt  = dh12wdt
      dh12dd  = dh12wdd



! intermediate and strong sceening regime
      if (gamef .gt. gamefx) then

       gamp14   = gamp**x14
       rr       = 1.0d0/gamp
       qq       = 0.25d0*gamp14*rr
       gamp14dt = qq * gampdt
       gamp14dd = qq * gampdd

       cc       =   0.896434d0 * gamp * zhat(jscreen) &
                  - 3.44740d0  * gamp14 * zhat2(jscreen) &
                  - 0.5551d0   * (log(gamp) + lzav(jscreen)) &
                  - 2.996d0

       dccdt    =   0.896434d0 * gampdt * zhat(jscreen) &
                  - 3.44740d0  * gamp14dt * zhat2(jscreen) &
                  - 0.5551d0*rr*gampdt

       dccdd    =   0.896434d0 * gampdd * zhat(jscreen) &
                  - 3.44740d0  * gamp14dd * zhat2(jscreen) &
                  - 0.5551d0*rr*gampdd

       a3     = alph12 * alph12 * alph12
       da3    = 3.0d0 * alph12 * alph12

       qq     = 0.014d0 + 0.0128d0*alph12
       dqqdt  = 0.0128d0*alph12dt
       dqqdd  = 0.0128d0*alph12dd

       rr     = x532 - alph12*qq
       drrdt  = -(alph12dt*qq + alph12*dqqdt)
       drrdd  = -(alph12dd*qq + alph12*dqqdd)

       ss     = tau12*rr
       dssdt  = tau12dt*rr + tau12*drrdt
       dssdd  = tau12*drrdd

       tt     =  -0.0098d0 + 0.0048d0*alph12
       dttdt  = 0.0048d0*alph12dt
       dttdd  = 0.0048d0*alph12dd

       uu     =  0.0055d0 + alph12*tt
       duudt  = alph12dt*tt + alph12*dttdt
       duudd  = alph12dd*tt + alph12*dttdd

       vv   = gamef * alph12 * uu
       dvvdt= gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt
       dvvdd= gamefdd*alph12*uu + gamef*alph12dd*uu + gamef*alph12*duudd

       h12     = cc - a3 * (ss + vv)
       rr      = da3 * (ss + vv)
       dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt)
       dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd)

       rr     =  1.0d0 - 0.0562d0*a3
       ss     =  -0.0562d0*da3
       drrdt  = ss*alph12dt
       drrdd  = ss*alph12dd

       if (rr .ge. 0.77d0) then
        xlgfac    = rr
        dxlgfacdt = drrdt
        dxlgfacdd = drrdd
       else
        xlgfac    = 0.77d0
        dxlgfacdt = 0.0d0
        dxlgfacdd = 0.0d0
       end if


       h12    = log(xlgfac) + h12
       rr     = 1.0d0/xlgfac
       dh12dt = rr*dxlgfacdt + dh12dt
       dh12dd = rr*dxlgfacdd + dh12dd


       if (gamef .le. gamefs) then
        rr     =  2.0d0*(gamefs - gamef)
        drrdt  = -2.0d0*gamefdt
        drrdd  = -2.0d0*gamefdd

        ss     = 2.0d0*(gamef - gamefx)
        dssdt  = 2.0d0*gamefdt
        dssdd  = 2.0d0*gamefdd


! store current values for possible blending
        h12x    = h12
        dh12xdt = dh12dt
        dh12xdd = dh12dd

        vv     = h12
        h12    = h12w*rr + vv*ss
        dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt
        dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd

! blend the transition region - from bill paxton
       if (gamefs - gamef .lt. blend_frac*(gamefs - gamefx)) then
         alfa   = (gamefs - gamef) / (blend_frac*(gamefs - gamefx))
         alfa   = 0.5d0 * (1d0 - cos(pi*alfa))
         beta   = 1.0d0 - alfa
         h12    = alfa * h12 + beta * h12x
         dh12dt = alfa * dh12dt + beta * dh12xdt
         dh12dd = alfa * dh12dd + beta * dh12xdd
        end if
       end if


! end of intermediate and strong screening if
      end if


! machine limit the output
      h12    = max(min(h12,300.0d0),0.0d0)
      scor   = exp(h12)
      if (h12 .eq. 300.0d0) then
       scordt = 0.0d0
       scordd = 0.0d0
      else
       scordt = scor * dh12dt
       scordd = scor * dh12dd
      end if

!      write(6,111) 'weak =',h12w,' total =',h12,
!     1             ' 1-ratio =',1.0d0-h12w/h12,' correction',scor
! 111  format(1x,4(a,1pe13.6))
!      read(5,*)

      return
      end



      subroutine sneut5(temp,den,abar,zbar, &
                        snu,dsnudt,dsnudd,dsnuda,dsnudz)
      include 'implno.dek'
      include 'const.dek'

! this routine computes neutrino losses from the analytic fits of
! itoh et al. apjs 102, 411, 1996, and also returns their derivatives.

! input:
! temp = temperature
! den  = density
! abar = mean atomic weight
! zbar = mean charge

! output:
! snu    = total neutrino loss rate in erg/g/sec
! dsnudt = derivative of snu with temperature
! dsnudd = derivative of snu with density
! dsnuda = derivative of snu with abar
! dsnudz = derivative of snu with zbar


! declare the pass
      double precision temp,den,abar,zbar, &
                       snu,dsnudt,dsnudd,dsnuda,dsnudz

! local variables
      double precision spair,spairdt,spairdd,spairda,spairdz, &
                       splas,splasdt,splasdd,splasda,splasdz, &
                       sphot,sphotdt,sphotdd,sphotda,sphotdz, &
                       sbrem,sbremdt,sbremdd,sbremda,sbremdz, &
                       sreco,srecodt,srecodd,srecoda,srecodz

      double precision t9,xl,xldt,xlp5,xl2,xl3,xl4,xl5,xl6,xl7,xl8,xl9, &
                       xlmp5,xlm1,xlm2,xlm3,xlm4,xlnt,cc,den6,tfermi, &
                       a0,a1,a2,a3,b1,b2,c00,c01,c02,c03,c04,c05,c06, &
                       c10,c11,c12,c13,c14,c15,c16,c20,c21,c22,c23,c24, &
                       c25,c26,dd00,dd01,dd02,dd03,dd04,dd05,dd11,dd12, &
                       dd13,dd14,dd15,dd21,dd22,dd23,dd24,dd25,b,c,d,f0, &
                       f1,deni,tempi,abari,zbari,f2,f3,z,xmue,ye, &
                       dum,dumdt,dumdd,dumda,dumdz, &
                       gum,gumdt,gumdd,gumda,gumdz


! pair production
      double precision rm,rmdd,rmda,rmdz,rmi,gl,gldt, &
                       zeta,zetadt,zetadd,zetada,zetadz,zeta2,zeta3, &
                       xnum,xnumdt,xnumdd,xnumda,xnumdz, &
                       xden,xdendt,xdendd,xdenda,xdendz, &
                       fpair,fpairdt,fpairdd,fpairda,fpairdz, &
                       qpair,qpairdt,qpairdd,qpairda,qpairdz

! plasma
      double precision gl2,gl2dt,gl2dd,gl2da,gl2dz,gl12,gl32,gl72,gl6, &
                       ft,ftdt,ftdd,ftda,ftdz,fl,fldt,fldd,flda,fldz, &
                       fxy,fxydt,fxydd,fxyda,fxydz

! photo
      double precision tau,taudt,cos1,cos2,cos3,cos4,cos5,sin1,sin2, &
                       sin3,sin4,sin5,last,xast, &
                       fphot,fphotdt,fphotdd,fphotda,fphotdz, &
                       qphot,qphotdt,qphotdd,qphotda,qphotdz

! brem
      double precision t8,t812,t832,t82,t83,t85,t86,t8m1,t8m2,t8m3,t8m5, &
                       t8m6, &
                       eta,etadt,etadd,etada,etadz,etam1,etam2,etam3, &
                       fbrem,fbremdt,fbremdd,fbremda,fbremdz, &
                       gbrem,gbremdt,gbremdd,gbremda,gbremdz, &
                       u,gm1,gm2,gm13,gm23,gm43,gm53,v,w,fb,gt,gb, &
                       fliq,fliqdt,fliqdd,fliqda,fliqdz, &
                       gliq,gliqdt,gliqdd,gliqda,gliqdz

! recomb
      double precision ifermi12,zfermim12,nu,nudt,nudd,nuda,nudz, &
                       nu2,nu3,bigj,bigjdt,bigjdd,bigjda,bigjdz



! numerical constants
      double precision fac1,fac2,fac3,oneth,twoth,con1,sixth,iln10
      parameter        (fac1   = 5.0d0 * pi / 3.0d0, &
                        fac2   = 10.0d0 * pi, &
                        fac3   = pi / 5.0d0, &
                        oneth  = 1.0d0/3.0d0, &
                        twoth  = 2.0d0/3.0d0, &
                        con1   = 1.0d0/5.9302d0, &
                        sixth  = 1.0d0/6.0d0, &
                        iln10  = 4.342944819032518d-1)


! theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
! xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
! change theta and xnufam if need be, and the changes will automatically
! propagate through the routine. cv and ca are the vector and axial currents.

      double precision theta,xnufam,cv,ca,cvp,cap,tfac1,tfac2,tfac3, &
                       tfac4,tfac5,tfac6
      parameter        (theta  = 0.2319d0, &
                        xnufam = 3.0d0, &
                        cv     = 0.5d0 + 2.0d0 * theta, &
                        cvp    = 1.0d0 - cv, &
                        ca     = 0.5d0, &
                        cap    = 1.0d0 - ca, &
                        tfac1  = cv*cv + ca*ca + &
                                 (xnufam-1.0d0) * (cvp*cvp+cap*cap), &
                        tfac2  = cv*cv - ca*ca + &
                                 (xnufam-1.0d0) * (cvp*cvp - cap*cap), &
                        tfac3  = tfac2/tfac1, &
                        tfac4  = 0.5d0 * tfac1, &
                        tfac5  = 0.5d0 * tfac2, &
                        tfac6  = cv*cv + 1.5d0*ca*ca + (xnufam - 1.0d0)* &
                                 (cvp*cvp + 1.5d0*cap*cap))



! initialize
      spair   = 0.0d0
      spairdt = 0.0d0
      spairdd = 0.0d0
      spairda = 0.0d0
      spairdz = 0.0d0

      splas   = 0.0d0
      splasdt = 0.0d0
      splasdd = 0.0d0
      splasda = 0.0d0
      splasdz = 0.0d0

      sphot   = 0.0d0
      sphotdt = 0.0d0
      sphotdd = 0.0d0
      sphotda = 0.0d0
      sphotdz = 0.0d0

      sbrem   = 0.0d0
      sbremdt = 0.0d0
      sbremdd = 0.0d0
      sbremda = 0.0d0
      sbremdz = 0.0d0

      sreco   = 0.0d0
      srecodt = 0.0d0
      srecodd = 0.0d0
      srecoda = 0.0d0
      srecodz = 0.0d0

      snu     = 0.0d0
      dsnudt  = 0.0d0
      dsnudd  = 0.0d0
      dsnuda  = 0.0d0
      dsnudz  = 0.0d0

      if (temp .lt. 1.0e7) return


! to avoid lots of divisions
      deni  = 1.0d0/den
      tempi = 1.0d0/temp
      abari = 1.0d0/abar
      zbari = 1.0d0/zbar


! some composition variables
      ye    = zbar*abari
      xmue  = abar*zbari




! some frequent factors
      t9     = temp * 1.0d-9
      xl     = t9 * con1
      xldt   = 1.0d-9 * con1
      xlp5   = sqrt(xl)
      xl2    = xl*xl
      xl3    = xl2*xl
      xl4    = xl3*xl
      xl5    = xl4*xl
      xl6    = xl5*xl
      xl7    = xl6*xl
      xl8    = xl7*xl
      xl9    = xl8*xl
      xlmp5  = 1.0d0/xlp5
      xlm1   = 1.0d0/xl
      xlm2   = xlm1*xlm1
      xlm3   = xlm1*xlm2
      xlm4   = xlm1*xlm3

      rm     = den*ye
      rmdd   = ye
      rmda   = -rm*abari
      rmdz   = den*abari
      rmi    = 1.0d0/rm

      a0     = rm * 1.0d-9
      a1     = a0**oneth
      zeta   = a1 * xlm1
      zetadt = -a1 * xlm2 * xldt
      a2     = oneth * a1*rmi * xlm1
      zetadd = a2 * rmdd
      zetada = a2 * rmda
      zetadz = a2 * rmdz

      zeta2 = zeta * zeta
      zeta3 = zeta2 * zeta




! pair neutrino section
! for reactions like e+ + e- => nu_e + nubar_e

! equation 2.8
      gl   = 1.0d0 - 13.04d0*xl2 +133.5d0*xl4 +1534.0d0*xl6 +918.6d0*xl8
      gldt = xldt*(-26.08d0*xl +534.0d0*xl3 +9204.0d0*xl5 +7348.8d0*xl7)

! equation 2.7

      a1     = 6.002d19 + 2.084d20*zeta + 1.872d21*zeta2
      a2     = 2.084d20 + 2.0d0*1.872d21*zeta

      if (t9 .lt. 10.0) then
       b1     = exp(-5.5924d0*zeta)
       b2     = -b1*5.5924d0
      else
       b1     = exp(-4.9924d0*zeta)
       b2     = -b1*4.9924d0
      end if

      xnum   = a1 * b1
      c      = a2*b1 + a1*b2
      xnumdt = c*zetadt
      xnumdd = c*zetadd
      xnumda = c*zetada
      xnumdz = c*zetadz

      if (t9 .lt. 10.0) then
       a1   = 9.383d-1*xlm1 - 4.141d-1*xlm2 + 5.829d-2*xlm3
       a2   = -9.383d-1*xlm2 + 2.0d0*4.141d-1*xlm3 - 3.0d0*5.829d-2*xlm4
      else
       a1   = 1.2383d0*xlm1 - 8.141d-1*xlm2
       a2   = -1.2383d0*xlm2 + 2.0d0*8.141d-1*xlm3
      end if

      b1   = 3.0d0*zeta2

      xden   = zeta3 + a1
      xdendt = b1*zetadt + a2*xldt
      xdendd = b1*zetadd
      xdenda = b1*zetada
      xdendz = b1*zetadz

      a1      = 1.0d0/xden
      fpair   = xnum*a1
      fpairdt = (xnumdt - fpair*xdendt)*a1
      fpairdd = (xnumdd - fpair*xdendd)*a1
      fpairda = (xnumda - fpair*xdenda)*a1
      fpairdz = (xnumdz - fpair*xdendz)*a1


! equation 2.6
      a1     = 10.7480d0*xl2 + 0.3967d0*xlp5 + 1.005d0
      a2     = xldt*(2.0d0*10.7480d0*xl + 0.5d0*0.3967d0*xlmp5)
      xnum   = 1.0d0/a1
      xnumdt = -xnum*xnum*a2

      a1     = 7.692d7*xl3 + 9.715d6*xlp5
      a2     = xldt*(3.0d0*7.692d7*xl2 + 0.5d0*9.715d6*xlmp5)

      c      = 1.0d0/a1
      b1     = 1.0d0 + rm*c

      xden   = b1**(-0.3d0)

      d      = -0.3d0*xden/b1
      xdendt = -d*rm*c*c*a2
      xdendd = d*rmdd*c
      xdenda = d*rmda*c
      xdendz = d*rmdz*c

      qpair   = xnum*xden
      qpairdt = xnumdt*xden + xnum*xdendt
      qpairdd = xnum*xdendd
      qpairda = xnum*xdenda
      qpairdz = xnum*xdendz



! equation 2.5
      a1    = exp(-2.0d0*xlm1)
      a2    = a1*2.0d0*xlm2*xldt

      spair   = a1*fpair
      spairdt = a2*fpair + a1*fpairdt
      spairdd = a1*fpairdd
      spairda = a1*fpairda
      spairdz = a1*fpairdz

      a1      = spair
      spair   = gl*a1
      spairdt = gl*spairdt + gldt*a1
      spairdd = gl*spairdd
      spairda = gl*spairda
      spairdz = gl*spairdz

      a1      = tfac4*(1.0d0 + tfac3 * qpair)
      a2      = tfac4*tfac3

      a3      = spair
      spair   = a1*a3
      spairdt = a1*spairdt + a2*qpairdt*a3
      spairdd = a1*spairdd + a2*qpairdd*a3
      spairda = a1*spairda + a2*qpairda*a3
      spairdz = a1*spairdz + a2*qpairdz*a3




! plasma neutrino section
! for collective reactions like gamma_plasmon => nu_e + nubar_e
! equation 4.6

      a1   = 1.019d-6*rm
      a2   = a1**twoth
      a3   = twoth*a2/a1

      b1   =  sqrt(1.0d0 + a2)
      b2   = 1.0d0/b1

      c00  = 1.0d0/(temp*temp*b1)

      gl2   = 1.1095d11 * rm * c00

      gl2dt = -2.0d0*gl2*tempi
      d     = rm*c00*b2*0.5d0*b2*a3*1.019d-6
      gl2dd = 1.1095d11 * (rmdd*c00  - d*rmdd)
      gl2da = 1.1095d11 * (rmda*c00  - d*rmda)
      gl2dz = 1.1095d11 * (rmdz*c00  - d*rmdz)


      gl    = sqrt(gl2)
      gl12  = sqrt(gl)
      gl32  = gl * gl12
      gl72  = gl2 * gl32
      gl6   = gl2 * gl2 * gl2


! equation 4.7
      ft   = 2.4d0 + 0.6d0*gl12 + 0.51d0*gl + 1.25d0*gl32
      gum  = 1.0d0/gl2
      a1   =(0.25d0*0.6d0*gl12 +0.5d0*0.51d0*gl +0.75d0*1.25d0*gl32)*gum
      ftdt = a1*gl2dt
      ftdd = a1*gl2dd
      ftda = a1*gl2da
      ftdz = a1*gl2dz


! equation 4.8
      a1   = 8.6d0*gl2 + 1.35d0*gl72
      a2   = 8.6d0 + 1.75d0*1.35d0*gl72*gum

      b1   = 225.0d0 - 17.0d0*gl + gl2
      b2   = -0.5d0*17.0d0*gl*gum + 1.0d0

      c    = 1.0d0/b1
      fl   = a1*c

      d    = (a2 - fl*b2)*c
      fldt = d*gl2dt
      fldd = d*gl2dd
      flda = d*gl2da
      fldz = d*gl2dz


! equation 4.9 and 4.10
      cc   = log10(2.0d0*rm)
      xlnt = log10(temp)

      xnum   = sixth * (17.5d0 + cc - 3.0d0*xlnt)
      xnumdt = -iln10*0.5d0*tempi
      a2     = iln10*sixth*rmi
      xnumdd = a2*rmdd
      xnumda = a2*rmda
      xnumdz = a2*rmdz

      xden   = sixth * (-24.5d0 + cc + 3.0d0*xlnt)
      xdendt = iln10*0.5d0*tempi
      xdendd = a2*rmdd
      xdenda = a2*rmda
      xdendz = a2*rmdz


! equation 4.11
      if (abs(xnum) .gt. 0.7d0  .or.  xden .lt. 0.0d0) then
       fxy   = 1.0d0
       fxydt = 0.0d0
       fxydd = 0.0d0
       fxydz = 0.0d0
       fxyda = 0.0d0

      else

       a1  = 0.39d0 - 1.25d0*xnum - 0.35d0*sin(4.5d0*xnum)
       a2  = -1.25d0 - 4.5d0*0.35d0*cos(4.5d0*xnum)

       b1  = 0.3d0 * exp(-1.0d0*(4.5d0*xnum + 0.9d0)**2)
       b2  = -b1*2.0d0*(4.5d0*xnum + 0.9d0)*4.5d0

       c   = min(0.0d0, xden - 1.6d0 + 1.25d0*xnum)
       if (c .eq. 0.0) then
        dumdt = 0.0d0
        dumdd = 0.0d0
        dumda = 0.0d0
        dumdz = 0.0d0
       else
        dumdt = xdendt + 1.25d0*xnumdt
        dumdd = xdendd + 1.25d0*xnumdd
        dumda = xdenda + 1.25d0*xnumda
        dumdz = xdendz + 1.25d0*xnumdz
       end if

       d   = 0.57d0 - 0.25d0*xnum
       a3  = c/d
       c00 = exp(-1.0d0*a3**2)

       f1  = -c00*2.0d0*a3/d
       c01 = f1*(dumdt + a3*0.25d0*xnumdt)
       c02 = f1*(dumdd + a3*0.25d0*xnumdd)
       c03 = f1*(dumda + a3*0.25d0*xnumda)
       c04 = f1*(dumdz + a3*0.25d0*xnumdz)

       fxy   = 1.05d0 + (a1 - b1)*c00
       fxydt = (a2*xnumdt -  b2*xnumdt)*c00 + (a1-b1)*c01
       fxydd = (a2*xnumdd -  b2*xnumdd)*c00 + (a1-b1)*c02
       fxyda = (a2*xnumda -  b2*xnumda)*c00 + (a1-b1)*c03
       fxydz = (a2*xnumdz -  b2*xnumdz)*c00 + (a1-b1)*c04

      end if



! equation 4.1 and 4.5
      splas   = (ft + fl) * fxy
      splasdt = (ftdt + fldt)*fxy + (ft+fl)*fxydt
      splasdd = (ftdd + fldd)*fxy + (ft+fl)*fxydd
      splasda = (ftda + flda)*fxy + (ft+fl)*fxyda
      splasdz = (ftdz + fldz)*fxy + (ft+fl)*fxydz

      a2      = exp(-gl)
      a3      = -0.5d0*a2*gl*gum

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*gl2dt*a1
      splasdd = a2*splasdd + a3*gl2dd*a1
      splasda = a2*splasda + a3*gl2da*a1
      splasdz = a2*splasdz + a3*gl2dz*a1

      a2      = gl6
      a3      = 3.0d0*gl6*gum

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*gl2dt*a1
      splasdd = a2*splasdd + a3*gl2dd*a1
      splasda = a2*splasda + a3*gl2da*a1
      splasdz = a2*splasdz + a3*gl2dz*a1


      a2      = 0.93153d0 * 3.0d21 * xl9
      a3      = 0.93153d0 * 3.0d21 * 9.0d0*xl8*xldt

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*a1
      splasdd = a2*splasdd
      splasda = a2*splasda
      splasdz = a2*splasdz




! photoneutrino process section
! for reactions like e- + gamma => e- + nu_e + nubar_e
!                    e+ + gamma => e+ + nu_e + nubar_e
! equation 3.8 for tau, equation 3.6 for cc,
! and table 2 written out for speed
      if (temp .ge. 1.0d7  .and. temp .lt. 1.0d8) then
       tau  =  log10(temp * 1.0d-7)
       cc   =  0.5654d0 + tau
       c00  =  1.008d11
       c01  =  0.0d0
       c02  =  0.0d0
       c03  =  0.0d0
       c04  =  0.0d0
       c05  =  0.0d0
       c06  =  0.0d0
       c10  =  8.156d10
       c11  =  9.728d8
       c12  = -3.806d9
       c13  = -4.384d9
       c14  = -5.774d9
       c15  = -5.249d9
       c16  = -5.153d9
       c20  =  1.067d11
       c21  = -9.782d9
       c22  = -7.193d9
       c23  = -6.936d9
       c24  = -6.893d9
       c25  = -7.041d9
       c26  = -7.193d9
       dd01 =  0.0d0
       dd02 =  0.0d0
       dd03 =  0.0d0
       dd04 =  0.0d0
       dd05 =  0.0d0
       dd11 = -1.879d10
       dd12 = -9.667d9
       dd13 = -5.602d9
       dd14 = -3.370d9
       dd15 = -1.825d9
       dd21 = -2.919d10
       dd22 = -1.185d10
       dd23 = -7.270d9
       dd24 = -4.222d9
       dd25 = -1.560d9

      else if (temp .ge. 1.0d8  .and. temp .lt. 1.0d9) then
       tau   =  log10(temp * 1.0d-8)
       cc   =  1.5654d0
       c00  =  9.889d10
       c01  = -4.524d8
       c02  = -6.088d6
       c03  =  4.269d7
       c04  =  5.172d7
       c05  =  4.910d7
       c06  =  4.388d7
       c10  =  1.813d11
       c11  = -7.556d9
       c12  = -3.304d9
       c13  = -1.031d9
       c14  = -1.764d9
       c15  = -1.851d9
       c16  = -1.928d9
       c20  =  9.750d10
       c21  =  3.484d10
       c22  =  5.199d9
       c23  = -1.695d9
       c24  = -2.865d9
       c25  = -3.395d9
       c26  = -3.418d9
       dd01 = -1.135d8
       dd02 =  1.256d8
       dd03 =  5.149d7
       dd04 =  3.436d7
       dd05 =  1.005d7
       dd11 =  1.652d9
       dd12 = -3.119d9
       dd13 = -1.839d9
       dd14 = -1.458d9
       dd15 = -8.956d8
       dd21 = -1.549d10
       dd22 = -9.338d9
       dd23 = -5.899d9
       dd24 = -3.035d9
       dd25 = -1.598d9

      else if (temp .ge. 1.0d9) then
       tau  =  log10(t9)
       cc   =  1.5654d0
       c00  =  9.581d10
       c01  =  4.107d8
       c02  =  2.305d8
       c03  =  2.236d8
       c04  =  1.580d8
       c05  =  2.165d8
       c06  =  1.721d8
       c10  =  1.459d12
       c11  =  1.314d11
       c12  = -1.169d11
       c13  = -1.765d11
       c14  = -1.867d11
       c15  = -1.983d11
       c16  = -1.896d11
       c20  =  2.424d11
       c21  = -3.669d9
       c22  = -8.691d9
       c23  = -7.967d9
       c24  = -7.932d9
       c25  = -7.987d9
       c26  = -8.333d9
       dd01 =  4.724d8
       dd02 =  2.976d8
       dd03 =  2.242d8
       dd04 =  7.937d7
       dd05 =  4.859d7
       dd11 = -7.094d11
       dd12 = -3.697d11
       dd13 = -2.189d11
       dd14 = -1.273d11
       dd15 = -5.705d10
       dd21 = -2.254d10
       dd22 = -1.551d10
       dd23 = -7.793d9
       dd24 = -4.489d9
       dd25 = -2.185d9
      end if

      taudt = iln10*tempi


! equation 3.7, compute the expensive trig functions only one time
      cos1 = cos(fac1*tau)
      cos2 = cos(fac1*2.0d0*tau)
      cos3 = cos(fac1*3.0d0*tau)
      cos4 = cos(fac1*4.0d0*tau)
      cos5 = cos(fac1*5.0d0*tau)
      last = cos(fac2*tau)

      sin1 = sin(fac1*tau)
      sin2 = sin(fac1*2.0d0*tau)
      sin3 = sin(fac1*3.0d0*tau)
      sin4 = sin(fac1*4.0d0*tau)
      sin5 = sin(fac1*5.0d0*tau)
      xast = sin(fac2*tau)

      a0 = 0.5d0*c00 &
           + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2 &
           + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4 &
           + c05*cos5 + dd05*sin5 + 0.5d0*c06*last

      f0 =  taudt*fac1*(-c01*sin1 + dd01*cos1 - c02*sin2*2.0d0 &
           + dd02*cos2*2.0d0 - c03*sin3*3.0d0 + dd03*cos3*3.0d0 &
           - c04*sin4*4.0d0 + dd04*cos4*4.0d0 &
           - c05*sin5*5.0d0 + dd05*cos5*5.0d0) &
           - 0.5d0*c06*xast*fac2*taudt

      a1 = 0.5d0*c10 &
           + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2 &
           + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4 &
           + c15*cos5 + dd15*sin5 + 0.5d0*c16*last

      f1 = taudt*fac1*(-c11*sin1 + dd11*cos1 - c12*sin2*2.0d0 &
           + dd12*cos2*2.0d0 - c13*sin3*3.0d0 + dd13*cos3*3.0d0 &
           - c14*sin4*4.0d0 + dd14*cos4*4.0d0 - c15*sin5*5.0d0 &
           + dd15*cos5*5.0d0) - 0.5d0*c16*xast*fac2*taudt

      a2 = 0.5d0*c20 &
           + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2 &
           + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4 &
           + c25*cos5 + dd25*sin5 + 0.5d0*c26*last

      f2 = taudt*fac1*(-c21*sin1 + dd21*cos1 - c22*sin2*2.0d0 &
           + dd22*cos2*2.0d0 - c23*sin3*3.0d0 + dd23*cos3*3.0d0 &
           - c24*sin4*4.0d0 + dd24*cos4*4.0d0 - c25*sin5*5.0d0 &
           + dd25*cos5*5.0d0) - 0.5d0*c26*xast*fac2*taudt

! equation 3.4
      dum   = a0 + a1*zeta + a2*zeta2
      dumdt = f0 + f1*zeta + a1*zetadt + f2*zeta2 + 2.0d0*a2*zeta*zetadt
      dumdd = a1*zetadd + 2.0d0*a2*zeta*zetadd
      dumda = a1*zetada + 2.0d0*a2*zeta*zetada
      dumdz = a1*zetadz + 2.0d0*a2*zeta*zetadz

      z      = exp(-cc*zeta)

      xnum   = dum*z
      xnumdt = dumdt*z - dum*z*cc*zetadt
      xnumdd = dumdd*z - dum*z*cc*zetadd
      xnumda = dumda*z - dum*z*cc*zetada
      xnumdz = dumdz*z - dum*z*cc*zetadz

      xden   = zeta3 + 6.290d-3*xlm1 + 7.483d-3*xlm2 + 3.061d-4*xlm3

      dum    = 3.0d0*zeta2
      xdendt = dum*zetadt - xldt*(6.290d-3*xlm2 &
               + 2.0d0*7.483d-3*xlm3 + 3.0d0*3.061d-4*xlm4)
      xdendd = dum*zetadd
      xdenda = dum*zetada
      xdendz = dum*zetadz

      dum      = 1.0d0/xden
      fphot   = xnum*dum
      fphotdt = (xnumdt - fphot*xdendt)*dum
      fphotdd = (xnumdd - fphot*xdendd)*dum
      fphotda = (xnumda - fphot*xdenda)*dum
      fphotdz = (xnumdz - fphot*xdendz)*dum


! equation 3.3
      a0     = 1.0d0 + 2.045d0 * xl
      xnum   = 0.666d0*a0**(-2.066d0)
      xnumdt = -2.066d0*xnum/a0 * 2.045d0*xldt

      dum    = 1.875d8*xl + 1.653d8*xl2 + 8.449d8*xl3 - 1.604d8*xl4
      dumdt  = xldt*(1.875d8 + 2.0d0*1.653d8*xl + 3.0d0*8.449d8*xl2 &
               - 4.0d0*1.604d8*xl3)

      z      = 1.0d0/dum
      xden   = 1.0d0 + rm*z
      xdendt =  -rm*z*z*dumdt
      xdendd =  rmdd*z
      xdenda =  rmda*z
      xdendz =  rmdz*z

      z      = 1.0d0/xden
      qphot = xnum*z
      qphotdt = (xnumdt - qphot*xdendt)*z
      dum      = -qphot*z
      qphotdd = dum*xdendd
      qphotda = dum*xdenda
      qphotdz = dum*xdendz

! equation 3.2
      sphot   = xl5 * fphot
      sphotdt = 5.0d0*xl4*xldt*fphot + xl5*fphotdt
      sphotdd = xl5*fphotdd
      sphotda = xl5*fphotda
      sphotdz = xl5*fphotdz

      a1      = sphot
      sphot   = rm*a1
      sphotdt = rm*sphotdt
      sphotdd = rm*sphotdd + rmdd*a1
      sphotda = rm*sphotda + rmda*a1
      sphotdz = rm*sphotdz + rmdz*a1

      a1      = tfac4*(1.0d0 - tfac3 * qphot)
      a2      = -tfac4*tfac3

      a3      = sphot
      sphot   = a1*a3
      sphotdt = a1*sphotdt + a2*qphotdt*a3
      sphotdd = a1*sphotdd + a2*qphotdd*a3
      sphotda = a1*sphotda + a2*qphotda*a3
      sphotdz = a1*sphotdz + a2*qphotdz*a3

      if (sphot .le. 0.0) then
       sphot   = 0.0d0
       sphotdt = 0.0d0
       sphotdd = 0.0d0
       sphotda = 0.0d0
       sphotdz = 0.0d0
      end if





! bremsstrahlung neutrino section
! for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
!                    n  + n     => n + n + nu + nubar
!                    n  + p     => n + p + nu + nubar
! equation 4.3

      den6   = den * 1.0d-6
      t8     = temp * 1.0d-8
      t812   = sqrt(t8)
      t832   = t8 * t812
      t82    = t8*t8
      t83    = t82*t8
      t85    = t82*t83
      t86    = t85*t8
      t8m1   = 1.0d0/t8
      t8m2   = t8m1*t8m1
      t8m3   = t8m2*t8m1
      t8m5   = t8m3*t8m2
      t8m6   = t8m5*t8m1


      tfermi = 5.9302d9*(sqrt(1.0d0+1.018d0*(den6*ye)**twoth)-1.0d0)

! "weak" degenerate electrons only
      if (temp .gt. 0.3d0 * tfermi) then

! equation 5.3
       dum   = 7.05d6 * t832 + 5.12d4 * t83
       dumdt = (1.5d0*7.05d6*t812 + 3.0d0*5.12d4*t82)*1.0d-8

       z     = 1.0d0/dum
       eta   = rm*z
       etadt = -rm*z*z*dumdt
       etadd = rmdd*z
       etada = rmda*z
       etadz = rmdz*z

       etam1 = 1.0d0/eta
       etam2 = etam1 * etam1
       etam3 = etam2 * etam1


! equation 5.2
       a0    = 23.5d0 + 6.83d4*t8m2 + 7.81d8*t8m5
       f0    = (-2.0d0*6.83d4*t8m3 - 5.0d0*7.81d8*t8m6)*1.0d-8
       xnum  = 1.0d0/a0

       dum   = 1.0d0 + 1.47d0*etam1 + 3.29d-2*etam2
       z     = -1.47d0*etam2 - 2.0d0*3.29d-2*etam3
       dumdt = z*etadt
       dumdd = z*etadd
       dumda = z*etada
       dumdz = z*etadz

       c00   = 1.26d0 * (1.0d0+etam1)
       z     = -1.26d0*etam2
       c01   = z*etadt
       c02   = z*etadd
       c03   = z*etada
       c04   = z*etadz

       z      = 1.0d0/dum
       xden   = c00*z
       xdendt = (c01 - xden*dumdt)*z
       xdendd = (c02 - xden*dumdd)*z
       xdenda = (c03 - xden*dumda)*z
       xdendz = (c04 - xden*dumdz)*z

       fbrem   = xnum + xden
       fbremdt = -xnum*xnum*f0 + xdendt
       fbremdd = xdendd
       fbremda = xdenda
       fbremdz = xdendz


! equation 5.9
       a0    = 230.0d0 + 6.7d5*t8m2 + 7.66d9*t8m5
       f0    = (-2.0d0*6.7d5*t8m3 - 5.0d0*7.66d9*t8m6)*1.0d-8

       z     = 1.0d0 + rm*1.0d-9
       dum   = a0*z
       dumdt = f0*z
       z     = a0*1.0d-9
       dumdd = z*rmdd
       dumda = z*rmda
       dumdz = z*rmdz

       xnum   = 1.0d0/dum
       z      = -xnum*xnum
       xnumdt = z*dumdt
       xnumdd = z*dumdd
       xnumda = z*dumda
       xnumdz = z*dumdz

       c00   = 7.75d5*t832 + 247.0d0*t8**(3.85d0)
       dd00  = (1.5d0*7.75d5*t812 + 3.85d0*247.0d0*t8**(2.85d0))*1.0d-8

       c01   = 4.07d0 + 0.0240d0 * t8**(1.4d0)
       dd01  = 1.4d0*0.0240d0*t8**(0.4d0)*1.0d-8

       c02   = 4.59d-5 * t8**(-0.110d0)
       dd02  = -0.11d0*4.59d-5 * t8**(-1.11d0)*1.0d-8

       z     = den**(0.656d0)
       dum   = c00*rmi  + c01  + c02*z
       dumdt = dd00*rmi + dd01 + dd02*z
       z     = -c00*rmi*rmi
       dumdd = z*rmdd + 0.656d0*c02*den**(-0.454d0)
       dumda = z*rmda
       dumdz = z*rmdz

       xden  = 1.0d0/dum
       z      = -xden*xden
       xdendt = z*dumdt
       xdendd = z*dumdd
       xdenda = z*dumda
       xdendz = z*dumdz

       gbrem   = xnum + xden
       gbremdt = xnumdt + xdendt
       gbremdd = xnumdd + xdendd
       gbremda = xnumda + xdenda
       gbremdz = xnumdz + xdendz


! equation 5.1
       dum    = 0.5738d0*zbar*ye*t86*den
       dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       dumdd  = 0.5738d0*zbar*ye*t86
       dumda  = -dum*abari
       dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fbrem - tfac5*gbrem
       sbrem   = dum * z
       sbremdt = dumdt*z + dum*(tfac4*fbremdt - tfac5*gbremdt)
       sbremdd = dumdd*z + dum*(tfac4*fbremdd - tfac5*gbremdd)
       sbremda = dumda*z + dum*(tfac4*fbremda - tfac5*gbremda)
       sbremdz = dumdz*z + dum*(tfac4*fbremdz - tfac5*gbremdz)




! liquid metal with c12 parameters (not too different for other elements)
! equation 5.18 and 5.16

      else
       u     = fac3 * (log10(den) - 3.0d0)
       a0    = iln10*fac3*deni

! compute the expensive trig functions of equation 5.21 only once
       cos1 = cos(u)
       cos2 = cos(2.0d0*u)
       cos3 = cos(3.0d0*u)
       cos4 = cos(4.0d0*u)
       cos5 = cos(5.0d0*u)

       sin1 = sin(u)
       sin2 = sin(2.0d0*u)
       sin3 = sin(3.0d0*u)
       sin4 = sin(4.0d0*u)
       sin5 = sin(5.0d0*u)

! equation 5.21
       fb =  0.5d0 * 0.17946d0  + 0.00945d0*u + 0.34529d0 &
             - 0.05821d0*cos1 - 0.04969d0*sin1 &
             - 0.01089d0*cos2 - 0.01584d0*sin2 &
             - 0.01147d0*cos3 - 0.00504d0*sin3 &
             - 0.00656d0*cos4 - 0.00281d0*sin4 &
             - 0.00519d0*cos5

       c00 =  a0*(0.00945d0 &
             + 0.05821d0*sin1       - 0.04969d0*cos1 &
             + 0.01089d0*sin2*2.0d0 - 0.01584d0*cos2*2.0d0 &
             + 0.01147d0*sin3*3.0d0 - 0.00504d0*cos3*3.0d0 &
             + 0.00656d0*sin4*4.0d0 - 0.00281d0*cos4*4.0d0 &
             + 0.00519d0*sin5*5.0d0)


! equation 5.22
       ft =  0.5d0 * 0.06781d0 - 0.02342d0*u + 0.24819d0 &
             - 0.00944d0*cos1 - 0.02213d0*sin1 &
             - 0.01289d0*cos2 - 0.01136d0*sin2 &
             - 0.00589d0*cos3 - 0.00467d0*sin3 &
             - 0.00404d0*cos4 - 0.00131d0*sin4 &
             - 0.00330d0*cos5

       c01 = a0*(-0.02342d0 &
             + 0.00944d0*sin1       - 0.02213d0*cos1 &
             + 0.01289d0*sin2*2.0d0 - 0.01136d0*cos2*2.0d0 &
             + 0.00589d0*sin3*3.0d0 - 0.00467d0*cos3*3.0d0 &
             + 0.00404d0*sin4*4.0d0 - 0.00131d0*cos4*4.0d0 &
             + 0.00330d0*sin5*5.0d0)


! equation 5.23
       gb =  0.5d0 * 0.00766d0 - 0.01259d0*u + 0.07917d0 &
             - 0.00710d0*cos1 + 0.02300d0*sin1 &
             - 0.00028d0*cos2 - 0.01078d0*sin2 &
             + 0.00232d0*cos3 + 0.00118d0*sin3 &
             + 0.00044d0*cos4 - 0.00089d0*sin4 &
             + 0.00158d0*cos5

       c02 = a0*(-0.01259d0 &
             + 0.00710d0*sin1       + 0.02300d0*cos1 &
             + 0.00028d0*sin2*2.0d0 - 0.01078d0*cos2*2.0d0 &
             - 0.00232d0*sin3*3.0d0 + 0.00118d0*cos3*3.0d0 &
             - 0.00044d0*sin4*4.0d0 - 0.00089d0*cos4*4.0d0 &
             - 0.00158d0*sin5*5.0d0)


! equation 5.24
       gt =  -0.5d0 * 0.00769d0  - 0.00829d0*u + 0.05211d0 &
             + 0.00356d0*cos1 + 0.01052d0*sin1 &
             - 0.00184d0*cos2 - 0.00354d0*sin2 &
             + 0.00146d0*cos3 - 0.00014d0*sin3 &
             + 0.00031d0*cos4 - 0.00018d0*sin4 &
             + 0.00069d0*cos5

       c03 = a0*(-0.00829d0 &
             - 0.00356d0*sin1       + 0.01052d0*cos1 &
             + 0.00184d0*sin2*2.0d0 - 0.00354d0*cos2*2.0d0 &
             - 0.00146d0*sin3*3.0d0 - 0.00014d0*cos3*3.0d0 &
             - 0.00031d0*sin4*4.0d0 - 0.00018d0*cos4*4.0d0 &
             - 0.00069d0*sin5*5.0d0)


       dum   = 2.275d-1 * zbar * zbar*t8m1 * (den6*abari)**oneth
       dumdt = -dum*tempi
       dumdd = oneth*dum*deni
       dumda = -oneth*dum*abari
       dumdz = 2.0d0*dum*zbari

       gm1   = 1.0d0/dum
       gm2   = gm1*gm1
       gm13  = gm1**oneth
       gm23  = gm13 * gm13
       gm43  = gm13*gm1
       gm53  = gm23*gm1


! equation 5.25 and 5.26
       v  = -0.05483d0 - 0.01946d0*gm13 + 1.86310d0*gm23 - 0.78873d0*gm1
       a0 = oneth*0.01946d0*gm43 - twoth*1.86310d0*gm53 + 0.78873d0*gm2

       w  = -0.06711d0 + 0.06859d0*gm13 + 1.74360d0*gm23 - 0.74498d0*gm1
       a1 = -oneth*0.06859d0*gm43 - twoth*1.74360d0*gm53 + 0.74498d0*gm2


! equation 5.19 and 5.20
       fliq   = v*fb + (1.0d0 - v)*ft
       fliqdt = a0*dumdt*(fb - ft)
       fliqdd = a0*dumdd*(fb - ft) + v*c00 + (1.0d0 - v)*c01
       fliqda = a0*dumda*(fb - ft)
       fliqdz = a0*dumdz*(fb - ft)

       gliq   = w*gb + (1.0d0 - w)*gt
       gliqdt = a1*dumdt*(gb - gt)
       gliqdd = a1*dumdd*(gb - gt) + w*c02 + (1.0d0 - w)*c03
       gliqda = a1*dumda*(gb - gt)
       gliqdz = a1*dumdz*(gb - gt)


! equation 5.17
       dum    = 0.5738d0*zbar*ye*t86*den
       dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       dumdd  = 0.5738d0*zbar*ye*t86
       dumda  = -dum*abari
       dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fliq - tfac5*gliq
       sbrem   = dum * z
       sbremdt = dumdt*z + dum*(tfac4*fliqdt - tfac5*gliqdt)
       sbremdd = dumdd*z + dum*(tfac4*fliqdd - tfac5*gliqdd)
       sbremda = dumda*z + dum*(tfac4*fliqda - tfac5*gliqda)
       sbremdz = dumdz*z + dum*(tfac4*fliqdz - tfac5*gliqdz)

      end if




! recombination neutrino section
! for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
! equation 6.11 solved for nu
      xnum   = 1.10520d8 * den * ye /(temp*sqrt(temp))
      xnumdt = -1.50d0*xnum*tempi
      xnumdd = xnum*deni
      xnumda = -xnum*abari
      xnumdz = xnum*zbari

! the chemical potential
      nu   = ifermi12(xnum)

! a0 is d(nu)/d(xnum)
      a0 = 1.0d0/(0.5d0*zfermim12(nu))
      nudt = a0*xnumdt
      nudd = a0*xnumdd
      nuda = a0*xnumda
      nudz = a0*xnumdz

      nu2  = nu * nu
      nu3  = nu2 * nu

! table 12
      if (nu .ge. -20.0  .and. nu .lt. 0.0) then
       a1 = 1.51d-2
       a2 = 2.42d-1
       a3 = 1.21d0
       b  = 3.71d-2
       c  = 9.06e-1
       d  = 9.28d-1
       f1 = 0.0d0
       f2 = 0.0d0
       f3 = 0.0d0
      else if (nu .ge. 0.0  .and. nu .le. 10.0) then
       a1 = 1.23d-2
       a2 = 2.66d-1
       a3 = 1.30d0
       b  = 1.17d-1
       c  = 8.97e-1
       d  = 1.77d-1
       f1 = -1.20d-2
       f2 = 2.29d-2
       f3 = -1.04d-3
      end if


! equation 6.7, 6.13 and 6.14
      if (nu .ge. -20.0  .and.  nu .le. 10.0) then

       zeta   = 1.579d5*zbar*zbar*tempi
       zetadt = -zeta*tempi
       zetadd = 0.0d0
       zetada = 0.0d0
       zetadz = 2.0d0*zeta*zbari

       c00    = 1.0d0/(1.0d0 + f1*nu + f2*nu2 + f3*nu3)
       c01    = f1 + f2*2.0d0*nu + f3*3.0d0*nu2
       dum    = zeta*c00
       dumdt  = zetadt*c00 + zeta*c01*nudt
       dumdd  = zeta*c01*nudd
       dumda  = zeta*c01*nuda
       dumdz  = zetadz*c00 + zeta*c01*nudz


       z      = 1.0d0/dum
       dd00   = dum**(-2.25)
       dd01   = dum**(-4.55)
       c00    = a1*z + a2*dd00 + a3*dd01
       c01    = -(a1*z + 2.25*a2*dd00 + 4.55*a3*dd01)*z


       z      = exp(c*nu)
       dd00   = b*z*(1.0d0 + d*dum)
       gum    = 1.0d0 + dd00
       gumdt  = dd00*c*nudt + b*z*d*dumdt
       gumdd  = dd00*c*nudd + b*z*d*dumdd
       gumda  = dd00*c*nuda + b*z*d*dumda
       gumdz  = dd00*c*nudz + b*z*d*dumdz


       z   = exp(nu)
       a1  = 1.0d0/gum

       bigj   = c00 * z * a1
       bigjdt = c01*dumdt*z*a1 + c00*z*nudt*a1 - c00*z*a1*a1 * gumdt
       bigjdd = c01*dumdd*z*a1 + c00*z*nudd*a1 - c00*z*a1*a1 * gumdd
       bigjda = c01*dumda*z*a1 + c00*z*nuda*a1 - c00*z*a1*a1 * gumda
       bigjdz = c01*dumdz*z*a1 + c00*z*nudz*a1 - c00*z*a1*a1 * gumdz


! equation 6.5
       z     = exp(zeta + nu)
       dum   = 1.0d0 + z
       a1    = 1.0d0/dum
       a2    = 1.0d0/bigj

       sreco   = tfac6 * 2.649d-18 * ye * zbar**13 * den * bigj*a1
       srecodt = sreco*(bigjdt*a2 - z*(zetadt + nudt)*a1)
       srecodd = sreco*(1.0d0*deni + bigjdd*a2 - z*(zetadd + nudd)*a1)
       srecoda = sreco*(-1.0d0*abari + bigjda*a2 - z*(zetada+nuda)*a1)
       srecodz = sreco*(14.0d0*zbari + bigjdz*a2 - z*(zetadz+nudz)*a1)

      end if


! convert from erg/cm^3/s to erg/g/s
! comment these out to duplicate the itoh et al plots

      spair   = spair*deni
      spairdt = spairdt*deni
      spairdd = spairdd*deni - spair*deni
      spairda = spairda*deni
      spairdz = spairdz*deni

      splas   = splas*deni
      splasdt = splasdt*deni
      splasdd = splasdd*deni - splas*deni
      splasda = splasda*deni
      splasdz = splasdz*deni

      sphot   = sphot*deni
      sphotdt = sphotdt*deni
      sphotdd = sphotdd*deni - sphot*deni
      sphotda = sphotda*deni
      sphotdz = sphotdz*deni

      sbrem   = sbrem*deni
      sbremdt = sbremdt*deni
      sbremdd = sbremdd*deni - sbrem*deni
      sbremda = sbremda*deni
      sbremdz = sbremdz*deni

      sreco   = sreco*deni
      srecodt = srecodt*deni
      srecodd = srecodd*deni - sreco*deni
      srecoda = srecoda*deni
      srecodz = srecodz*deni


! the total neutrino loss rate
      snu    =  splas + spair + sphot + sbrem + sreco
      dsnudt =  splasdt + spairdt + sphotdt + sbremdt + srecodt
      dsnudd =  splasdd + spairdd + sphotdd + sbremdd + srecodd
      dsnuda =  splasda + spairda + sphotda + sbremda + srecoda
      dsnudz =  splasdz + spairdz + sphotdz + sbremdz + srecodz

      return
      end






      double precision function ifermi12(f)
      include 'implno.dek'

! this routine applies a rational function expansion to get the inverse
! fermi-dirac integral of order 1/2 when it is equal to f.
! maximum error is 4.19d-9.   reference: antia apjs 84,101 1993

! declare
      integer          i,m1,k1,m2,k2
      double precision f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff


! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
      data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3, &
           6.610132843877d2,   3.818838129486d1, &
           1.0d0/
      data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3, &
           9.130355392717d1,  -1.670718177489d0/
      data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2, &
                          -4.262314235106d-1,  4.997559426872d-1, &
                          -1.285579118012d0,  -3.930805454272d-1, &
           1.0d0/
      data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2, &
                          -3.299466243260d-1,  4.077841975923d-1, &
                          -1.145531476975d0,  -6.067091689181d-2/


      if (f .lt. 4.0d0) then
       rn  = f + a1(m1)
       do i=m1-1,1,-1
        rn  = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi12 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi12 = rn/(den*ff)
      end if
      return
      end






      double precision function zfermim12(x)
      include 'implno.dek'

! this routine applies a rational function expansion to get the fermi-dirac
! integral of order -1/2 evaluated at x. maximum error is 1.23d-12.
! reference: antia apjs 84,101 1993

! declare
      integer          i,m1,k1,m2,k2
      double precision x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /-0.5d0, 7, 7, 11, 11/
      data  (a1(i),i=1,8)/ 1.71446374704454d7,    3.88148302324068d7, &
                           3.16743385304962d7,    1.14587609192151d7, &
                           1.83696370756153d6,    1.14980998186874d5, &
                           1.98276889924768d3,    1.0d0/
      data  (b1(i),i=1,8)/ 9.67282587452899d6,    2.87386436731785d7, &
                           3.26070130734158d7,    1.77657027846367d7, &
                           4.81648022267831d6,    6.13709569333207d5, &
                           3.13595854332114d4,    4.35061725080755d2/
      data (a2(i),i=1,12)/-4.46620341924942d-15, -1.58654991146236d-12, &
                          -4.44467627042232d-10, -6.84738791621745d-8, &
                          -6.64932238528105d-6,  -3.69976170193942d-4, &
                          -1.12295393687006d-2,  -1.60926102124442d-1, &
                          -8.52408612877447d-1,  -7.45519953763928d-1, &
                           2.98435207466372d0,    1.0d0/
      data (b2(i),i=1,12)/-2.23310170962369d-15, -7.94193282071464d-13, &
                          -2.22564376956228d-10, -3.43299431079845d-8, &
                          -3.33919612678907d-6,  -1.86432212187088d-4, &
                          -5.69764436880529d-3,  -8.34904593067194d-2, &
                          -4.78770844009440d-1,  -4.99759250374148d-1, &
                           1.86795964993052d0,    4.16485970495288d-1/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermim12 = xx * rn/den
!
      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermim12 = sqrt(x)*rn/den
      end if
      return
      end




!---------------------------------------------------------------------
! reaction rates


      subroutine tfactors(temp)
      include 'implno.dek'
      include 'tfactors.dek'

! sets various popular temperature factors into common block
! this routine must be called before any of the rates are called

! declare the pass
      double precision temp

! all these are in common block

      t9    = temp * 1.0d-9
      t92   = t9*t9
      t93   = t9*t92
      t94   = t9*t93
      t95   = t9*t94
      t96   = t9*t95

      t912  = sqrt(t9)
      t932  = t9*t912
      t952  = t9*t932
      t972  = t9*t952

      t913  = t9**oneth
      t923  = t913*t913
      t943  = t9*t913
      t953  = t9*t923
      t973  = t953*t923
      t9113 = t973*t943

      t914  = t9**(0.25d0)
      t934  = t914*t914*t914
      t954  = t9*t914
      t974  = t9*t934

      t915  = t9**onefif
      t935  = t915*t915*t915
      t945  = t915 * t935
      t965  = t9 * t915

      t916  = t9**onesix
      t976  = t9 * t916
      t9i76 = 1.0d0/t976

      t917  = t9**onesev
      t927  = t917*t917
      t947  = t927*t927

      t918  = sqrt(t914)
      t938  = t918*t918*t918
      t958  = t938*t918*t918

      t9i   = 1.0d0/t9
      t9i2  = t9i*t9i
      t9i3  = t9i2*t9i

      t9i12 = 1.0d0/t912
      t9i32 = t9i*t9i12
      t9i52 = t9i*t9i32
      t9i72 = t9i*t9i52

      t9i13 = 1.0d0/t913
      t9i23 = t9i13*t9i13
      t9i43 = t9i*t9i13
      t9i53 = t9i*t9i23

      t9i14 = 1.0d0/t914
      t9i34 = t9i14*t9i14*t9i14
      t9i54 = t9i*t9i14

      t9i15 = 1.0d0/t915
      t9i35 = t9i15*t9i15*t9i15
      t9i45 = t9i15 * t9i35
      t9i65 = t9i*t9i15

      t9i17 = 1.0d0/t917
      t9i27 = t9i17*t9i17
      t9i47 = t9i27*t9i27

      t9i18 = 1.0d0/t918
      t9i38 = t9i18*t9i18*t9i18
      t9i58 = t9i38*t9i18*t9i18

      return
      end





      subroutine rate_c12ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,f1,df1,f2,df2, &
                       zz,q1
      parameter        (q1 = 1.0d0/12.222016d0)


! c12(a,g)o16
      aa   = 1.0d0 + 0.0489d0*t9i23
      daa  = -twoth*0.0489d0*t9i53

      bb   = t92*aa*aa
      dbb  = 2.0d0*(bb*t9i + t92*aa*daa)

      cc   = exp(-32.120d0*t9i13 - t92*q1)
      dcc  = cc * (oneth*32.120d0*t9i43 - 2.0d0*t9*q1)

      dd   = 1.0d0 + 0.2654d0*t9i23
      ddd  = -twoth*0.2654d0*t9i53

      ee   = t92*dd*dd
      dee  = 2.0d0*(ee*t9i + t92*dd*ddd)

      ff   = exp(-32.120d0*t9i13)
      dff  = ff * oneth*32.120d0*t9i43

      gg   = 1.25d3 * t9i32 * exp(-27.499*t9i)
      dgg  = gg*(-1.5d0*t9i + 27.499*t9i2)

      hh   = 1.43d-2 * t95 * exp(-15.541*t9i)
      dhh  = hh*(5.0d0*t9i + 15.541*t9i2)

      zz   = 1.0d0/bb
      f1   = cc*zz
      df1  = (dcc - f1*dbb)*zz

      zz   = 1.0d0/ee
      f2   = ff*zz
      df2  = (dff - f2*dee)*zz

      term    = 1.04d8*f1  + 1.76d8*f2 + gg + hh
      dtermdt = 1.04d8*df1 + 1.76d8*df2 + dgg + dhh


! 1.7 times cf88 value
      term     = 1.7d0 * term
      dtermdt  = 1.7d0 * dtermdt

      fr    = term * den
      dfrdt = dtermdt * den * 1.0d-9
      dfrdd = term

      rev    = 5.13d10 * t932 * exp(-83.111*t9i)
      drevdt = rev*(1.5d0*t9i + 83.111*t9i2)

      rr     = rev * term
      drrdt  = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd  = 0.0d0

      return
      end






      subroutine rate_tripalf(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,r2abe,dr2abedt,rbeac, &
                       drbeacdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,xx,dxx,yy,dyy,zz,dzz,uu,vv,f1,df1,rc28, &
                       q1,q2
      parameter        (rc28   = 0.1d0, &
                        q1     = 1.0d0/0.009604d0, &
                        q2     = 1.0d0/0.055225d0)



! triple alfa to c12
! this is a(a,g)be8
      aa    = 7.40d+05 * t9i32 * exp(-1.0663*t9i)
      daa   = aa*(-1.5d0*t9i  + 1.0663*t9i2)

      bb    = 4.164d+09 * t9i23 * exp(-13.49*t9i13 - t92*q1)
      dbb   = bb*(-twoth*t9i + oneth*13.49*t9i43 - 2.0d0*t9*q1)

      cc    = 1.0d0 + 0.031*t913 + 8.009*t923 + 1.732*t9 &
              + 49.883*t943 + 27.426*t953
      dcc   = oneth*0.031*t9i23 + twoth*8.009*t9i13 + 1.732 &
              + fourth*49.883*t913 + fiveth*27.426*t923

      r2abe    = aa + bb * cc
      dr2abedt = daa + dbb*cc + bb*dcc


! this is be8(a,g)c12
      dd    = 130.0d0 * t9i32 * exp(-3.3364*t9i)
      ddd   = dd*(-1.5d0*t9i + 3.3364*t9i2)

      ee    = 2.510d+07 * t9i23 * exp(-23.57*t9i13 - t92*q2)
      dee   = ee*(-twoth*t9i + oneth*23.57*t9i43 - 2.0d0*t9*q2)

      ff    = 1.0d0 + 0.018*t913 + 5.249*t923 + 0.650*t9 + &
              19.176*t943 + 6.034*t953
      dff   = oneth*0.018*t9i23 + twoth*5.249*t9i13 + 0.650 &
              + fourth*19.176*t913 + fiveth*6.034*t923

      rbeac    = dd + ee * ff
      drbeacdt = ddd + dee * ff + ee * dff


! a factor
      xx    = rc28 * 1.35d-07 * t9i32 * exp(-24.811*t9i)
      dxx   = xx*(-1.5d0*t9i + 24.811*t9i2)


! high temperature rate
      if (t9.gt.0.08) then
       term    = 2.90d-16 * r2abe * rbeac + xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac &
                 + 2.90d-16 * r2abe * drbeacdt &
                 + dxx


! low temperature rate
      else
       uu   = 0.8d0*exp(-(0.025*t9i)**3.263)
       yy   = 0.01 + 0.2d0 + uu
       dyy  = uu * 3.263*(0.025*t9i)**2.263 * (0.025*t9i2)
       vv   = 4.0d0*exp(-(t9/0.025)**9.227)
       zz   = 1.0d0 + vv
       dzz  = vv * 9.227*(t9/0.025)**8.227 * 40.0d0
       aa   = 1.0d0/zz
       f1   = yy * aa
       df1  = (dyy - f1*dzz)*aa
       term = 2.90d-16 * r2abe * rbeac * f1 +  xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac * f1 &
                 + 2.90d-16 * r2abe * drbeacdt * f1 &
                 + 2.90d-16 * r2abe * rbeac * df1 &
                 + dxx
      end if


! rates
!      term    = 1.2d0 * term
!      dtermdt = 1.2d0 * term

      fr    = term * den * den
      dfrdt = dtermdt * den * den * 1.0d-9
      dfrdd = 2.0d0 * term * den

      rev    = 2.00d+20*t93*exp(-84.424*t9i)
      drevdt = rev*(3.0d0*t9i + 84.424*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_c12c12(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       aa,zz


! c12 + c12 reaction
      aa      = 1.0d0 + 0.0396*t9
      zz      = 1.0d0/aa

      t9a     = t9*zz
      dt9a    = (1.0d0 -  t9a*0.0396)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      term    = 4.27d+26 * t9a56 * t9i32 * &
                exp(-84.165/t9a13 - 2.12d-03*t93)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 84.165/t9a13**2*dt9a13 - 6.36d-3*t92)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_c12o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a23,dt9a23, &
                       t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,zz


! c12 + o16 reaction; see cf88 references 47-4

      if (t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.055)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0d0 * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa      = exp(-0.18*t9a*t9a)
       daa     = -aa * 0.36 * t9a * dt9a

       bb      = 1.06d-03*exp(2.562*t9a23)
       dbb     = bb * 2.562 * dt9a23

       cc      = aa + bb
       dcc     = daa + dbb

       zz      = 1.0d0/cc
       term    = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) * zz
       dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                       + 106.594/t9a23*dt9a13 - zz*dcc)

      else
!       term    = 2.6288035d-29
       term    = 0.0d0
       dtermdt = 0.0d0
      endif


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_o16o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt


! o16 + o16
      term  = 7.10d36 * t9i23 * &
              exp(-135.93 * t9i13 - 0.629*t923 &
                   - 0.445*t943 + 0.0103*t9*t9)

      dtermdt = -twoth*term*t9i &
                + term * (oneth*135.93*t9i43 - twoth*0.629*t9i13 &
                          - fourth*0.445*t913 + 0.0206*t9)


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_o16ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                       cc,dcc,term2,dterm2,rev,drevdt,q1
      parameter        (q1 = 1.0d0/2.515396d0)



! o16(a,g)ne20
      term1   = 9.37d9 * t9i23 * exp(-39.757*t9i13 - t92*q1)
      dterm1  = term1*(-twoth*t9i + oneth*39.757*t9i43 - 2.0d0*t9*q1)

      aa      = 62.1 * t9i32 * exp(-10.297*t9i)
      daa     = aa*(-1.5d0*t9i + 10.297*t9i2)

      bb      = 538.0d0 * t9i32 * exp(-12.226*t9i)
      dbb     = bb*(-1.5d0*t9i + 12.226*t9i2)

      cc      = 13.0d0 * t92 * exp(-20.093*t9i)
      dcc     = cc*(2.0d0*t9i + 20.093*t9i2)

      term2   = aa + bb + cc
      dterm2  = daa + dbb + dcc

      term    = term1 + term2
      dtermdt = dterm1 + dterm2


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.65d+10*t932*exp(-54.937*t9i)
      drevdt   = rev*(1.5d0*t9i + 54.937*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne20ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                       term2,dterm2,term3,dterm3,rev,drevdt,zz,rc102,q1
      parameter        (rc102 = 0.1d0, &
                        q1    = 1.0d0/4.923961d0)



! ne20(a,g)mg24

      aa   = 4.11d+11 * t9i23 * exp(-46.766*t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*46.766*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.009*t913 + 0.882*t923 + 0.055*t9 &
             + 0.749*t943 + 0.119*t953
      dbb  = oneth*0.009*t9i23 + twoth*0.882*t9i13 + 0.055 &
             + fourth*0.749*t913 + fiveth*0.119*t923

      term1  = aa * bb
      dterm1 = daa * bb + aa * dbb


      aa   = 5.27d+03 * t9i32 * exp(-15.869*t9i)
      daa  = aa*(-1.5d0*t9i + 15.869*t9i2)

      bb   = 6.51d+03 * t912 * exp(-16.223*t9i)
      dbb  = bb*(0.5d0*t9i + 16.223*t9i2)

      term2  = aa + bb
      dterm2 = daa + dbb


      aa   = 42.1 * t9i32 * exp(-9.115*t9i)
      daa  = aa*(-1.5d0*t9i + 9.115*t9i2)

      bb   =  32.0 * t9i23 * exp(-9.383*t9i)
      dbb  = bb*(-twoth*t9i + 9.383*t9i2)

      term3  = rc102 * (aa + bb)
      dterm3 = rc102 * (daa + dbb)


      aa  = 5.0d0*exp(-18.960*t9i)
      daa = aa*18.960*t9i2

      bb  = 1.0d0 + aa
      dbb = daa

      zz      = 1.0d0/bb
      term    = (term1 + term2 + term3)*zz
      dtermdt = ((dterm1 + dterm2 + dterm3) - term*dbb)*zz


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.01d+10 * t932 * exp(-108.059*t9i)
      drevdt   = rev*(1.5d0*t9i + 108.059*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg24ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,hh,hhi,rev,drevdt,rc121
      parameter        (rc121 = 0.1d0)



! 24mg(a,g)28si

      aa    = 4.78d+01 * t9i32 * exp(-13.506*t9i)
      daa   = aa*(-1.5d0*t9i + 13.506*t9i2)

      bb    =  2.38d+03 * t9i32 * exp(-15.218*t9i)
      dbb   = bb*(-1.5d0*t9i + 15.218*t9i2)

      cc    = 2.47d+02 * t932 * exp(-15.147*t9i)
      dcc   = cc*(1.5d0*t9i + 15.147*t9i2)

      dd    = rc121 * 1.72d-09 * t9i32 * exp(-5.028*t9i)
      ddd   = dd*(-1.5d0*t9i + 5.028*t9i2)

      ee    = rc121* 1.25d-03 * t9i32 * exp(-7.929*t9i)
      dee   = ee*(-1.5d0*t9i + 7.929*t9i2)

      ff    = rc121 * 2.43d+01 * t9i * exp(-11.523*t9i)
      dff   = ff*(-t9i + 11.523*t9i2)

      gg    = 5.0d0*exp(-15.882*t9i)
      dgg   = gg*15.882*t9i2

      hh    = 1.0d0 + gg
      hhi   = 1.0d0/hh

      term    = (aa + bb + cc + dd + ee + ff) * hhi
      dtermdt = (daa + dbb + dcc + ddd + dee + dff - term*dgg) * hhi


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.27d+10 * t932 * exp(-115.862*t9i)
      drevdt   = rev*(1.5d0*t9i + 115.862*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_mg24ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,term1,dterm1,term2,dterm2, &
                       rev,drevdt,rc148,q1
      parameter        (rc148 = 0.1d0, &
                        q1    = 1.0d0/0.024649d0)



! 24mg(a,p)al27
      aa     = 1.10d+08 * t9i23 * exp(-23.261*t9i13 - t92*q1)
      daa    = -twoth*aa*t9i + aa*(23.261*t9i43 - 2.0d0*t9*q1)

      bb     =  1.0d0 + 0.018*t913 + 12.85*t923 + 1.61*t9 &
               + 89.87*t943 + 28.66*t953
      dbb    = oneth*0.018*t9i23 + twoth*12.85*t9i13 + 1.61 &
                + fourth*89.87*t913 + fiveth*28.66*t923

      term1  = aa * bb
      dterm1 = daa * bb + aa * dbb

      aa     = 129.0d0 * t9i32 * exp(-2.517*t9i)
      daa    = -1.5d0*aa*t9i + aa*2.517*t9i2

      bb     = 5660.0d0 * t972 * exp(-3.421*t9i)
      dbb    = 3.5d0*bb*t9i +  bb*3.421*t9i2

      cc     = rc148 * 3.89d-08 * t9i32 * exp(-0.853*t9i)
      dcc    = -1.5d0*cc*t9i + cc*0.853*t9i2

      dd     = rc148 * 8.18d-09 * t9i32 * exp(-1.001*t9i)
      ddd    = -1.5d0*dd*t9i + dd*1.001*t9i2

      term2  = aa + bb + cc + dd
      dterm2 = daa + dbb + dcc + ddd

      ee     = oneth*exp(-9.792*t9i)
      dee    = ee*9.792*t9i2

      ff     =  twoth * exp(-11.773*t9i)
      dff    = ff*11.773*t9i2

      gg     = 1.0d0 + ee + ff
      dgg    = dee + dff

      term    = (term1 + term2)/gg
      dtermdt = ((dterm1 + dterm2) - term*dgg)/gg


! the rates
      rev      = 1.81 * exp(-18.572*t9i)
      drevdt   = rev*18.572*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end






      subroutine rate_al27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg


! al27(p,g)si28
! champagne 1996

      aa  = 1.32d+09 * t9i23 * exp(-23.26*t9i13)
      daa = aa*(-twoth*t9i + oneth*23.26*t9i43)

      bb  = 3.22d-10 * t9i32 * exp(-0.836*t9i)*0.17
      dbb = bb*(-1.5d0*t9i + 0.836*t9i2)

      cc  = 1.74d+00 * t9i32 * exp(-2.269*t9i)
      dcc = cc*(-1.5d0*t9i + 2.269*t9i2)

      dd  = 9.92d+00 * t9i32 * exp(-2.492*t9i)
      ddd = dd*(-1.5d0*t9i + 2.492*t9i2)

      ee  = 4.29d+01 * t9i32 * exp(-3.273*t9i)
      dee = ee*(-1.5d0*t9i + 3.273*t9i2)

      ff  = 1.34d+02 * t9i32 * exp(-3.654*t9i)
      dff = ff*(-1.5d0*t9i + 3.654*t9i2)

      gg  = 1.77d+04 * (t9**0.53) * exp(-4.588*t9i)
      dgg = gg*(0.53*t9i + 4.588*t9i2)

      term    = aa + bb + cc + dd + ee + ff + gg
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.13d+11 * t932 * exp(-134.434*t9i)
      drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si28ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! si28(a,g)s32
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 6.340d-2*z + 2.541d-3*z2 - 2.900d-4*z3
      if (z .eq. 10.0) then
       daa = 0
      else
       daa   = 6.340d-2 + 2.0d0*2.541d-3*t9 - 3.0d0*2.900d-4*t92
      end if

      term    = 4.82d+22 * t9i23 * exp(-61.015 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 61.015*t9i13*(oneth*t9i*aa - daa))

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.461d+10 * t932 * exp(-80.643*t9i)
      drevdt   = rev*(1.5d0*t9i + 80.643*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si28ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! si28(a,p)p31
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 2.798d-3*z + 2.763d-3*z2 - 2.341d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 2.798d-3 + 2.0d0*2.763d-3*t9 - 3.0d0*2.341d-4*t92
      end if

      term    = 4.16d+13 * t9i23 * exp(-25.631 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*25.631*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 0.5825d0 * exp(-22.224*t9i)
      drevdt   = rev*22.224*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_p31pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! p31(p,g)s32
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.928d-1*z - 1.540d-2*z2 + 6.444d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.928d-1 - 2.0d0*1.540d-2*t9 + 3.0d0*6.444d-4*t92
      end if

      term    = 1.08d+16 * t9i23 * exp(-27.042 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 27.042*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.764d+10 * t932 * exp(-102.865*t9i)
      drevdt   = rev*(1.5d0*t9i + 102.865*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_s32ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! s32(a,g)ar36
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 4.913d-2*z + 4.637d-3*z2 - 4.067d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 4.913d-2 + 2.0d0*4.637d-3*t9 - 3.0d0*4.067d-4*t92
      end if

      term    = 1.16d+24 * t9i23 * exp(-66.690 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 66.690*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.616d+10 * t932 * exp(-77.080*t9i)
      drevdt   = rev*(1.5d0*t9i + 77.080*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_s32ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! s32(a,p)cl35
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.041d-1*z - 1.368d-2*z2 + 6.969d-4*z3
      if (z .eq. 10) then
       daa = 0.0d0
      else
       daa   = 1.041d-1 - 2.0d0*1.368d-2*t9 + 3.0d0*6.969d-4*t92
      end if

      term    = 1.27d+16 * t9i23 * exp(-31.044 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*31.044*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.144 * exp(-21.643*t9i)
      drevdt   = rev*21.643*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_cl35pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt


! cl35(p,g)ar36
      aa    = 1.0d0 + 1.761d-1*t9 - 1.322d-2*t92 + 5.245d-4*t93
      daa   = 1.761d-1 - 2.0d0*1.322d-2*t9 + 3.0d0*5.245d-4*t92


      term    =  4.48d+16 * t9i23 * exp(-29.483 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 29.483*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.568d+10*t932*exp(-98.722*t9i)
      drevdt   = rev*(1.5d0*t9i + 98.722*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ar36ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ar36(a,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.458d-1*z - 1.069d-2*z2 + 3.790d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.458d-1 - 2.0d0*1.069d-2*t9 + 3.0d0*3.790d-4*t92
      end if

      term    = 2.81d+30 * t9i23 * exp(-78.271 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 78.271*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.740d+10 * t932 * exp(-81.711*t9i)
      drevdt   = rev*(1.5d0*t9i + 81.711*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ar36ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ar36(a,p)k39
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 4.826d-3*z - 5.534d-3*z2 + 4.021d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 4.826d-3 - 2.0d0*5.534d-3*t9 + 3.0d0*4.021d-4*t92
      end if

      term    = 2.76d+13 * t9i23 * exp(-34.922 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*34.922*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.128*exp(-14.959*t9i)
      drevdt   = rev*14.959*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_k39pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! k39(p,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.622d-1*z - 1.119d-2*z2 + 3.910d-4*z3
      if (z .eq. 10) then
       daa = 0.0d0
      else
       daa   = 1.622d-1 - 2.0d0*1.119d-2*t9 + 3.0d0*3.910d-4*t92
      end if

      term    = 4.09d+16 * t9i23 * exp(-31.727 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 31.727*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.600d+10 * t932 * exp(-96.657*t9i)
      drevdt   = rev*(1.5d0*t9i + 96.657*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ca40ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ca40(a,g)ti44
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.650d-2*z + 5.973d-3*z2 - 3.889d-04*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.650d-2 + 2.0d0*5.973d-3*t9 - 3.0d0*3.889d-4*t92
      end if

      term    = 4.66d+24 * t9i23 * exp(-76.435 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 76.435*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.843d+10 * t932 * exp(-59.510*t9i)
      drevdt   = rev*(1.5d0*t9i + 59.510*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ca40ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ca40(a,p)sc43
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 - 1.206d-2*z + 7.753d-3*z2 - 5.071d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = -1.206d-2 + 2.0d0*7.753d-3*t9 - 3.0d0*5.071d-4*t92
      end if

      term    = 4.54d+14 * t9i23 * exp(-32.177 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*32.177*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 2.229 * exp(-40.966*t9i)
      drevdt   = rev*40.966*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_sc43pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! sc43(p,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.023d-1*z - 2.242d-3*z2 - 5.463d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.023d-1 - 2.0d0*2.242d-3*t9 - 3.0d0*5.463d-5*t92
      end if

      term    = 3.85d+16 * t9i23 * exp(-33.234 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 33.234*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.525d+11 * t932 * exp(-100.475*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.475*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ti44ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ti44(a,g)cr48
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.066d-1*z - 1.102d-2*z2 + 5.324d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.066d-1 - 2.0d0*1.102d-2*t9 + 3.0d0*5.324d-4*t92
      end if

      term    = 1.37d+26 * t9i23 * exp(-81.227 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 81.227*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.928d+10*t932*exp(-89.289*t9i)
      drevdt   = rev*(1.5d0*t9i + 89.289*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ti44ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ti44(a,p)v47
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 2.655d-2*z - 3.947d-3*z2 + 2.522d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 2.655d-2 - 2.0d0*3.947d-3*t9 + 3.0d0*2.522d-4*t92
      end if

      term    = 6.54d+20 * t9i23 * exp(-66.678 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*66.678*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.104 * exp(-4.723*t9i)
      drevdt   = rev*4.723*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_v47pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! v47(p,g)cr48
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.979d-2*z - 2.269d-3*z2 - 6.662d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.979d-2 - 2.0d0*2.269d-3*t9 - 3.0d0*6.662d-5*t92
      end if

      term    = 2.05d+17 * t9i23 * exp(-35.568 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 35.568*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.649d+10*t932*exp(-93.999*t9i)
      drevdt   = rev*(1.5d0*t9i + 93.999*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_cr48ag(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! cr48(a,g)fe52
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 6.325d-2*z - 5.671d-3*z2 + 2.848d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 6.325d-2 - 2.0d0*5.671d-3*t9 + 3.0d0*2.848d-4*t92
      end if

      term    = 1.04d+23 * t9i23 * exp(-81.420 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 81.420*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.001d+10 * t932 * exp(-92.177*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.177*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_cr48ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! cr48(a,p)mn51
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.384d-2*z + 1.081d-3*z2 - 5.933d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.384d-2 + 2.0d0*1.081d-3*t9 - 3.0d0*5.933d-5*t92
      end if

      term    = 1.83d+26 * t9i23 * exp(-86.741 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*86.741*t9i13*(oneth*t9i*aa - daa)


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.6087*exp(-6.510*t9i)
      drevdt   = rev*6.510*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_mn51pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! mn51(p,g)fe52
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 8.922d-2*z - 1.256d-3*z2 - 9.453d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 8.922d-2 - 2.0d0*1.256d-3*t9 - 3.0d0*9.453d-5*t92
      end if

      term    = 3.77d+17 * t9i23 * exp(-37.516 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 37.516*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.150d+11*t932*exp(-85.667*t9i)
      drevdt   = rev*(1.5d0*t9i + 85.667*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_fe52ag(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! fe52(a,g)ni56
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 7.846d-2*z - 7.430d-3*z2 + 3.723d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 7.846d-2 - 2.0d0*7.430d-3*t9 + 3.0d0*3.723d-4*t92
      end if

      term    = 1.05d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 91.674*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.064d+10*t932*exp(-92.850*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.850*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_fe52ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! fe52(a,p)co55
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.367d-2*z + 7.428d-4*z2 - 3.050d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.367d-2 + 2.0d0*7.428d-4*t9 - 3.0d0*3.050d-5*t92
      end if

      term    = 1.30d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*91.674*t9i13*(oneth*t9i*aa - daa)


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.4597*exp(-9.470*t9i)
      drevdt   = rev*9.470*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_co55pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! co55(p,g)ni56
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.894d-2*z - 3.131d-3*z2 - 2.160d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.894d-2 - 2.0d0*3.131d-3*t9 - 3.0d0*2.160d-5*t92
      end if

      term    = 1.21d+18 * t9i23 * exp(-39.604 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 39.604*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.537d+11*t932*exp(-83.382*t9i)
      drevdt   = rev*(1.5d0*t9i + 83.382*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_fe52ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,tq2


! fe52(n,g)fe53
      tq2     = t9 - 0.348d0
      term    = 9.604d+05 * exp(-0.0626*tq2)
      dtermdt = -term*0.0626

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.43d+09 * t932 * exp(-123.951*t9i)
      drevdt   = rev*(1.5d0*t9i + 123.951*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_fe53ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,tq1,tq10,dtq10,tq2


! fe53(n,g)fe54
      tq1   = t9/0.348
      tq10  = tq1**0.10
      dtq10 = 0.1d0*tq10/(0.348*tq1)
      tq2   = t9 - 0.348d0

      term    = 1.817d+06 * tq10 * exp(-0.06319*tq2)
      dtermdt = term/tq10*dtq10 - term*0.06319

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.56d+11 * t932 * exp(-155.284*t9i)
      drevdt   = rev*(1.5d0*t9i + 155.284*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_fe54pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,z,z2,z3


! fe54(p,g)co55
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.593d-2*z - 3.445d-3*z2 + 8.594d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.593d-2 - 2.0d0*3.445d-3*t9 + 3.0d0*8.594d-5*t92
      end if

      term    = 4.51d+17 * t9i23 * exp(-38.483 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 38.483*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.400d+09 * t932 * exp(-58.605*t9i)
      drevdt   = rev*(1.5d0*t9i + 58.605*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





!---------------------------------------------------------------------
! this file contains routines that sort, search and select parts of arrays:
!
! index and rank makers:
! routine indexx constructs a sort index for a real array



      subroutine indexx(n,arr,indx)
      include 'implno.dek'
!
! indexes an array arr(1:n). that is it outputs the array indx(1:n) such
! that arr(indx(j)) is in ascending order for j=1...n. the input quantities
! are not changed.
!
! declare
      integer          n,indx(n),m,nstack
      parameter        (m=7, nstack = 50)
      integer          i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      double precision arr(n),a
!
! initialize
      do 11 j=1,n
       indx(j) = j
11    continue
      jstack = 0
      l      = 1
      ir     = n
!
! insertion sort when subbarray small enough
1     if (ir - l .lt. m) then
       do 13 j=l+1,ir
        indxt = indx(j)
        a     = arr(indxt)
        do 12 i=j-1,l,-1
         if (arr(indx(i)) .le. a) go to 2
         indx(i+1) = indx(i)
12      continue
        i = l - 1
2       indx(i+1) = indxt
13     continue
!
! pop stack and begin a new round of partitioning
       if (jstack .eq. 0) return
       ir     = istack(jstack)
       l      = istack(jstack-1)
       jstack = jstack - 2
!
! choose median of left, center and right elements as partitioning element
! also rearrange so that a(l+1) < a(l) < a(ir)
      else
       k         = (l + ir)/2
       itemp     = indx(k)
       indx(k)   = indx(l+1)
       indx(l+1) = itemp

       if (arr(indx(l)) .gt. arr(indx(ir))) then
        itemp    = indx(l)
        indx(l)  = indx(ir)
        indx(ir) = itemp
       end if


       if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
       endif
       if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
       endif

!
! initialize pointers for partitioning
       i     = l + 1
       j     = ir
       indxt = indx(l+1)
       a     = arr(indxt)
3      continue
       i = i + 1
       if (arr(indx(i)) .lt. a) go to 3
4      continue
       j = j - 1
       if (arr(indx(j)) .gt. a) go to 4
       if (j .lt. i) go to 5
       itemp   = indx(i)
       indx(i) = indx(j)
       indx(j) = itemp
       go to 3
!
5      indx(l+1) = indx(j)
       indx(j)   = indxt
       jstack    = jstack + 2
!
! push pointers to larger subarray on stack
       if (jstack .gt. nstack) stop 'jstack > nstack in routine indexx'
       if (ir - i + 1  .ge.  j - l) then
        istack(jstack)   = ir
        istack(jstack-1) = i
        ir               = j - 1
       else
        istack(jstack)   = j-1
        istack(jstack-1) = l
        l                = i
       end if
      end if
      go to 1
      end
!---------------------------------------------------------------------




      subroutine net_pzextr(iest,xest,yest,yz,dy,nv)
      include 'implno.dek'

! use polynomial extrapolation to evaluate nv functions at x=0 by fitting
! a polynomial to a sequence of estimates with progressively smaller values
! x=xest, and corresponding function vectors yest(1:nv). the call is number
! iest in the sequence of calls. extrapolated function values are output as
! yz(1:nv), and their estimated error is output as dy(1:nv)


! declare
      integer          iest,nv,j,k1,nmax,imax
      parameter        (nmax=3500, imax=13)
      double precision xest,dy(nv),yest(nv),yz(nv),delta,f1,f2,q, &
                       d(nmax),qcol(nmax,imax),x(imax)


! sanity checks

      if (iest .gt. imax) stop 'iest > imax in net_pzextr'
      if (nv .gt. nmax) stop 'nv > nmax in net_pzextr'


! save current independent variables
      x(iest) = xest
      do j=1,nv
       dy(j) = yest(j)
       yz(j) = yest(j)
      enddo

! store first estimate in first column
      if (iest .eq. 1) then
       do j=1,nv
        qcol(j,1) = yest(j)
       enddo
      else
       do j=1,nv
        d(j) = yest(j)
       enddo
       do k1=1,iest-1
        delta = 1.0d0/(x(iest-k1) - xest)
        f1    = xest * delta
        f2    = x(iest-k1) * delta

! propagate tableu 1 diagonal more
        do j=1,nv
         q          = qcol(j,k1)
         qcol(j,k1) = dy(j)
         delta      = d(j) - q
         dy(j)      = f1*delta
         d(j)       = f2*delta
         yz(j)      = yz(j) + dy(j)
        enddo
       enddo
       do j=1,nv
        qcol(j,iest) = dy(j)
       enddo
      end if
      return
      end
!---------------------------------------------------------------------


      subroutine stifbs_gift(y,dydx,nv,x,htry,eps,yscal,hdid,hnext, &
                              derivs,jakob,ierr)
      include 'implno.dek'
      include 'dense_matrix.dek'

! for dense analytic jacobians, lu decomposition linear algebra
!
! semi-implicit extrapolation step for integrating stiff odes with monitoring
! of local truncation error to adjust stepsize. inputs are the dependent
! variable vector y(1:nv) and its derivative dydx(1:nv) at the start of the
! independent variable x. also input are the stepsize to be attempted htry,
! the required accuracy eps, and the vector yscal against which the error is
! scaled. on output, y and x are replaced by their new values, hdid is the
! stepsize actually accomplished, and hnext is the estimated next stepsize.
! dervs is a user supplied function that computes the right hand side of
! the equations.

! declare
      external         derivs,jakob
      logical          first,reduct
      integer          nv,nmax,kmaxx,imax,ierr
      parameter        (nmax  = iodemax, &
                        kmaxx = 7, &
                        imax  = kmaxx+1)
      integer          i,iq,k,kk,km,kmax,kopt,nvold,nseq(imax)
      double precision y(nv),dydx(nv),x,htry,eps,yscal(nv),hdid,hnext, &
                       eps1,epsold,errmax,fact,h,red,scale,work,wrkmin, &
                       xest,xnew,a(imax),alf(kmaxx,kmaxx),err(kmaxx), &
                       yerr(nmax),ysav(nmax),yseq(nmax),safe1,safe2, &
                       redmax,redmin,tiny,scalmx,dum
      parameter        (safe1 = 0.25d0, safe2 = 0.7d0, redmax=1.0d-5, &
                        redmin = 0.7d0, tiny = 1.0d-30, scalmx = 0.1d0)


! assume that the independent variable is not explicit in the odes
      data             first/.true./, epsold/-1.0d0/, nvold/-1/
      data             nseq /2, 6, 10, 14, 22, 34, 50, 70/



! a new tolerance or a new number, so reinitialize
      if (eps .ne. epsold  .or.  nv .ne. nvold) then
       hnext = -1.0e29
       xnew  = -1.0e29
       eps1  = safe1 * eps

! compute the work coefficients a_k
       a(1)  = nseq(1) + 1
       do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
       enddo

! compute alf(k,q)
       do iq=2,kmaxx
        do k=1,iq-1
         alf(k,iq) = eps1**((a(k+1) - a(iq+1)) / &
                     ((a(iq+1) - a(1) + 1.0d0) * (2*k + 1)))
        enddo
       enddo
       epsold = eps
       nvold  = nv

! add cost of jacobians to work coefficients
       a(1)   = nv + a(1)
       do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
       enddo

! determine optimal row number for convergence
       do kopt=2,kmaxx-1
        if (a(kopt+1) .gt. a(kopt)*alf(kopt-1,kopt)) go to 01
       enddo
01     kmax = kopt
      end if

! save the starting values
      h    = htry
      do i=1,nv
       ysav(i)  = y(i)
      enddo

! get the dense jacobian in dens_dfdy
      call jakob(x,y,dens_dfdy,nv,iodemax)


! a new stepsize or a new integration, re-establish the order window
      if (h .ne. hnext  .or.  x .ne. xnew) then
       first = .true.
       kopt = kmax
      end if
      reduct = .false.

! evaluate the sequence of semi implicit midpoint rules
02    do 18 k=1,kmax
       xnew = x + h
       if (xnew .eq. x) stop 'stepsize too small in routine stiffbs'
       call simpr_gift13(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
       xest = (h/nseq(k))**2
       call net_pzextr(k,xest,yseq,y,yerr,nv)


! compute normalized error estimate
       if (k .ne. 1) then
        errmax = tiny
        ierr   = 0
        do i=1,nv
!        errmax = max(errmax,abs(yerr(i)/yscal(i)))
         dum = abs(yerr(i)/yscal(i))
         if (dum .ge. errmax) then
          errmax = dum
          ierr = i
         end if
        enddo

        errmax   = errmax/eps
        km = k - 1
        err(km) = (errmax/safe1)**(1.0d0/(2*km+1))
       end if

! in order window
       if (k .ne. 1  .and. (k .ge. kopt-1  .or. first)) then

! converged
        if (errmax .lt. 1.0) go to 04

! possible step size reductions
        if (k .eq. kmax  .or.  k .eq. kopt + 1) then
         red = safe2/err(km)
         go to 03
        else if (k .eq. kopt) then
         if (alf(kopt-1,kopt) .lt. err(km)) then
          red = 1.0d0/err(km)
          go to 03
         end if
        else if (kopt .eq. kmax) then
         if (alf(km,kmax-1) .lt. err(km)) then
          red = alf(km,kmax-1) * safe2/err(km)
          go to 03
         end if
        else if (alf(km,kopt) .lt. err(km)) then
         red = alf(km,kopt-1)/err(km)
         go to 03
        end if
       end if
18    continue

! reduce stepsize by at least redmin and at most redmax
03    red    = min(red,redmin)
      red    = max(red,redmax)
      h      = h * red
      ierr   = 0
      reduct = .true.
      go to 2

! successful step; get optimal row for convergence and corresponding stepsize
04    x = xnew
      hdid = h
      first = .false.
      wrkmin = 1.0e35
      do kk=1,km
       fact = max(err(kk),scalmx)
       work = fact * a(kk+1)
       if (work .lt. wrkmin) then
        scale  = fact
        wrkmin = work
        kopt   = kk + 1
       end if
      enddo

! check for possible order increase, but not if stepsize was just reduced
      hnext = h/scale
      if (kopt .ge. k  .and.  kopt .ne. kmax  .and.  .not.reduct) then
       fact = max(scale/alf(kopt-1,kopt),scalmx)
       if (a(kopt+1)*fact .le. wrkmin) then
        hnext = h/fact
        kopt = kopt + 1
       end if
      end if
      return
      end


      subroutine simpr_gift13(y,dydx,n,xs,htot,nstep,yout,derivs)
      include 'implno.dek'
      include 'dense_matrix.dek'
!
! an implicit midpoint stepper, for lu decomp dense linear algebra.
!

! declare
      external         derivs
      integer          n,nstep,nmaxx
      parameter        (nmaxx=iodemax)
      integer          i,j,nn
      double precision y(n),dydx(n),xs,htot, &
                       yout(n),h,x,del(nmaxx),ytemp(nmaxx)

! for the gift linear algebra
      integer          nmaxp1
      parameter        (nmaxp1 = iodemax + 1)
      double precision av(iodemax,nmaxp1)
      common /giftav/  av


! stepsize this trip, and make the a matrix
      h = htot/nstep
      do j=1,n
       do i=1,n
        dmat(i,j) = -h * dens_dfdy(i,j)
       enddo
      enddo
      do i=1,n
       dmat(i,i) = 1.0d0 + dmat(i,i)
      end do

! use yout as temporary storage; the first step
      do i=1,n
       yout(i) = h * dydx(i)
      enddo

      do j = 1, n
       do i = 1, n
        av(i,j) = dmat(i,j)
       enddo
      enddo
      do i = 1,n
       av(i,n+1) = yout(i)
      enddo
      call gift(av,iodemax,nmaxp1)
      do i = 1, n
       yout(i) = av(i,n+1)
      enddo

      do i=1,n
       del(i)   = yout(i)
       ytemp(i) = y(i) + del(i)
      enddo
      x = xs + h
      call derivs(x,ytemp,yout)

! use yout as temporary storage; general step
      do nn=2,nstep
       do i=1,n
        yout(i) = h*yout(i) - del(i)
       enddo

       do j = 1, n
        do i = 1, n
         av(i,j) = dmat(i,j)
        enddo
       enddo
       do i = 1,n
        av(i,n+1) = yout(i)
       enddo
       call gift(av,iodemax,nmaxp1)
       do i = 1, n
        yout(i) = av(i,n+1)
       enddo

       do i=1,n
        del(i)   = del(i) + 2.0d0 * yout(i)
        ytemp(i) = ytemp(i) + del(i)
       enddo
       x = x + h
       call derivs(x,ytemp,yout)
      enddo


! take the last step
      do i=1,n
       yout(i) = h * yout(i) - del(i)
      enddo

      do j = 1, n
       do i = 1, n
        av(i,j) = dmat(i,j)
       enddo
      enddo
      do i = 1,n
       av(i,n+1) = yout(i)
      enddo
      call gift(av,iodemax,nmaxp1)
      do i = 1, n
       yout(i) = av(i,n+1)
      enddo

      do i=1,n
       yout(i) = ytemp(i) + yout(i)
      enddo
      return
      end

!---------------------------------------------------------------------



      subroutine gift(ab,n1,n2)
      implicit none
      integer n1,n2
      double precision ab(n1,n2)
      call gsub_aprox13et_100(ab,n1,n2)
      return
      end




      subroutine gsub_aprox13et_100(ab,n1,n2)
      implicit none
      integer n1,n2
      double precision ab(n1,n2),c,tmp( 15)

        tmp( 15)   = 1.0d0/ab( 15, 15)
        c           = ab( 14, 15) * tmp(  15)
        ab( 14, 16) = ab( 14, 16) - ab( 15, 16) * c
        ab( 14,  1) = ab( 14,  1) - ab( 15,  1) * c
        ab( 14,  2) = ab( 14,  2) - ab( 15,  2) * c
        ab( 14,  3) = ab( 14,  3) - ab( 15,  3) * c
        ab( 14,  4) = ab( 14,  4) - ab( 15,  4) * c
        ab( 14,  5) = ab( 14,  5) - ab( 15,  5) * c
        ab( 14,  6) = ab( 14,  6) - ab( 15,  6) * c
        ab( 14,  7) = ab( 14,  7) - ab( 15,  7) * c
        ab( 14,  8) = ab( 14,  8) - ab( 15,  8) * c
        ab( 14,  9) = ab( 14,  9) - ab( 15,  9) * c
        ab( 14, 10) = ab( 14, 10) - ab( 15, 10) * c
        ab( 14, 11) = ab( 14, 11) - ab( 15, 11) * c
        ab( 14, 12) = ab( 14, 12) - ab( 15, 12) * c
        ab( 14, 13) = ab( 14, 13) - ab( 15, 13) * c
        c           = ab( 13, 15) * tmp(  15)
        ab( 13, 16) = ab( 13, 16) - ab( 15, 16) * c
        ab( 13,  1) = ab( 13,  1) - ab( 15,  1) * c
        ab( 13,  2) = ab( 13,  2) - ab( 15,  2) * c
        ab( 13,  3) = ab( 13,  3) - ab( 15,  3) * c
        ab( 13,  4) = ab( 13,  4) - ab( 15,  4) * c
        ab( 13,  5) = ab( 13,  5) - ab( 15,  5) * c
        ab( 13,  6) = ab( 13,  6) - ab( 15,  6) * c
        ab( 13,  7) = ab( 13,  7) - ab( 15,  7) * c
        ab( 13,  8) = ab( 13,  8) - ab( 15,  8) * c
        ab( 13,  9) = ab( 13,  9) - ab( 15,  9) * c
        ab( 13, 10) = ab( 13, 10) - ab( 15, 10) * c
        ab( 13, 11) = ab( 13, 11) - ab( 15, 11) * c
        ab( 13, 12) = ab( 13, 12) - ab( 15, 12) * c
        ab( 13, 13) = ab( 13, 13) - ab( 15, 13) * c
        c           = ab( 12, 15) * tmp(  15)
        ab( 12, 16) = ab( 12, 16) - ab( 15, 16) * c
        ab( 12,  1) = ab( 12,  1) - ab( 15,  1) * c
        ab( 12,  2) = ab( 12,  2) - ab( 15,  2) * c
        ab( 12,  3) = ab( 12,  3) - ab( 15,  3) * c
        ab( 12,  4) = ab( 12,  4) - ab( 15,  4) * c
        ab( 12,  5) = ab( 12,  5) - ab( 15,  5) * c
        ab( 12,  6) = ab( 12,  6) - ab( 15,  6) * c
        ab( 12,  7) = ab( 12,  7) - ab( 15,  7) * c
        ab( 12,  8) = ab( 12,  8) - ab( 15,  8) * c
        ab( 12,  9) = ab( 12,  9) - ab( 15,  9) * c
        ab( 12, 10) = ab( 12, 10) - ab( 15, 10) * c
        ab( 12, 11) = ab( 12, 11) - ab( 15, 11) * c
        ab( 12, 12) = ab( 12, 12) - ab( 15, 12) * c
        ab( 12, 13) = ab( 12, 13) - ab( 15, 13) * c
        c           = ab( 11, 15) * tmp(  15)
        ab( 11, 16) = ab( 11, 16) - ab( 15, 16) * c
        ab( 11,  1) = ab( 11,  1) - ab( 15,  1) * c
        ab( 11,  2) = ab( 11,  2) - ab( 15,  2) * c
        ab( 11,  3) = ab( 11,  3) - ab( 15,  3) * c
        ab( 11,  4) = ab( 11,  4) - ab( 15,  4) * c
        ab( 11,  5) = ab( 11,  5) - ab( 15,  5) * c
        ab( 11,  6) = ab( 11,  6) - ab( 15,  6) * c
        ab( 11,  7) = ab( 11,  7) - ab( 15,  7) * c
        ab( 11,  8) = ab( 11,  8) - ab( 15,  8) * c
        ab( 11,  9) = ab( 11,  9) - ab( 15,  9) * c
        ab( 11, 10) = ab( 11, 10) - ab( 15, 10) * c
        ab( 11, 11) = ab( 11, 11) - ab( 15, 11) * c
        ab( 11, 12) = ab( 11, 12) - ab( 15, 12) * c
        ab( 11, 13) = ab( 11, 13) - ab( 15, 13) * c
        c           = ab( 10, 15) * tmp(  15)
        ab( 10, 16) = ab( 10, 16) - ab( 15, 16) * c
        ab( 10,  1) = ab( 10,  1) - ab( 15,  1) * c
        ab( 10,  2) = ab( 10,  2) - ab( 15,  2) * c
        ab( 10,  3) = ab( 10,  3) - ab( 15,  3) * c
        ab( 10,  4) = ab( 10,  4) - ab( 15,  4) * c
        ab( 10,  5) = ab( 10,  5) - ab( 15,  5) * c
        ab( 10,  6) = ab( 10,  6) - ab( 15,  6) * c
        ab( 10,  7) = ab( 10,  7) - ab( 15,  7) * c
        ab( 10,  8) = ab( 10,  8) - ab( 15,  8) * c
        ab( 10,  9) = ab( 10,  9) - ab( 15,  9) * c
        ab( 10, 10) = ab( 10, 10) - ab( 15, 10) * c
        ab( 10, 11) = ab( 10, 11) - ab( 15, 11) * c
        ab( 10, 12) = ab( 10, 12) - ab( 15, 12) * c
        ab( 10, 13) = ab( 10, 13) - ab( 15, 13) * c
        c           = ab(  9, 15) * tmp(  15)
        ab(  9, 16) = ab(  9, 16) - ab( 15, 16) * c
        ab(  9,  1) = ab(  9,  1) - ab( 15,  1) * c
        ab(  9,  2) = ab(  9,  2) - ab( 15,  2) * c
        ab(  9,  3) = ab(  9,  3) - ab( 15,  3) * c
        ab(  9,  4) = ab(  9,  4) - ab( 15,  4) * c
        ab(  9,  5) = ab(  9,  5) - ab( 15,  5) * c
        ab(  9,  6) = ab(  9,  6) - ab( 15,  6) * c
        ab(  9,  7) = ab(  9,  7) - ab( 15,  7) * c
        ab(  9,  8) = ab(  9,  8) - ab( 15,  8) * c
        ab(  9,  9) = ab(  9,  9) - ab( 15,  9) * c
        ab(  9, 10) = ab(  9, 10) - ab( 15, 10) * c
        ab(  9, 11) = ab(  9, 11) - ab( 15, 11) * c
        ab(  9, 12) = ab(  9, 12) - ab( 15, 12) * c
        ab(  9, 13) = ab(  9, 13) - ab( 15, 13) * c
        c           = ab(  8, 15) * tmp(  15)
        ab(  8, 16) = ab(  8, 16) - ab( 15, 16) * c
        ab(  8,  1) = ab(  8,  1) - ab( 15,  1) * c
        ab(  8,  2) = ab(  8,  2) - ab( 15,  2) * c
        ab(  8,  3) = ab(  8,  3) - ab( 15,  3) * c
        ab(  8,  4) = ab(  8,  4) - ab( 15,  4) * c
        ab(  8,  5) = ab(  8,  5) - ab( 15,  5) * c
        ab(  8,  6) = ab(  8,  6) - ab( 15,  6) * c
        ab(  8,  7) = ab(  8,  7) - ab( 15,  7) * c
        ab(  8,  8) = ab(  8,  8) - ab( 15,  8) * c
        ab(  8,  9) = ab(  8,  9) - ab( 15,  9) * c
        ab(  8, 10) = ab(  8, 10) - ab( 15, 10) * c
        ab(  8, 11) = ab(  8, 11) - ab( 15, 11) * c
        ab(  8, 12) = ab(  8, 12) - ab( 15, 12) * c
        ab(  8, 13) = ab(  8, 13) - ab( 15, 13) * c
        c           = ab(  7, 15) * tmp(  15)
        ab(  7, 16) = ab(  7, 16) - ab( 15, 16) * c
        ab(  7,  1) = ab(  7,  1) - ab( 15,  1) * c
        ab(  7,  2) = ab(  7,  2) - ab( 15,  2) * c
        ab(  7,  3) = ab(  7,  3) - ab( 15,  3) * c
        ab(  7,  4) = ab(  7,  4) - ab( 15,  4) * c
        ab(  7,  5) = ab(  7,  5) - ab( 15,  5) * c
        ab(  7,  6) = ab(  7,  6) - ab( 15,  6) * c
        ab(  7,  7) = ab(  7,  7) - ab( 15,  7) * c
        ab(  7,  8) = ab(  7,  8) - ab( 15,  8) * c
        ab(  7,  9) = ab(  7,  9) - ab( 15,  9) * c
        ab(  7, 10) = ab(  7, 10) - ab( 15, 10) * c
        ab(  7, 11) = ab(  7, 11) - ab( 15, 11) * c
        ab(  7, 12) = ab(  7, 12) - ab( 15, 12) * c
        ab(  7, 13) = ab(  7, 13) - ab( 15, 13) * c
        c           = ab(  6, 15) * tmp(  15)
        ab(  6, 16) = ab(  6, 16) - ab( 15, 16) * c
        ab(  6,  1) = ab(  6,  1) - ab( 15,  1) * c
        ab(  6,  2) = ab(  6,  2) - ab( 15,  2) * c
        ab(  6,  3) = ab(  6,  3) - ab( 15,  3) * c
        ab(  6,  4) = ab(  6,  4) - ab( 15,  4) * c
        ab(  6,  5) = ab(  6,  5) - ab( 15,  5) * c
        ab(  6,  6) = ab(  6,  6) - ab( 15,  6) * c
        ab(  6,  7) = ab(  6,  7) - ab( 15,  7) * c
        ab(  6,  8) = ab(  6,  8) - ab( 15,  8) * c
        ab(  6,  9) = ab(  6,  9) - ab( 15,  9) * c
        ab(  6, 10) = ab(  6, 10) - ab( 15, 10) * c
        ab(  6, 11) = ab(  6, 11) - ab( 15, 11) * c
        ab(  6, 12) = ab(  6, 12) - ab( 15, 12) * c
        ab(  6, 13) = ab(  6, 13) - ab( 15, 13) * c
        c           = ab(  5, 15) * tmp(  15)
        ab(  5, 16) = ab(  5, 16) - ab( 15, 16) * c
        ab(  5,  1) = ab(  5,  1) - ab( 15,  1) * c
        ab(  5,  2) = ab(  5,  2) - ab( 15,  2) * c
        ab(  5,  3) = ab(  5,  3) - ab( 15,  3) * c
        ab(  5,  4) = ab(  5,  4) - ab( 15,  4) * c
        ab(  5,  5) = ab(  5,  5) - ab( 15,  5) * c
        ab(  5,  6) = ab(  5,  6) - ab( 15,  6) * c
        ab(  5,  7) = ab(  5,  7) - ab( 15,  7) * c
        ab(  5,  8) = ab(  5,  8) - ab( 15,  8) * c
        ab(  5,  9) = ab(  5,  9) - ab( 15,  9) * c
        ab(  5, 10) = ab(  5, 10) - ab( 15, 10) * c
        ab(  5, 11) = ab(  5, 11) - ab( 15, 11) * c
        ab(  5, 12) = ab(  5, 12) - ab( 15, 12) * c
        ab(  5, 13) = ab(  5, 13) - ab( 15, 13) * c
        c           = ab(  4, 15) * tmp(  15)
        ab(  4, 16) = ab(  4, 16) - ab( 15, 16) * c
        ab(  4,  1) = ab(  4,  1) - ab( 15,  1) * c
        ab(  4,  2) = ab(  4,  2) - ab( 15,  2) * c
        ab(  4,  3) = ab(  4,  3) - ab( 15,  3) * c
        ab(  4,  4) = ab(  4,  4) - ab( 15,  4) * c
        ab(  4,  5) = ab(  4,  5) - ab( 15,  5) * c
        ab(  4,  6) = ab(  4,  6) - ab( 15,  6) * c
        ab(  4,  7) = ab(  4,  7) - ab( 15,  7) * c
        ab(  4,  8) = ab(  4,  8) - ab( 15,  8) * c
        ab(  4,  9) = ab(  4,  9) - ab( 15,  9) * c
        ab(  4, 10) = ab(  4, 10) - ab( 15, 10) * c
        ab(  4, 11) = ab(  4, 11) - ab( 15, 11) * c
        ab(  4, 12) = ab(  4, 12) - ab( 15, 12) * c
        ab(  4, 13) = ab(  4, 13) - ab( 15, 13) * c
        c           = ab(  3, 15) * tmp(  15)
        ab(  3, 16) = ab(  3, 16) - ab( 15, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab( 15,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab( 15,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab( 15,  3) * c
        ab(  3,  4) = ab(  3,  4) - ab( 15,  4) * c
        ab(  3,  5) = ab(  3,  5) - ab( 15,  5) * c
        ab(  3,  6) = ab(  3,  6) - ab( 15,  6) * c
        ab(  3,  7) = ab(  3,  7) - ab( 15,  7) * c
        ab(  3,  8) = ab(  3,  8) - ab( 15,  8) * c
        ab(  3,  9) = ab(  3,  9) - ab( 15,  9) * c
        ab(  3, 10) = ab(  3, 10) - ab( 15, 10) * c
        ab(  3, 11) = ab(  3, 11) - ab( 15, 11) * c
        ab(  3, 12) = ab(  3, 12) - ab( 15, 12) * c
        ab(  3, 13) = ab(  3, 13) - ab( 15, 13) * c
        c           = ab(  2, 15) * tmp(  15)
        ab(  2, 16) = ab(  2, 16) - ab( 15, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab( 15,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab( 15,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab( 15,  3) * c
        ab(  2,  4) = ab(  2,  4) - ab( 15,  4) * c
        ab(  2,  5) = ab(  2,  5) - ab( 15,  5) * c
        ab(  2,  6) = ab(  2,  6) - ab( 15,  6) * c
        ab(  2,  7) = ab(  2,  7) - ab( 15,  7) * c
        ab(  2,  8) = ab(  2,  8) - ab( 15,  8) * c
        ab(  2,  9) = ab(  2,  9) - ab( 15,  9) * c
        ab(  2, 10) = ab(  2, 10) - ab( 15, 10) * c
        ab(  2, 11) = ab(  2, 11) - ab( 15, 11) * c
        ab(  2, 12) = ab(  2, 12) - ab( 15, 12) * c
        ab(  2, 13) = ab(  2, 13) - ab( 15, 13) * c
        c           = ab(  1, 15) * tmp(  15)
        ab(  1, 16) = ab(  1, 16) - ab( 15, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab( 15,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab( 15,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab( 15,  3) * c
        ab(  1,  4) = ab(  1,  4) - ab( 15,  4) * c
        ab(  1,  5) = ab(  1,  5) - ab( 15,  5) * c
        ab(  1,  6) = ab(  1,  6) - ab( 15,  6) * c
        ab(  1,  7) = ab(  1,  7) - ab( 15,  7) * c
        ab(  1,  8) = ab(  1,  8) - ab( 15,  8) * c
        ab(  1,  9) = ab(  1,  9) - ab( 15,  9) * c
        ab(  1, 10) = ab(  1, 10) - ab( 15, 10) * c
        ab(  1, 11) = ab(  1, 11) - ab( 15, 11) * c
        ab(  1, 12) = ab(  1, 12) - ab( 15, 12) * c
        ab(  1, 13) = ab(  1, 13) - ab( 15, 13) * c

        tmp( 14)   = 1.0d0/ab( 14, 14)

        tmp( 13)   = 1.0d0/ab( 13, 13)
        c           = ab( 12, 13) * tmp(  13)
        ab( 12, 16) = ab( 12, 16) - ab( 13, 16) * c
        ab( 12,  1) = ab( 12,  1) - ab( 13,  1) * c
        ab( 12,  2) = ab( 12,  2) - ab( 13,  2) * c
        ab( 12,  3) = ab( 12,  3) - ab( 13,  3) * c
        ab( 12,  4) = ab( 12,  4) - ab( 13,  4) * c
        ab( 12,  5) = ab( 12,  5) - ab( 13,  5) * c
        ab( 12,  6) = ab( 12,  6) - ab( 13,  6) * c
        ab( 12,  7) = ab( 12,  7) - ab( 13,  7) * c
        ab( 12,  8) = ab( 12,  8) - ab( 13,  8) * c
        ab( 12,  9) = ab( 12,  9) - ab( 13,  9) * c
        ab( 12, 10) = ab( 12, 10) - ab( 13, 10) * c
        ab( 12, 11) = ab( 12, 11) - ab( 13, 11) * c
        ab( 12, 12) = ab( 12, 12) - ab( 13, 12) * c
        c           = ab( 11, 13) * tmp(  13)
        ab( 11, 16) = ab( 11, 16) - ab( 13, 16) * c
        ab( 11,  1) = ab( 11,  1) - ab( 13,  1) * c
        ab( 11,  2) = ab( 11,  2) - ab( 13,  2) * c
        ab( 11,  3) = ab( 11,  3) - ab( 13,  3) * c
        ab( 11,  4) = ab( 11,  4) - ab( 13,  4) * c
        ab( 11,  5) = ab( 11,  5) - ab( 13,  5) * c
        ab( 11,  6) = ab( 11,  6) - ab( 13,  6) * c
        ab( 11,  7) = ab( 11,  7) - ab( 13,  7) * c
        ab( 11,  8) = ab( 11,  8) - ab( 13,  8) * c
        ab( 11,  9) = ab( 11,  9) - ab( 13,  9) * c
        ab( 11, 10) = ab( 11, 10) - ab( 13, 10) * c
        ab( 11, 11) = ab( 11, 11) - ab( 13, 11) * c
        ab( 11, 12) = ab( 11, 12) - ab( 13, 12) * c
        c           = ab( 10, 13) * tmp(  13)
        ab( 10, 16) = ab( 10, 16) - ab( 13, 16) * c
        ab( 10,  1) = ab( 10,  1) - ab( 13,  1) * c
        ab( 10,  2) = ab( 10,  2) - ab( 13,  2) * c
        ab( 10,  3) = ab( 10,  3) - ab( 13,  3) * c
        ab( 10,  4) = ab( 10,  4) - ab( 13,  4) * c
        ab( 10,  5) = ab( 10,  5) - ab( 13,  5) * c
        ab( 10,  6) = ab( 10,  6) - ab( 13,  6) * c
        ab( 10,  7) = ab( 10,  7) - ab( 13,  7) * c
        ab( 10,  8) = ab( 10,  8) - ab( 13,  8) * c
        ab( 10,  9) = ab( 10,  9) - ab( 13,  9) * c
        ab( 10, 10) = ab( 10, 10) - ab( 13, 10) * c
        ab( 10, 11) = ab( 10, 11) - ab( 13, 11) * c
        ab( 10, 12) = ab( 10, 12) - ab( 13, 12) * c
        c           = ab(  9, 13) * tmp(  13)
        ab(  9, 16) = ab(  9, 16) - ab( 13, 16) * c
        ab(  9,  1) = ab(  9,  1) - ab( 13,  1) * c
        ab(  9,  2) = ab(  9,  2) - ab( 13,  2) * c
        ab(  9,  3) = ab(  9,  3) - ab( 13,  3) * c
        ab(  9,  4) = ab(  9,  4) - ab( 13,  4) * c
        ab(  9,  5) = ab(  9,  5) - ab( 13,  5) * c
        ab(  9,  6) = ab(  9,  6) - ab( 13,  6) * c
        ab(  9,  7) = ab(  9,  7) - ab( 13,  7) * c
        ab(  9,  8) = ab(  9,  8) - ab( 13,  8) * c
        ab(  9,  9) = ab(  9,  9) - ab( 13,  9) * c
        ab(  9, 10) = ab(  9, 10) - ab( 13, 10) * c
        ab(  9, 11) = ab(  9, 11) - ab( 13, 11) * c
        ab(  9, 12) = ab(  9, 12) - ab( 13, 12) * c
        c           = ab(  8, 13) * tmp(  13)
        ab(  8, 16) = ab(  8, 16) - ab( 13, 16) * c
        ab(  8,  1) = ab(  8,  1) - ab( 13,  1) * c
        ab(  8,  2) = ab(  8,  2) - ab( 13,  2) * c
        ab(  8,  3) = ab(  8,  3) - ab( 13,  3) * c
        ab(  8,  4) = ab(  8,  4) - ab( 13,  4) * c
        ab(  8,  5) = ab(  8,  5) - ab( 13,  5) * c
        ab(  8,  6) = ab(  8,  6) - ab( 13,  6) * c
        ab(  8,  7) = ab(  8,  7) - ab( 13,  7) * c
        ab(  8,  8) = ab(  8,  8) - ab( 13,  8) * c
        ab(  8,  9) = ab(  8,  9) - ab( 13,  9) * c
        ab(  8, 10) = ab(  8, 10) - ab( 13, 10) * c
        ab(  8, 11) = ab(  8, 11) - ab( 13, 11) * c
        ab(  8, 12) = ab(  8, 12) - ab( 13, 12) * c
        c           = ab(  7, 13) * tmp(  13)
        ab(  7, 16) = ab(  7, 16) - ab( 13, 16) * c
        ab(  7,  1) = ab(  7,  1) - ab( 13,  1) * c
        ab(  7,  2) = ab(  7,  2) - ab( 13,  2) * c
        ab(  7,  3) = ab(  7,  3) - ab( 13,  3) * c
        ab(  7,  4) = ab(  7,  4) - ab( 13,  4) * c
        ab(  7,  5) = ab(  7,  5) - ab( 13,  5) * c
        ab(  7,  6) = ab(  7,  6) - ab( 13,  6) * c
        ab(  7,  7) = ab(  7,  7) - ab( 13,  7) * c
        ab(  7,  8) = ab(  7,  8) - ab( 13,  8) * c
        ab(  7,  9) = ab(  7,  9) - ab( 13,  9) * c
        ab(  7, 10) = ab(  7, 10) - ab( 13, 10) * c
        ab(  7, 11) = ab(  7, 11) - ab( 13, 11) * c
        ab(  7, 12) = ab(  7, 12) - ab( 13, 12) * c
        c           = ab(  6, 13) * tmp(  13)
        ab(  6, 16) = ab(  6, 16) - ab( 13, 16) * c
        ab(  6,  1) = ab(  6,  1) - ab( 13,  1) * c
        ab(  6,  2) = ab(  6,  2) - ab( 13,  2) * c
        ab(  6,  3) = ab(  6,  3) - ab( 13,  3) * c
        ab(  6,  4) = ab(  6,  4) - ab( 13,  4) * c
        ab(  6,  5) = ab(  6,  5) - ab( 13,  5) * c
        ab(  6,  6) = ab(  6,  6) - ab( 13,  6) * c
        ab(  6,  7) = ab(  6,  7) - ab( 13,  7) * c
        ab(  6,  8) = ab(  6,  8) - ab( 13,  8) * c
        ab(  6,  9) = ab(  6,  9) - ab( 13,  9) * c
        ab(  6, 10) = ab(  6, 10) - ab( 13, 10) * c
        ab(  6, 11) = ab(  6, 11) - ab( 13, 11) * c
        ab(  6, 12) = ab(  6, 12) - ab( 13, 12) * c
        c           = ab(  5, 13) * tmp(  13)
        ab(  5, 16) = ab(  5, 16) - ab( 13, 16) * c
        ab(  5,  1) = ab(  5,  1) - ab( 13,  1) * c
        ab(  5,  2) = ab(  5,  2) - ab( 13,  2) * c
        ab(  5,  3) = ab(  5,  3) - ab( 13,  3) * c
        ab(  5,  4) = ab(  5,  4) - ab( 13,  4) * c
        ab(  5,  5) = ab(  5,  5) - ab( 13,  5) * c
        ab(  5,  6) = ab(  5,  6) - ab( 13,  6) * c
        ab(  5,  7) = ab(  5,  7) - ab( 13,  7) * c
        ab(  5,  8) = ab(  5,  8) - ab( 13,  8) * c
        ab(  5,  9) = ab(  5,  9) - ab( 13,  9) * c
        ab(  5, 10) = ab(  5, 10) - ab( 13, 10) * c
        ab(  5, 11) = ab(  5, 11) - ab( 13, 11) * c
        ab(  5, 12) = ab(  5, 12) - ab( 13, 12) * c
        c           = ab(  4, 13) * tmp(  13)
        ab(  4, 16) = ab(  4, 16) - ab( 13, 16) * c
        ab(  4,  1) = ab(  4,  1) - ab( 13,  1) * c
        ab(  4,  2) = ab(  4,  2) - ab( 13,  2) * c
        ab(  4,  3) = ab(  4,  3) - ab( 13,  3) * c
        ab(  4,  4) = ab(  4,  4) - ab( 13,  4) * c
        ab(  4,  5) = ab(  4,  5) - ab( 13,  5) * c
        ab(  4,  6) = ab(  4,  6) - ab( 13,  6) * c
        ab(  4,  7) = ab(  4,  7) - ab( 13,  7) * c
        ab(  4,  8) = ab(  4,  8) - ab( 13,  8) * c
        ab(  4,  9) = ab(  4,  9) - ab( 13,  9) * c
        ab(  4, 10) = ab(  4, 10) - ab( 13, 10) * c
        ab(  4, 11) = ab(  4, 11) - ab( 13, 11) * c
        ab(  4, 12) = ab(  4, 12) - ab( 13, 12) * c
        c           = ab(  3, 13) * tmp(  13)
        ab(  3, 16) = ab(  3, 16) - ab( 13, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab( 13,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab( 13,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab( 13,  3) * c
        ab(  3,  4) = ab(  3,  4) - ab( 13,  4) * c
        ab(  3,  5) = ab(  3,  5) - ab( 13,  5) * c
        ab(  3,  6) = ab(  3,  6) - ab( 13,  6) * c
        ab(  3,  7) = ab(  3,  7) - ab( 13,  7) * c
        ab(  3,  8) = ab(  3,  8) - ab( 13,  8) * c
        ab(  3,  9) = ab(  3,  9) - ab( 13,  9) * c
        ab(  3, 10) = ab(  3, 10) - ab( 13, 10) * c
        ab(  3, 11) = ab(  3, 11) - ab( 13, 11) * c
        ab(  3, 12) = ab(  3, 12) - ab( 13, 12) * c
        c           = ab(  2, 13) * tmp(  13)
        ab(  2, 16) = ab(  2, 16) - ab( 13, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab( 13,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab( 13,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab( 13,  3) * c
        ab(  2,  4) = ab(  2,  4) - ab( 13,  4) * c
        ab(  2,  5) = ab(  2,  5) - ab( 13,  5) * c
        ab(  2,  6) = ab(  2,  6) - ab( 13,  6) * c
        ab(  2,  7) = ab(  2,  7) - ab( 13,  7) * c
        ab(  2,  8) = ab(  2,  8) - ab( 13,  8) * c
        ab(  2,  9) = ab(  2,  9) - ab( 13,  9) * c
        ab(  2, 10) = ab(  2, 10) - ab( 13, 10) * c
        ab(  2, 11) = ab(  2, 11) - ab( 13, 11) * c
        ab(  2, 12) = ab(  2, 12) - ab( 13, 12) * c
        c           = ab(  1, 13) * tmp(  13)
        ab(  1, 16) = ab(  1, 16) - ab( 13, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab( 13,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab( 13,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab( 13,  3) * c
        ab(  1,  4) = ab(  1,  4) - ab( 13,  4) * c
        ab(  1,  5) = ab(  1,  5) - ab( 13,  5) * c
        ab(  1,  6) = ab(  1,  6) - ab( 13,  6) * c
        ab(  1,  7) = ab(  1,  7) - ab( 13,  7) * c
        ab(  1,  8) = ab(  1,  8) - ab( 13,  8) * c
        ab(  1,  9) = ab(  1,  9) - ab( 13,  9) * c
        ab(  1, 10) = ab(  1, 10) - ab( 13, 10) * c
        ab(  1, 11) = ab(  1, 11) - ab( 13, 11) * c
        ab(  1, 12) = ab(  1, 12) - ab( 13, 12) * c

        tmp( 12)   = 1.0d0/ab( 12, 12)
        c           = ab( 11, 12) * tmp(  12)
        ab( 11, 16) = ab( 11, 16) - ab( 12, 16) * c
        ab( 11,  1) = ab( 11,  1) - ab( 12,  1) * c
        ab( 11,  2) = ab( 11,  2) - ab( 12,  2) * c
        ab( 11,  3) = ab( 11,  3) - ab( 12,  3) * c
        ab( 11,  4) = ab( 11,  4) - ab( 12,  4) * c
        ab( 11,  5) = ab( 11,  5) - ab( 12,  5) * c
        ab( 11,  6) = ab( 11,  6) - ab( 12,  6) * c
        ab( 11,  7) = ab( 11,  7) - ab( 12,  7) * c
        ab( 11,  8) = ab( 11,  8) - ab( 12,  8) * c
        ab( 11,  9) = ab( 11,  9) - ab( 12,  9) * c
        ab( 11, 10) = ab( 11, 10) - ab( 12, 10) * c
        ab( 11, 11) = ab( 11, 11) - ab( 12, 11) * c
        c           = ab( 10, 12) * tmp(  12)
        ab( 10, 16) = ab( 10, 16) - ab( 12, 16) * c
        ab( 10,  1) = ab( 10,  1) - ab( 12,  1) * c
        ab( 10,  2) = ab( 10,  2) - ab( 12,  2) * c
        ab( 10,  3) = ab( 10,  3) - ab( 12,  3) * c
        ab( 10,  4) = ab( 10,  4) - ab( 12,  4) * c
        ab( 10,  5) = ab( 10,  5) - ab( 12,  5) * c
        ab( 10,  6) = ab( 10,  6) - ab( 12,  6) * c
        ab( 10,  7) = ab( 10,  7) - ab( 12,  7) * c
        ab( 10,  8) = ab( 10,  8) - ab( 12,  8) * c
        ab( 10,  9) = ab( 10,  9) - ab( 12,  9) * c
        ab( 10, 10) = ab( 10, 10) - ab( 12, 10) * c
        ab( 10, 11) = ab( 10, 11) - ab( 12, 11) * c
        c           = ab(  9, 12) * tmp(  12)
        ab(  9, 16) = ab(  9, 16) - ab( 12, 16) * c
        ab(  9,  1) = ab(  9,  1) - ab( 12,  1) * c
        ab(  9,  2) = ab(  9,  2) - ab( 12,  2) * c
        ab(  9,  3) = ab(  9,  3) - ab( 12,  3) * c
        ab(  9,  4) = ab(  9,  4) - ab( 12,  4) * c
        ab(  9,  5) = ab(  9,  5) - ab( 12,  5) * c
        ab(  9,  6) = ab(  9,  6) - ab( 12,  6) * c
        ab(  9,  7) = ab(  9,  7) - ab( 12,  7) * c
        ab(  9,  8) = ab(  9,  8) - ab( 12,  8) * c
        ab(  9,  9) = ab(  9,  9) - ab( 12,  9) * c
        ab(  9, 10) = ab(  9, 10) - ab( 12, 10) * c
        ab(  9, 11) = ab(  9, 11) - ab( 12, 11) * c
        c           = ab(  8, 12) * tmp(  12)
        ab(  8, 16) = ab(  8, 16) - ab( 12, 16) * c
        ab(  8,  1) = ab(  8,  1) - ab( 12,  1) * c
        ab(  8,  2) = ab(  8,  2) - ab( 12,  2) * c
        ab(  8,  3) = ab(  8,  3) - ab( 12,  3) * c
        ab(  8,  4) = ab(  8,  4) - ab( 12,  4) * c
        ab(  8,  5) = ab(  8,  5) - ab( 12,  5) * c
        ab(  8,  6) = ab(  8,  6) - ab( 12,  6) * c
        ab(  8,  7) = ab(  8,  7) - ab( 12,  7) * c
        ab(  8,  8) = ab(  8,  8) - ab( 12,  8) * c
        ab(  8,  9) = ab(  8,  9) - ab( 12,  9) * c
        ab(  8, 10) = ab(  8, 10) - ab( 12, 10) * c
        ab(  8, 11) = ab(  8, 11) - ab( 12, 11) * c
        c           = ab(  7, 12) * tmp(  12)
        ab(  7, 16) = ab(  7, 16) - ab( 12, 16) * c
        ab(  7,  1) = ab(  7,  1) - ab( 12,  1) * c
        ab(  7,  2) = ab(  7,  2) - ab( 12,  2) * c
        ab(  7,  3) = ab(  7,  3) - ab( 12,  3) * c
        ab(  7,  4) = ab(  7,  4) - ab( 12,  4) * c
        ab(  7,  5) = ab(  7,  5) - ab( 12,  5) * c
        ab(  7,  6) = ab(  7,  6) - ab( 12,  6) * c
        ab(  7,  7) = ab(  7,  7) - ab( 12,  7) * c
        ab(  7,  8) = ab(  7,  8) - ab( 12,  8) * c
        ab(  7,  9) = ab(  7,  9) - ab( 12,  9) * c
        ab(  7, 10) = ab(  7, 10) - ab( 12, 10) * c
        ab(  7, 11) = ab(  7, 11) - ab( 12, 11) * c
        c           = ab(  6, 12) * tmp(  12)
        ab(  6, 16) = ab(  6, 16) - ab( 12, 16) * c
        ab(  6,  1) = ab(  6,  1) - ab( 12,  1) * c
        ab(  6,  2) = ab(  6,  2) - ab( 12,  2) * c
        ab(  6,  3) = ab(  6,  3) - ab( 12,  3) * c
        ab(  6,  4) = ab(  6,  4) - ab( 12,  4) * c
        ab(  6,  5) = ab(  6,  5) - ab( 12,  5) * c
        ab(  6,  6) = ab(  6,  6) - ab( 12,  6) * c
        ab(  6,  7) = ab(  6,  7) - ab( 12,  7) * c
        ab(  6,  8) = ab(  6,  8) - ab( 12,  8) * c
        ab(  6,  9) = ab(  6,  9) - ab( 12,  9) * c
        ab(  6, 10) = ab(  6, 10) - ab( 12, 10) * c
        ab(  6, 11) = ab(  6, 11) - ab( 12, 11) * c
        c           = ab(  5, 12) * tmp(  12)
        ab(  5, 16) = ab(  5, 16) - ab( 12, 16) * c
        ab(  5,  1) = ab(  5,  1) - ab( 12,  1) * c
        ab(  5,  2) = ab(  5,  2) - ab( 12,  2) * c
        ab(  5,  3) = ab(  5,  3) - ab( 12,  3) * c
        ab(  5,  4) = ab(  5,  4) - ab( 12,  4) * c
        ab(  5,  5) = ab(  5,  5) - ab( 12,  5) * c
        ab(  5,  6) = ab(  5,  6) - ab( 12,  6) * c
        ab(  5,  7) = ab(  5,  7) - ab( 12,  7) * c
        ab(  5,  8) = ab(  5,  8) - ab( 12,  8) * c
        ab(  5,  9) = ab(  5,  9) - ab( 12,  9) * c
        ab(  5, 10) = ab(  5, 10) - ab( 12, 10) * c
        ab(  5, 11) = ab(  5, 11) - ab( 12, 11) * c
        c           = ab(  4, 12) * tmp(  12)
        ab(  4, 16) = ab(  4, 16) - ab( 12, 16) * c
        ab(  4,  1) = ab(  4,  1) - ab( 12,  1) * c
        ab(  4,  2) = ab(  4,  2) - ab( 12,  2) * c
        ab(  4,  3) = ab(  4,  3) - ab( 12,  3) * c
        ab(  4,  4) = ab(  4,  4) - ab( 12,  4) * c
        ab(  4,  5) = ab(  4,  5) - ab( 12,  5) * c
        ab(  4,  6) = ab(  4,  6) - ab( 12,  6) * c
        ab(  4,  7) = ab(  4,  7) - ab( 12,  7) * c
        ab(  4,  8) = ab(  4,  8) - ab( 12,  8) * c
        ab(  4,  9) = ab(  4,  9) - ab( 12,  9) * c
        ab(  4, 10) = ab(  4, 10) - ab( 12, 10) * c
        ab(  4, 11) = ab(  4, 11) - ab( 12, 11) * c
        c           = ab(  3, 12) * tmp(  12)
        ab(  3, 16) = ab(  3, 16) - ab( 12, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab( 12,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab( 12,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab( 12,  3) * c
        ab(  3,  4) = ab(  3,  4) - ab( 12,  4) * c
        ab(  3,  5) = ab(  3,  5) - ab( 12,  5) * c
        ab(  3,  6) = ab(  3,  6) - ab( 12,  6) * c
        ab(  3,  7) = ab(  3,  7) - ab( 12,  7) * c
        ab(  3,  8) = ab(  3,  8) - ab( 12,  8) * c
        ab(  3,  9) = ab(  3,  9) - ab( 12,  9) * c
        ab(  3, 10) = ab(  3, 10) - ab( 12, 10) * c
        ab(  3, 11) = ab(  3, 11) - ab( 12, 11) * c
        c           = ab(  2, 12) * tmp(  12)
        ab(  2, 16) = ab(  2, 16) - ab( 12, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab( 12,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab( 12,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab( 12,  3) * c
        ab(  2,  4) = ab(  2,  4) - ab( 12,  4) * c
        ab(  2,  5) = ab(  2,  5) - ab( 12,  5) * c
        ab(  2,  6) = ab(  2,  6) - ab( 12,  6) * c
        ab(  2,  7) = ab(  2,  7) - ab( 12,  7) * c
        ab(  2,  8) = ab(  2,  8) - ab( 12,  8) * c
        ab(  2,  9) = ab(  2,  9) - ab( 12,  9) * c
        ab(  2, 10) = ab(  2, 10) - ab( 12, 10) * c
        ab(  2, 11) = ab(  2, 11) - ab( 12, 11) * c
        c           = ab(  1, 12) * tmp(  12)
        ab(  1, 16) = ab(  1, 16) - ab( 12, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab( 12,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab( 12,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab( 12,  3) * c
        ab(  1,  4) = ab(  1,  4) - ab( 12,  4) * c
        ab(  1,  5) = ab(  1,  5) - ab( 12,  5) * c
        ab(  1,  6) = ab(  1,  6) - ab( 12,  6) * c
        ab(  1,  7) = ab(  1,  7) - ab( 12,  7) * c
        ab(  1,  8) = ab(  1,  8) - ab( 12,  8) * c
        ab(  1,  9) = ab(  1,  9) - ab( 12,  9) * c
        ab(  1, 10) = ab(  1, 10) - ab( 12, 10) * c
        ab(  1, 11) = ab(  1, 11) - ab( 12, 11) * c

        tmp( 11)   = 1.0d0/ab( 11, 11)
        c           = ab( 10, 11) * tmp(  11)
        ab( 10, 16) = ab( 10, 16) - ab( 11, 16) * c
        ab( 10,  1) = ab( 10,  1) - ab( 11,  1) * c
        ab( 10,  2) = ab( 10,  2) - ab( 11,  2) * c
        ab( 10,  3) = ab( 10,  3) - ab( 11,  3) * c
        ab( 10,  4) = ab( 10,  4) - ab( 11,  4) * c
        ab( 10,  5) = ab( 10,  5) - ab( 11,  5) * c
        ab( 10,  6) = ab( 10,  6) - ab( 11,  6) * c
        ab( 10,  7) = ab( 10,  7) - ab( 11,  7) * c
        ab( 10,  8) = ab( 10,  8) - ab( 11,  8) * c
        ab( 10,  9) = ab( 10,  9) - ab( 11,  9) * c
        ab( 10, 10) = ab( 10, 10) - ab( 11, 10) * c
        c           = ab(  9, 11) * tmp(  11)
        ab(  9, 16) = ab(  9, 16) - ab( 11, 16) * c
        ab(  9,  1) = ab(  9,  1) - ab( 11,  1) * c
        ab(  9,  2) = ab(  9,  2) - ab( 11,  2) * c
        ab(  9,  3) = ab(  9,  3) - ab( 11,  3) * c
        ab(  9,  4) = ab(  9,  4) - ab( 11,  4) * c
        ab(  9,  5) = ab(  9,  5) - ab( 11,  5) * c
        ab(  9,  6) = ab(  9,  6) - ab( 11,  6) * c
        ab(  9,  7) = ab(  9,  7) - ab( 11,  7) * c
        ab(  9,  8) = ab(  9,  8) - ab( 11,  8) * c
        ab(  9,  9) = ab(  9,  9) - ab( 11,  9) * c
        ab(  9, 10) = ab(  9, 10) - ab( 11, 10) * c
        c           = ab(  8, 11) * tmp(  11)
        ab(  8, 16) = ab(  8, 16) - ab( 11, 16) * c
        ab(  8,  1) = ab(  8,  1) - ab( 11,  1) * c
        ab(  8,  2) = ab(  8,  2) - ab( 11,  2) * c
        ab(  8,  3) = ab(  8,  3) - ab( 11,  3) * c
        ab(  8,  4) = ab(  8,  4) - ab( 11,  4) * c
        ab(  8,  5) = ab(  8,  5) - ab( 11,  5) * c
        ab(  8,  6) = ab(  8,  6) - ab( 11,  6) * c
        ab(  8,  7) = ab(  8,  7) - ab( 11,  7) * c
        ab(  8,  8) = ab(  8,  8) - ab( 11,  8) * c
        ab(  8,  9) = ab(  8,  9) - ab( 11,  9) * c
        ab(  8, 10) = ab(  8, 10) - ab( 11, 10) * c
        c           = ab(  7, 11) * tmp(  11)
        ab(  7, 16) = ab(  7, 16) - ab( 11, 16) * c
        ab(  7,  1) = ab(  7,  1) - ab( 11,  1) * c
        ab(  7,  2) = ab(  7,  2) - ab( 11,  2) * c
        ab(  7,  3) = ab(  7,  3) - ab( 11,  3) * c
        ab(  7,  4) = ab(  7,  4) - ab( 11,  4) * c
        ab(  7,  5) = ab(  7,  5) - ab( 11,  5) * c
        ab(  7,  6) = ab(  7,  6) - ab( 11,  6) * c
        ab(  7,  7) = ab(  7,  7) - ab( 11,  7) * c
        ab(  7,  8) = ab(  7,  8) - ab( 11,  8) * c
        ab(  7,  9) = ab(  7,  9) - ab( 11,  9) * c
        ab(  7, 10) = ab(  7, 10) - ab( 11, 10) * c
        c           = ab(  6, 11) * tmp(  11)
        ab(  6, 16) = ab(  6, 16) - ab( 11, 16) * c
        ab(  6,  1) = ab(  6,  1) - ab( 11,  1) * c
        ab(  6,  2) = ab(  6,  2) - ab( 11,  2) * c
        ab(  6,  3) = ab(  6,  3) - ab( 11,  3) * c
        ab(  6,  4) = ab(  6,  4) - ab( 11,  4) * c
        ab(  6,  5) = ab(  6,  5) - ab( 11,  5) * c
        ab(  6,  6) = ab(  6,  6) - ab( 11,  6) * c
        ab(  6,  7) = ab(  6,  7) - ab( 11,  7) * c
        ab(  6,  8) = ab(  6,  8) - ab( 11,  8) * c
        ab(  6,  9) = ab(  6,  9) - ab( 11,  9) * c
        ab(  6, 10) = ab(  6, 10) - ab( 11, 10) * c
        c           = ab(  5, 11) * tmp(  11)
        ab(  5, 16) = ab(  5, 16) - ab( 11, 16) * c
        ab(  5,  1) = ab(  5,  1) - ab( 11,  1) * c
        ab(  5,  2) = ab(  5,  2) - ab( 11,  2) * c
        ab(  5,  3) = ab(  5,  3) - ab( 11,  3) * c
        ab(  5,  4) = ab(  5,  4) - ab( 11,  4) * c
        ab(  5,  5) = ab(  5,  5) - ab( 11,  5) * c
        ab(  5,  6) = ab(  5,  6) - ab( 11,  6) * c
        ab(  5,  7) = ab(  5,  7) - ab( 11,  7) * c
        ab(  5,  8) = ab(  5,  8) - ab( 11,  8) * c
        ab(  5,  9) = ab(  5,  9) - ab( 11,  9) * c
        ab(  5, 10) = ab(  5, 10) - ab( 11, 10) * c
        c           = ab(  4, 11) * tmp(  11)
        ab(  4, 16) = ab(  4, 16) - ab( 11, 16) * c
        ab(  4,  1) = ab(  4,  1) - ab( 11,  1) * c
        ab(  4,  2) = ab(  4,  2) - ab( 11,  2) * c
        ab(  4,  3) = ab(  4,  3) - ab( 11,  3) * c
        ab(  4,  4) = ab(  4,  4) - ab( 11,  4) * c
        ab(  4,  5) = ab(  4,  5) - ab( 11,  5) * c
        ab(  4,  6) = ab(  4,  6) - ab( 11,  6) * c
        ab(  4,  7) = ab(  4,  7) - ab( 11,  7) * c
        ab(  4,  8) = ab(  4,  8) - ab( 11,  8) * c
        ab(  4,  9) = ab(  4,  9) - ab( 11,  9) * c
        ab(  4, 10) = ab(  4, 10) - ab( 11, 10) * c
        c           = ab(  3, 11) * tmp(  11)
        ab(  3, 16) = ab(  3, 16) - ab( 11, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab( 11,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab( 11,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab( 11,  3) * c
        ab(  3,  4) = ab(  3,  4) - ab( 11,  4) * c
        ab(  3,  5) = ab(  3,  5) - ab( 11,  5) * c
        ab(  3,  6) = ab(  3,  6) - ab( 11,  6) * c
        ab(  3,  7) = ab(  3,  7) - ab( 11,  7) * c
        ab(  3,  8) = ab(  3,  8) - ab( 11,  8) * c
        ab(  3,  9) = ab(  3,  9) - ab( 11,  9) * c
        ab(  3, 10) = ab(  3, 10) - ab( 11, 10) * c
        c           = ab(  2, 11) * tmp(  11)
        ab(  2, 16) = ab(  2, 16) - ab( 11, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab( 11,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab( 11,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab( 11,  3) * c
        ab(  2,  4) = ab(  2,  4) - ab( 11,  4) * c
        ab(  2,  5) = ab(  2,  5) - ab( 11,  5) * c
        ab(  2,  6) = ab(  2,  6) - ab( 11,  6) * c
        ab(  2,  7) = ab(  2,  7) - ab( 11,  7) * c
        ab(  2,  8) = ab(  2,  8) - ab( 11,  8) * c
        ab(  2,  9) = ab(  2,  9) - ab( 11,  9) * c
        ab(  2, 10) = ab(  2, 10) - ab( 11, 10) * c
        c           = ab(  1, 11) * tmp(  11)
        ab(  1, 16) = ab(  1, 16) - ab( 11, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab( 11,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab( 11,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab( 11,  3) * c
        ab(  1,  4) = ab(  1,  4) - ab( 11,  4) * c
        ab(  1,  5) = ab(  1,  5) - ab( 11,  5) * c
        ab(  1,  6) = ab(  1,  6) - ab( 11,  6) * c
        ab(  1,  7) = ab(  1,  7) - ab( 11,  7) * c
        ab(  1,  8) = ab(  1,  8) - ab( 11,  8) * c
        ab(  1,  9) = ab(  1,  9) - ab( 11,  9) * c
        ab(  1, 10) = ab(  1, 10) - ab( 11, 10) * c

        tmp( 10)   = 1.0d0/ab( 10, 10)
        c           = ab(  9, 10) * tmp(  10)
        ab(  9, 16) = ab(  9, 16) - ab( 10, 16) * c
        ab(  9,  1) = ab(  9,  1) - ab( 10,  1) * c
        ab(  9,  2) = ab(  9,  2) - ab( 10,  2) * c
        ab(  9,  3) = ab(  9,  3) - ab( 10,  3) * c
        ab(  9,  4) = ab(  9,  4) - ab( 10,  4) * c
        ab(  9,  5) = ab(  9,  5) - ab( 10,  5) * c
        ab(  9,  6) = ab(  9,  6) - ab( 10,  6) * c
        ab(  9,  7) = ab(  9,  7) - ab( 10,  7) * c
        ab(  9,  8) = ab(  9,  8) - ab( 10,  8) * c
        ab(  9,  9) = ab(  9,  9) - ab( 10,  9) * c
        c           = ab(  8, 10) * tmp(  10)
        ab(  8, 16) = ab(  8, 16) - ab( 10, 16) * c
        ab(  8,  1) = ab(  8,  1) - ab( 10,  1) * c
        ab(  8,  2) = ab(  8,  2) - ab( 10,  2) * c
        ab(  8,  3) = ab(  8,  3) - ab( 10,  3) * c
        ab(  8,  4) = ab(  8,  4) - ab( 10,  4) * c
        ab(  8,  5) = ab(  8,  5) - ab( 10,  5) * c
        ab(  8,  6) = ab(  8,  6) - ab( 10,  6) * c
        ab(  8,  7) = ab(  8,  7) - ab( 10,  7) * c
        ab(  8,  8) = ab(  8,  8) - ab( 10,  8) * c
        ab(  8,  9) = ab(  8,  9) - ab( 10,  9) * c
        c           = ab(  7, 10) * tmp(  10)
        ab(  7, 16) = ab(  7, 16) - ab( 10, 16) * c
        ab(  7,  1) = ab(  7,  1) - ab( 10,  1) * c
        ab(  7,  2) = ab(  7,  2) - ab( 10,  2) * c
        ab(  7,  3) = ab(  7,  3) - ab( 10,  3) * c
        ab(  7,  4) = ab(  7,  4) - ab( 10,  4) * c
        ab(  7,  5) = ab(  7,  5) - ab( 10,  5) * c
        ab(  7,  6) = ab(  7,  6) - ab( 10,  6) * c
        ab(  7,  7) = ab(  7,  7) - ab( 10,  7) * c
        ab(  7,  8) = ab(  7,  8) - ab( 10,  8) * c
        ab(  7,  9) = ab(  7,  9) - ab( 10,  9) * c
        c           = ab(  6, 10) * tmp(  10)
        ab(  6, 16) = ab(  6, 16) - ab( 10, 16) * c
        ab(  6,  1) = ab(  6,  1) - ab( 10,  1) * c
        ab(  6,  2) = ab(  6,  2) - ab( 10,  2) * c
        ab(  6,  3) = ab(  6,  3) - ab( 10,  3) * c
        ab(  6,  4) = ab(  6,  4) - ab( 10,  4) * c
        ab(  6,  5) = ab(  6,  5) - ab( 10,  5) * c
        ab(  6,  6) = ab(  6,  6) - ab( 10,  6) * c
        ab(  6,  7) = ab(  6,  7) - ab( 10,  7) * c
        ab(  6,  8) = ab(  6,  8) - ab( 10,  8) * c
        ab(  6,  9) = ab(  6,  9) - ab( 10,  9) * c
        c           = ab(  5, 10) * tmp(  10)
        ab(  5, 16) = ab(  5, 16) - ab( 10, 16) * c
        ab(  5,  1) = ab(  5,  1) - ab( 10,  1) * c
        ab(  5,  2) = ab(  5,  2) - ab( 10,  2) * c
        ab(  5,  3) = ab(  5,  3) - ab( 10,  3) * c
        ab(  5,  4) = ab(  5,  4) - ab( 10,  4) * c
        ab(  5,  5) = ab(  5,  5) - ab( 10,  5) * c
        ab(  5,  6) = ab(  5,  6) - ab( 10,  6) * c
        ab(  5,  7) = ab(  5,  7) - ab( 10,  7) * c
        ab(  5,  8) = ab(  5,  8) - ab( 10,  8) * c
        ab(  5,  9) = ab(  5,  9) - ab( 10,  9) * c
        c           = ab(  4, 10) * tmp(  10)
        ab(  4, 16) = ab(  4, 16) - ab( 10, 16) * c
        ab(  4,  1) = ab(  4,  1) - ab( 10,  1) * c
        ab(  4,  2) = ab(  4,  2) - ab( 10,  2) * c
        ab(  4,  3) = ab(  4,  3) - ab( 10,  3) * c
        ab(  4,  4) = ab(  4,  4) - ab( 10,  4) * c
        ab(  4,  5) = ab(  4,  5) - ab( 10,  5) * c
        ab(  4,  6) = ab(  4,  6) - ab( 10,  6) * c
        ab(  4,  7) = ab(  4,  7) - ab( 10,  7) * c
        ab(  4,  8) = ab(  4,  8) - ab( 10,  8) * c
        ab(  4,  9) = ab(  4,  9) - ab( 10,  9) * c
        c           = ab(  3, 10) * tmp(  10)
        ab(  3, 16) = ab(  3, 16) - ab( 10, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab( 10,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab( 10,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab( 10,  3) * c
        ab(  3,  4) = ab(  3,  4) - ab( 10,  4) * c
        ab(  3,  5) = ab(  3,  5) - ab( 10,  5) * c
        ab(  3,  6) = ab(  3,  6) - ab( 10,  6) * c
        ab(  3,  7) = ab(  3,  7) - ab( 10,  7) * c
        ab(  3,  8) = ab(  3,  8) - ab( 10,  8) * c
        ab(  3,  9) = ab(  3,  9) - ab( 10,  9) * c
        c           = ab(  2, 10) * tmp(  10)
        ab(  2, 16) = ab(  2, 16) - ab( 10, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab( 10,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab( 10,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab( 10,  3) * c
        ab(  2,  4) = ab(  2,  4) - ab( 10,  4) * c
        ab(  2,  5) = ab(  2,  5) - ab( 10,  5) * c
        ab(  2,  6) = ab(  2,  6) - ab( 10,  6) * c
        ab(  2,  7) = ab(  2,  7) - ab( 10,  7) * c
        ab(  2,  8) = ab(  2,  8) - ab( 10,  8) * c
        ab(  2,  9) = ab(  2,  9) - ab( 10,  9) * c
        c           = ab(  1, 10) * tmp(  10)
        ab(  1, 16) = ab(  1, 16) - ab( 10, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab( 10,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab( 10,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab( 10,  3) * c
        ab(  1,  4) = ab(  1,  4) - ab( 10,  4) * c
        ab(  1,  5) = ab(  1,  5) - ab( 10,  5) * c
        ab(  1,  6) = ab(  1,  6) - ab( 10,  6) * c
        ab(  1,  7) = ab(  1,  7) - ab( 10,  7) * c
        ab(  1,  8) = ab(  1,  8) - ab( 10,  8) * c
        ab(  1,  9) = ab(  1,  9) - ab( 10,  9) * c

        tmp(  9)   = 1.0d0/ab(  9,  9)
        c           = ab(  8,  9) * tmp(   9)
        ab(  8, 16) = ab(  8, 16) - ab(  9, 16) * c
        ab(  8,  1) = ab(  8,  1) - ab(  9,  1) * c
        ab(  8,  2) = ab(  8,  2) - ab(  9,  2) * c
        ab(  8,  3) = ab(  8,  3) - ab(  9,  3) * c
        ab(  8,  4) = ab(  8,  4) - ab(  9,  4) * c
        ab(  8,  5) = ab(  8,  5) - ab(  9,  5) * c
        ab(  8,  6) = ab(  8,  6) - ab(  9,  6) * c
        ab(  8,  7) = ab(  8,  7) - ab(  9,  7) * c
        ab(  8,  8) = ab(  8,  8) - ab(  9,  8) * c
        c           = ab(  7,  9) * tmp(   9)
        ab(  7, 16) = ab(  7, 16) - ab(  9, 16) * c
        ab(  7,  1) = ab(  7,  1) - ab(  9,  1) * c
        ab(  7,  2) = ab(  7,  2) - ab(  9,  2) * c
        ab(  7,  3) = ab(  7,  3) - ab(  9,  3) * c
        ab(  7,  4) = ab(  7,  4) - ab(  9,  4) * c
        ab(  7,  5) = ab(  7,  5) - ab(  9,  5) * c
        ab(  7,  6) = ab(  7,  6) - ab(  9,  6) * c
        ab(  7,  7) = ab(  7,  7) - ab(  9,  7) * c
        ab(  7,  8) = ab(  7,  8) - ab(  9,  8) * c
        c           = ab(  6,  9) * tmp(   9)
        ab(  6, 16) = ab(  6, 16) - ab(  9, 16) * c
        ab(  6,  1) = ab(  6,  1) - ab(  9,  1) * c
        ab(  6,  2) = ab(  6,  2) - ab(  9,  2) * c
        ab(  6,  3) = ab(  6,  3) - ab(  9,  3) * c
        ab(  6,  4) = ab(  6,  4) - ab(  9,  4) * c
        ab(  6,  5) = ab(  6,  5) - ab(  9,  5) * c
        ab(  6,  6) = ab(  6,  6) - ab(  9,  6) * c
        ab(  6,  7) = ab(  6,  7) - ab(  9,  7) * c
        ab(  6,  8) = ab(  6,  8) - ab(  9,  8) * c
        c           = ab(  5,  9) * tmp(   9)
        ab(  5, 16) = ab(  5, 16) - ab(  9, 16) * c
        ab(  5,  1) = ab(  5,  1) - ab(  9,  1) * c
        ab(  5,  2) = ab(  5,  2) - ab(  9,  2) * c
        ab(  5,  3) = ab(  5,  3) - ab(  9,  3) * c
        ab(  5,  4) = ab(  5,  4) - ab(  9,  4) * c
        ab(  5,  5) = ab(  5,  5) - ab(  9,  5) * c
        ab(  5,  6) = ab(  5,  6) - ab(  9,  6) * c
        ab(  5,  7) = ab(  5,  7) - ab(  9,  7) * c
        ab(  5,  8) = ab(  5,  8) - ab(  9,  8) * c
        c           = ab(  4,  9) * tmp(   9)
        ab(  4, 16) = ab(  4, 16) - ab(  9, 16) * c
        ab(  4,  1) = ab(  4,  1) - ab(  9,  1) * c
        ab(  4,  2) = ab(  4,  2) - ab(  9,  2) * c
        ab(  4,  3) = ab(  4,  3) - ab(  9,  3) * c
        ab(  4,  4) = ab(  4,  4) - ab(  9,  4) * c
        ab(  4,  5) = ab(  4,  5) - ab(  9,  5) * c
        ab(  4,  6) = ab(  4,  6) - ab(  9,  6) * c
        ab(  4,  7) = ab(  4,  7) - ab(  9,  7) * c
        ab(  4,  8) = ab(  4,  8) - ab(  9,  8) * c
        c           = ab(  3,  9) * tmp(   9)
        ab(  3, 16) = ab(  3, 16) - ab(  9, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab(  9,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab(  9,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab(  9,  3) * c
        ab(  3,  4) = ab(  3,  4) - ab(  9,  4) * c
        ab(  3,  5) = ab(  3,  5) - ab(  9,  5) * c
        ab(  3,  6) = ab(  3,  6) - ab(  9,  6) * c
        ab(  3,  7) = ab(  3,  7) - ab(  9,  7) * c
        ab(  3,  8) = ab(  3,  8) - ab(  9,  8) * c
        c           = ab(  2,  9) * tmp(   9)
        ab(  2, 16) = ab(  2, 16) - ab(  9, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab(  9,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab(  9,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab(  9,  3) * c
        ab(  2,  4) = ab(  2,  4) - ab(  9,  4) * c
        ab(  2,  5) = ab(  2,  5) - ab(  9,  5) * c
        ab(  2,  6) = ab(  2,  6) - ab(  9,  6) * c
        ab(  2,  7) = ab(  2,  7) - ab(  9,  7) * c
        ab(  2,  8) = ab(  2,  8) - ab(  9,  8) * c
        c           = ab(  1,  9) * tmp(   9)
        ab(  1, 16) = ab(  1, 16) - ab(  9, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab(  9,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab(  9,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab(  9,  3) * c
        ab(  1,  4) = ab(  1,  4) - ab(  9,  4) * c
        ab(  1,  5) = ab(  1,  5) - ab(  9,  5) * c
        ab(  1,  6) = ab(  1,  6) - ab(  9,  6) * c
        ab(  1,  7) = ab(  1,  7) - ab(  9,  7) * c
        ab(  1,  8) = ab(  1,  8) - ab(  9,  8) * c

        tmp(  8)   = 1.0d0/ab(  8,  8)
        c           = ab(  7,  8) * tmp(   8)
        ab(  7, 16) = ab(  7, 16) - ab(  8, 16) * c
        ab(  7,  1) = ab(  7,  1) - ab(  8,  1) * c
        ab(  7,  2) = ab(  7,  2) - ab(  8,  2) * c
        ab(  7,  3) = ab(  7,  3) - ab(  8,  3) * c
        ab(  7,  4) = ab(  7,  4) - ab(  8,  4) * c
        ab(  7,  5) = ab(  7,  5) - ab(  8,  5) * c
        ab(  7,  6) = ab(  7,  6) - ab(  8,  6) * c
        ab(  7,  7) = ab(  7,  7) - ab(  8,  7) * c
        c           = ab(  6,  8) * tmp(   8)
        ab(  6, 16) = ab(  6, 16) - ab(  8, 16) * c
        ab(  6,  1) = ab(  6,  1) - ab(  8,  1) * c
        ab(  6,  2) = ab(  6,  2) - ab(  8,  2) * c
        ab(  6,  3) = ab(  6,  3) - ab(  8,  3) * c
        ab(  6,  4) = ab(  6,  4) - ab(  8,  4) * c
        ab(  6,  5) = ab(  6,  5) - ab(  8,  5) * c
        ab(  6,  6) = ab(  6,  6) - ab(  8,  6) * c
        ab(  6,  7) = ab(  6,  7) - ab(  8,  7) * c
        c           = ab(  5,  8) * tmp(   8)
        ab(  5, 16) = ab(  5, 16) - ab(  8, 16) * c
        ab(  5,  1) = ab(  5,  1) - ab(  8,  1) * c
        ab(  5,  2) = ab(  5,  2) - ab(  8,  2) * c
        ab(  5,  3) = ab(  5,  3) - ab(  8,  3) * c
        ab(  5,  4) = ab(  5,  4) - ab(  8,  4) * c
        ab(  5,  5) = ab(  5,  5) - ab(  8,  5) * c
        ab(  5,  6) = ab(  5,  6) - ab(  8,  6) * c
        ab(  5,  7) = ab(  5,  7) - ab(  8,  7) * c
        c           = ab(  4,  8) * tmp(   8)
        ab(  4, 16) = ab(  4, 16) - ab(  8, 16) * c
        ab(  4,  1) = ab(  4,  1) - ab(  8,  1) * c
        ab(  4,  2) = ab(  4,  2) - ab(  8,  2) * c
        ab(  4,  3) = ab(  4,  3) - ab(  8,  3) * c
        ab(  4,  4) = ab(  4,  4) - ab(  8,  4) * c
        ab(  4,  5) = ab(  4,  5) - ab(  8,  5) * c
        ab(  4,  6) = ab(  4,  6) - ab(  8,  6) * c
        ab(  4,  7) = ab(  4,  7) - ab(  8,  7) * c
        c           = ab(  3,  8) * tmp(   8)
        ab(  3, 16) = ab(  3, 16) - ab(  8, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab(  8,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab(  8,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab(  8,  3) * c
        ab(  3,  4) = ab(  3,  4) - ab(  8,  4) * c
        ab(  3,  5) = ab(  3,  5) - ab(  8,  5) * c
        ab(  3,  6) = ab(  3,  6) - ab(  8,  6) * c
        ab(  3,  7) = ab(  3,  7) - ab(  8,  7) * c
        c           = ab(  2,  8) * tmp(   8)
        ab(  2, 16) = ab(  2, 16) - ab(  8, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab(  8,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab(  8,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab(  8,  3) * c
        ab(  2,  4) = ab(  2,  4) - ab(  8,  4) * c
        ab(  2,  5) = ab(  2,  5) - ab(  8,  5) * c
        ab(  2,  6) = ab(  2,  6) - ab(  8,  6) * c
        ab(  2,  7) = ab(  2,  7) - ab(  8,  7) * c
        c           = ab(  1,  8) * tmp(   8)
        ab(  1, 16) = ab(  1, 16) - ab(  8, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab(  8,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab(  8,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab(  8,  3) * c
        ab(  1,  4) = ab(  1,  4) - ab(  8,  4) * c
        ab(  1,  5) = ab(  1,  5) - ab(  8,  5) * c
        ab(  1,  6) = ab(  1,  6) - ab(  8,  6) * c
        ab(  1,  7) = ab(  1,  7) - ab(  8,  7) * c

        tmp(  7)   = 1.0d0/ab(  7,  7)
        c           = ab(  6,  7) * tmp(   7)
        ab(  6, 16) = ab(  6, 16) - ab(  7, 16) * c
        ab(  6,  1) = ab(  6,  1) - ab(  7,  1) * c
        ab(  6,  2) = ab(  6,  2) - ab(  7,  2) * c
        ab(  6,  3) = ab(  6,  3) - ab(  7,  3) * c
        ab(  6,  4) = ab(  6,  4) - ab(  7,  4) * c
        ab(  6,  5) = ab(  6,  5) - ab(  7,  5) * c
        ab(  6,  6) = ab(  6,  6) - ab(  7,  6) * c
        c           = ab(  5,  7) * tmp(   7)
        ab(  5, 16) = ab(  5, 16) - ab(  7, 16) * c
        ab(  5,  1) = ab(  5,  1) - ab(  7,  1) * c
        ab(  5,  2) = ab(  5,  2) - ab(  7,  2) * c
        ab(  5,  3) = ab(  5,  3) - ab(  7,  3) * c
        ab(  5,  4) = ab(  5,  4) - ab(  7,  4) * c
        ab(  5,  5) = ab(  5,  5) - ab(  7,  5) * c
        ab(  5,  6) = ab(  5,  6) - ab(  7,  6) * c
        c           = ab(  4,  7) * tmp(   7)
        ab(  4, 16) = ab(  4, 16) - ab(  7, 16) * c
        ab(  4,  1) = ab(  4,  1) - ab(  7,  1) * c
        ab(  4,  2) = ab(  4,  2) - ab(  7,  2) * c
        ab(  4,  3) = ab(  4,  3) - ab(  7,  3) * c
        ab(  4,  4) = ab(  4,  4) - ab(  7,  4) * c
        ab(  4,  5) = ab(  4,  5) - ab(  7,  5) * c
        ab(  4,  6) = ab(  4,  6) - ab(  7,  6) * c
        c           = ab(  3,  7) * tmp(   7)
        ab(  3, 16) = ab(  3, 16) - ab(  7, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab(  7,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab(  7,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab(  7,  3) * c
        ab(  3,  4) = ab(  3,  4) - ab(  7,  4) * c
        ab(  3,  5) = ab(  3,  5) - ab(  7,  5) * c
        ab(  3,  6) = ab(  3,  6) - ab(  7,  6) * c
        c           = ab(  2,  7) * tmp(   7)
        ab(  2, 16) = ab(  2, 16) - ab(  7, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab(  7,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab(  7,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab(  7,  3) * c
        ab(  2,  4) = ab(  2,  4) - ab(  7,  4) * c
        ab(  2,  5) = ab(  2,  5) - ab(  7,  5) * c
        ab(  2,  6) = ab(  2,  6) - ab(  7,  6) * c
        c           = ab(  1,  7) * tmp(   7)
        ab(  1, 16) = ab(  1, 16) - ab(  7, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab(  7,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab(  7,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab(  7,  3) * c
        ab(  1,  4) = ab(  1,  4) - ab(  7,  4) * c
        ab(  1,  5) = ab(  1,  5) - ab(  7,  5) * c
        ab(  1,  6) = ab(  1,  6) - ab(  7,  6) * c

        tmp(  6)   = 1.0d0/ab(  6,  6)
        c           = ab(  5,  6) * tmp(   6)
        ab(  5, 16) = ab(  5, 16) - ab(  6, 16) * c
        ab(  5,  1) = ab(  5,  1) - ab(  6,  1) * c
        ab(  5,  2) = ab(  5,  2) - ab(  6,  2) * c
        ab(  5,  3) = ab(  5,  3) - ab(  6,  3) * c
        ab(  5,  4) = ab(  5,  4) - ab(  6,  4) * c
        ab(  5,  5) = ab(  5,  5) - ab(  6,  5) * c
        c           = ab(  4,  6) * tmp(   6)
        ab(  4, 16) = ab(  4, 16) - ab(  6, 16) * c
        ab(  4,  1) = ab(  4,  1) - ab(  6,  1) * c
        ab(  4,  2) = ab(  4,  2) - ab(  6,  2) * c
        ab(  4,  3) = ab(  4,  3) - ab(  6,  3) * c
        ab(  4,  4) = ab(  4,  4) - ab(  6,  4) * c
        ab(  4,  5) = ab(  4,  5) - ab(  6,  5) * c
        c           = ab(  3,  6) * tmp(   6)
        ab(  3, 16) = ab(  3, 16) - ab(  6, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab(  6,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab(  6,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab(  6,  3) * c
        ab(  3,  4) = ab(  3,  4) - ab(  6,  4) * c
        ab(  3,  5) = ab(  3,  5) - ab(  6,  5) * c
        c           = ab(  2,  6) * tmp(   6)
        ab(  2, 16) = ab(  2, 16) - ab(  6, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab(  6,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab(  6,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab(  6,  3) * c
        ab(  2,  4) = ab(  2,  4) - ab(  6,  4) * c
        ab(  2,  5) = ab(  2,  5) - ab(  6,  5) * c
        c           = ab(  1,  6) * tmp(   6)
        ab(  1, 16) = ab(  1, 16) - ab(  6, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab(  6,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab(  6,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab(  6,  3) * c
        ab(  1,  4) = ab(  1,  4) - ab(  6,  4) * c
        ab(  1,  5) = ab(  1,  5) - ab(  6,  5) * c

        tmp(  5)   = 1.0d0/ab(  5,  5)
        c           = ab(  4,  5) * tmp(   5)
        ab(  4, 16) = ab(  4, 16) - ab(  5, 16) * c
        ab(  4,  1) = ab(  4,  1) - ab(  5,  1) * c
        ab(  4,  2) = ab(  4,  2) - ab(  5,  2) * c
        ab(  4,  3) = ab(  4,  3) - ab(  5,  3) * c
        ab(  4,  4) = ab(  4,  4) - ab(  5,  4) * c
        c           = ab(  3,  5) * tmp(   5)
        ab(  3, 16) = ab(  3, 16) - ab(  5, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab(  5,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab(  5,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab(  5,  3) * c
        ab(  3,  4) = ab(  3,  4) - ab(  5,  4) * c
        c           = ab(  2,  5) * tmp(   5)
        ab(  2, 16) = ab(  2, 16) - ab(  5, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab(  5,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab(  5,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab(  5,  3) * c
        ab(  2,  4) = ab(  2,  4) - ab(  5,  4) * c
        c           = ab(  1,  5) * tmp(   5)
        ab(  1, 16) = ab(  1, 16) - ab(  5, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab(  5,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab(  5,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab(  5,  3) * c
        ab(  1,  4) = ab(  1,  4) - ab(  5,  4) * c

        tmp(  4)   = 1.0d0/ab(  4,  4)
        c           = ab(  3,  4) * tmp(   4)
        ab(  3, 16) = ab(  3, 16) - ab(  4, 16) * c
        ab(  3,  1) = ab(  3,  1) - ab(  4,  1) * c
        ab(  3,  2) = ab(  3,  2) - ab(  4,  2) * c
        ab(  3,  3) = ab(  3,  3) - ab(  4,  3) * c
        c           = ab(  2,  4) * tmp(   4)
        ab(  2, 16) = ab(  2, 16) - ab(  4, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab(  4,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab(  4,  2) * c
        ab(  2,  3) = ab(  2,  3) - ab(  4,  3) * c
        c           = ab(  1,  4) * tmp(   4)
        ab(  1, 16) = ab(  1, 16) - ab(  4, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab(  4,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab(  4,  2) * c
        ab(  1,  3) = ab(  1,  3) - ab(  4,  3) * c

        tmp(  3)   = 1.0d0/ab(  3,  3)
        c           = ab(  2,  3) * tmp(   3)
        ab(  2, 16) = ab(  2, 16) - ab(  3, 16) * c
        ab(  2,  1) = ab(  2,  1) - ab(  3,  1) * c
        ab(  2,  2) = ab(  2,  2) - ab(  3,  2) * c
        c           = ab(  1,  3) * tmp(   3)
        ab(  1, 16) = ab(  1, 16) - ab(  3, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab(  3,  1) * c
        ab(  1,  2) = ab(  1,  2) - ab(  3,  2) * c

        tmp(  2)   = 1.0d0/ab(  2,  2)
        c           = ab(  1,  2) * tmp(   2)
        ab(  1, 16) = ab(  1, 16) - ab(  2, 16) * c
        ab(  1,  1) = ab(  1,  1) - ab(  2,  1) * c


! backsubstitution
        tmp(  1)     = 1.0d0/ab(  1,  1)
        ab(  1, 16) =   ab(  1, 16) * tmp(   1)
        ab(  2, 16) = ( ab(  2, 16) &
                        - ab(  2,  1) * ab(  1, 16) &
                                                    ) *tmp(   2)
        ab(  3, 16) = ( ab(  3, 16) &
                        - ab(  3,  1) * ab(  1, 16) &
                        - ab(  3,  2) * ab(  2, 16) &
                                                    ) *tmp(   3)
        ab(  4, 16) = ( ab(  4, 16) &
                        - ab(  4,  1) * ab(  1, 16) &
                        - ab(  4,  2) * ab(  2, 16) &
                        - ab(  4,  3) * ab(  3, 16) &
                                                    ) *tmp(   4)
        ab(  5, 16) = ( ab(  5, 16) &
                        - ab(  5,  1) * ab(  1, 16) &
                        - ab(  5,  2) * ab(  2, 16) &
                        - ab(  5,  3) * ab(  3, 16) &
                        - ab(  5,  4) * ab(  4, 16) &
                                                    ) *tmp(   5)
        ab(  6, 16) = ( ab(  6, 16) &
                        - ab(  6,  1) * ab(  1, 16) &
                        - ab(  6,  2) * ab(  2, 16) &
                        - ab(  6,  3) * ab(  3, 16) &
                        - ab(  6,  4) * ab(  4, 16) &
                        - ab(  6,  5) * ab(  5, 16) &
                                                    ) *tmp(   6)
        ab(  7, 16) = ( ab(  7, 16) &
                        - ab(  7,  1) * ab(  1, 16) &
                        - ab(  7,  2) * ab(  2, 16) &
                        - ab(  7,  3) * ab(  3, 16) &
                        - ab(  7,  4) * ab(  4, 16) &
                        - ab(  7,  5) * ab(  5, 16) &
                        - ab(  7,  6) * ab(  6, 16) &
                                                    ) *tmp(   7)
        ab(  8, 16) = ( ab(  8, 16) &
                        - ab(  8,  1) * ab(  1, 16) &
                        - ab(  8,  2) * ab(  2, 16) &
                        - ab(  8,  3) * ab(  3, 16) &
                        - ab(  8,  4) * ab(  4, 16) &
                        - ab(  8,  5) * ab(  5, 16) &
                        - ab(  8,  6) * ab(  6, 16) &
                        - ab(  8,  7) * ab(  7, 16) &
                                                    ) *tmp(   8)
        ab(  9, 16) = ( ab(  9, 16) &
                        - ab(  9,  1) * ab(  1, 16) &
                        - ab(  9,  2) * ab(  2, 16) &
                        - ab(  9,  3) * ab(  3, 16) &
                        - ab(  9,  4) * ab(  4, 16) &
                        - ab(  9,  5) * ab(  5, 16) &
                        - ab(  9,  6) * ab(  6, 16) &
                        - ab(  9,  7) * ab(  7, 16) &
                        - ab(  9,  8) * ab(  8, 16) &
                                                    ) *tmp(   9)
        ab( 10, 16) = ( ab( 10, 16) &
                        - ab( 10,  1) * ab(  1, 16) &
                        - ab( 10,  2) * ab(  2, 16) &
                        - ab( 10,  3) * ab(  3, 16) &
                        - ab( 10,  4) * ab(  4, 16) &
                        - ab( 10,  5) * ab(  5, 16) &
                        - ab( 10,  6) * ab(  6, 16) &
                        - ab( 10,  7) * ab(  7, 16) &
                        - ab( 10,  8) * ab(  8, 16) &
                        - ab( 10,  9) * ab(  9, 16) &
                                                    ) *tmp(  10)
        ab( 11, 16) = ( ab( 11, 16) &
                        - ab( 11,  1) * ab(  1, 16) &
                        - ab( 11,  2) * ab(  2, 16) &
                        - ab( 11,  3) * ab(  3, 16) &
                        - ab( 11,  4) * ab(  4, 16) &
                        - ab( 11,  5) * ab(  5, 16) &
                        - ab( 11,  6) * ab(  6, 16) &
                        - ab( 11,  7) * ab(  7, 16) &
                        - ab( 11,  8) * ab(  8, 16) &
                        - ab( 11,  9) * ab(  9, 16) &
                        - ab( 11, 10) * ab( 10, 16) &
                                                    ) *tmp(  11)
        ab( 12, 16) = ( ab( 12, 16) &
                        - ab( 12,  1) * ab(  1, 16) &
                        - ab( 12,  2) * ab(  2, 16) &
                        - ab( 12,  3) * ab(  3, 16) &
                        - ab( 12,  4) * ab(  4, 16) &
                        - ab( 12,  5) * ab(  5, 16) &
                        - ab( 12,  6) * ab(  6, 16) &
                        - ab( 12,  7) * ab(  7, 16) &
                        - ab( 12,  8) * ab(  8, 16) &
                        - ab( 12,  9) * ab(  9, 16) &
                        - ab( 12, 10) * ab( 10, 16) &
                        - ab( 12, 11) * ab( 11, 16) &
                                                    ) *tmp(  12)
        ab( 13, 16) = ( ab( 13, 16) &
                        - ab( 13,  1) * ab(  1, 16) &
                        - ab( 13,  2) * ab(  2, 16) &
                        - ab( 13,  3) * ab(  3, 16) &
                        - ab( 13,  4) * ab(  4, 16) &
                        - ab( 13,  5) * ab(  5, 16) &
                        - ab( 13,  6) * ab(  6, 16) &
                        - ab( 13,  7) * ab(  7, 16) &
                        - ab( 13,  8) * ab(  8, 16) &
                        - ab( 13,  9) * ab(  9, 16) &
                        - ab( 13, 10) * ab( 10, 16) &
                        - ab( 13, 11) * ab( 11, 16) &
                        - ab( 13, 12) * ab( 12, 16) &
                                                    ) *tmp(  13)
        ab( 14, 16) = ( ab( 14, 16) &
                        - ab( 14,  1) * ab(  1, 16) &
                        - ab( 14,  2) * ab(  2, 16) &
                        - ab( 14,  3) * ab(  3, 16) &
                        - ab( 14,  4) * ab(  4, 16) &
                        - ab( 14,  5) * ab(  5, 16) &
                        - ab( 14,  6) * ab(  6, 16) &
                        - ab( 14,  7) * ab(  7, 16) &
                        - ab( 14,  8) * ab(  8, 16) &
                        - ab( 14,  9) * ab(  9, 16) &
                        - ab( 14, 10) * ab( 10, 16) &
                        - ab( 14, 11) * ab( 11, 16) &
                        - ab( 14, 12) * ab( 12, 16) &
                        - ab( 14, 13) * ab( 13, 16) &
                                                    ) *tmp(  14)
        ab( 15, 16) = ( ab( 15, 16) &
                        - ab( 15,  1) * ab(  1, 16) &
                        - ab( 15,  2) * ab(  2, 16) &
                        - ab( 15,  3) * ab(  3, 16) &
                        - ab( 15,  4) * ab(  4, 16) &
                        - ab( 15,  5) * ab(  5, 16) &
                        - ab( 15,  6) * ab(  6, 16) &
                        - ab( 15,  7) * ab(  7, 16) &
                        - ab( 15,  8) * ab(  8, 16) &
                        - ab( 15,  9) * ab(  9, 16) &
                        - ab( 15, 10) * ab( 10, 16) &
                        - ab( 15, 11) * ab( 11, 16) &
                        - ab( 15, 12) * ab( 12, 16) &
                        - ab( 15, 13) * ab( 13, 16) &
                                                    ) *tmp(  15)
      return
      end



!---------------------------------------------------------------------


      subroutine azbar_burn(xmass,aion,zion,wion,ionmax, &
                       ymass,abar,zbar,wbar,ye,nxcess)
      include 'implno.dek'

! this routine calculates composition variables

! input:
! mass fractions               = xmass(1:ionmax)  dimensionless
! number of nucleons           = aion(1:ionmax)   dimensionless
! charge of nucleus            = zion(1:ionmax)   dimensionless
! atomic weight or molar mass  = wion(1:ionmax)    g/mole
! number of isotopes           = ionmax
!
! output:
! molar abundances        = ymass(1:ionmax)   mole/g
! mean number of nucleons = abar              dimensionless
! mean nucleon charge     = zbar              dimensionless
! mean weight             = wbar              g/mole
! electron fraction       = ye                mole/g
! neutron excess          = xcess


! declare the pass
      integer          ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax), &
                       wion(ionmax),ymass(ionmax),abar,zbar,wbar, &
                       ye,nxcess

! local variables
      double precision asum,sum1


! molar abundances
      ymass(1:ionmax) = xmass(1:ionmax)/wion(1:ionmax)

! mean molar mass
      wbar  = 1.0d0/sum(ymass(1:ionmax))

! mean number of nucleons
      sum1  = sum(aion(1:ionmax)*ymass(1:ionmax))
      abar  = wbar * sum1

! mean charge
      ye  = sum(zion(1:ionmax)*ymass(1:ionmax))
      zbar  = wbar * ye

! neutron excess
      nxcess = sum1 - 2.0d0 * ye

      return
      end





      subroutine ener_gener_rate(dydt,enuc)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! computes the instantaneous energy generation rate

! declare the pass
       double precision dydt(1),enuc

! local variables
       integer          i

! conversion factors for the nuclear energy generation rate
! detlap is the mass excess of the proton in amu
! detlan is the mass excess of the neutron in amu

      double precision enuc_conv,enuc_conv2,deltap,deltan
      parameter        (enuc_conv  = ev2erg*1.0d6*avo, &
                        enuc_conv2 = -avo*clight*clight, &
                        deltap     = 7.288969d0, &
                        deltan     = 8.071323d0)


! instantaneous energy generation rate
!      enuc = enuc_conv * sum(dydt(ionbeg:ionend)*bion(ionbeg:ionend))
!      enuc = enuc_conv * sum(dydt(ionbeg:ionend) * (bion(ionbeg:ionend) - zion(ionbeg:ionend)*deltap - nion(ionbeg:ionend)*deltan))
      enuc = enuc_conv2 * sum(dydt(ionbeg:ionend)*mion(ionbeg:ionend))

      return
      end


