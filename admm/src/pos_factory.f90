module pos_factory

	!--------------------------------------------------------
	! Factory for proximal operators 
	! Uses submodules instead of procedural pointer.
	!
	! Laurence Yang, UCSD
	!
	! 30 Aug 2018: first port from po_factory.f90
	! 09 Oct 2018: error stop --> stop to support PGI compiler
	! 03 Feb 2019: fixed preconditioner bug in bilinear PO
	!--------------------------------------------------------
	!
	! 1:	bounds (l <= x <= u)
	! 2:	bilinear (x'y <=> d)
	! 3:	L1
	! 4:    Linear system
	!
	!--------------------------------------------------------

	!--------------------------------------------------------
	! The collection of proximal operators: 
	!--------------------------------------------------------
	use po_constants
	use po_bounds, 			only: opt_boundss,		&
								  lb, ub
	use po_bilinear, 		only: opt_bilinears,	&
								  d, csense, nx, xstart, xend, ystart, yend 

	use po_bilinear_multi, 	only: ds, csenses, xstarts, xends, ystarts, yends, nxs, &
								  opt_bilinear_multi

	use po_L1, 				only: opt_L1s, calc_obj_L1s
	use po_card,			only: opt_cards
	use po_linsys, 			only: opt_linsyss, calc_obj_linsys,	&
								  linsys_mF, cleanup_linsys
	use po_linsys_ma57,		only: opt_linsys57, calc_obj_linsys57
	use po_bounds_L1,		only: opt_bounds_L1s
	use po_bounds_card,		only: opt_bounds_cards
	use po_linsys_fg,		only: opt_linsys_fg


	implicit none

	integer, parameter 		:: ip = 4, sp = 4, dp = 8

	public opt_po, calc_objval, cleanup_PO

	private

	!--------------------------------------------------------
	! Factory pattern to return prox op.
	! Common interface for all prox ops.
	!--------------------------------------------------------
contains

	subroutine opt_po(po_ind, inform, objval, &
		   			  xopt, Drho, lenx, n, rho0, SD)
		!  !$acc routine seq

        integer, 	intent(in)				:: po_ind
		integer,	intent(out)				:: inform
		real(dp),	intent(out)				:: objval
        integer, 	intent(in)				:: lenx
		real(dp),	intent(out)				:: xopt(lenx)
        real(dp), 	intent(in)				:: Drho(lenx)
        real(dp), 	intent(in)				:: rho0
        real(dp), 	intent(in)				:: n(lenx)
		real(dp),	intent(in), optional	:: SD(lenx)

		integer						:: i,j,k

		select case(po_ind)
		case (PO_IND_BOUNDS)
			if (allocated(lb) .and. allocated(ub)) then
				call opt_boundss(inform, objval, xopt,		&
					Drho, lenx, n)
			else
				stop  "Must directly set lb, ub from po_bounds and l1weights from po_L1"
			endif

		case (PO_IND_BILINEAR)

			!!  $acc data copyout(n), copyin(xopt)
			!------------------------------------------------
			! These data transfers potentially problematic
			! Copy n from accel to CPU
			!  !$acc update self(n)
			!------------------------------------------------
			! 03 Feb 2019: if scaled, unscale n, and then send out scaled xopt
			!------------------------------------------------
			if (present(SD)) then
				! Unscale n
				call opt_bilinears(inform, objval, xopt, rho0, lenx, SD*n)
				! Scale x
                xopt = (1/SD)*xopt
			else
				call opt_bilinears(inform, objval, xopt, rho0, lenx, n)
			endif
			! 03 Feb 2019: and send out scaled xopt
			!------------------------------------------------
			! Copy xopt from CPU to accel
			!  !$acc update device(xopt)
			!------------------------------------------------
			!  !$acc end data

		case (PO_IND_L1)

			call opt_L1s(inform, objval, xopt, rho0, lenx, n)!, l1weight)

		case (PO_IND_CARD)

			call opt_cards(inform, objval, xopt, rho0, lenx, n)

		case (PO_IND_LINSYS)

			call opt_linsyss(inform, objval, xopt, Drho, lenx, n, SD=SD)

		case (PO_IND_LINSYS_MA57)

			call opt_linsys57(inform, objval, xopt, Drho, lenx, n, SD=SD)

		case (PO_IND_BOUNDS_L1)
			!------------------------------------------------
			! Combine bounds and L1 into same elemental PO
			!------------------------------------------------
			if (allocated(lb) .and. allocated(ub)) then
				call opt_bounds_L1s(inform, objval, xopt, Drho, lenx, n)
			else
				stop  "Must directly set lb, ub from po_bounds and l1weights from po_L1"
			endif

		case (PO_IND_BOUNDS_CARD)

			if (allocated(lb) .and. allocated(ub)) then
				call opt_bounds_cards(inform, objval, xopt, Drho, lenx, n)
			else
				stop  "Must directly set lb, ub from po_bounds and l0max from po_card"
			endif

		case (PO_IND_BILINEAR_MULT)
			! Copy n from accel to CPU
			!  !$acc update self(n)
			!------------------------------------------------
			! 03 Feb 2019: if scaled, unscale n, and then send out scaled xopt
			!------------------------------------------------
			if (present(SD)) then
				! Unscale n
				call opt_bilinear_multi(inform, objval, xopt, rho0, lenx, SD*n)
				! Scale x
                xopt = (1/SD)*xopt
			else
				call opt_bilinear_multi(inform, objval, xopt, rho0, lenx, n)
			endif
			! 03 Feb 2019: and send out scaled xopt
			!------------------------------------------------
			! Copy xopt from CPU to accel
			!  !$acc update device(xopt)

		case (PO_IND_LINSYS_FG)

			call opt_linsys_fg(inform, objval, xopt, Drho, lenx, n, SD=SD)
		

		case DEFAULT
			write (*, *) "proximal operator ind=",po_ind," not recognized!"
			write (*, *) "ind must be in [1:4, 13]"
!			error stop 'Unrecognized proximal operator ind'
			stop 'Unrecognized proximal operator ind'
		end select
	end subroutine

	function calc_objval(po_ind, xopt) result(objval)

		!----------------------------------------------------
		! Compute objval
		!----------------------------------------------------

        integer, 	intent(in)		:: po_ind
		real(dp),	intent(in)		:: xopt(:)
		real(dp) 	            	:: objval

		select case(po_ind)
		case (PO_IND_BOUNDS)
			objval = 0.0_dp

		case (PO_IND_BILINEAR)
			objval = 0.0_dp

		case (PO_IND_L1)
			!------------------------------------------------
			! objval = delta * sum(abs(x))
			!------------------------------------------------
			objval = calc_obj_L1s(xopt)

		case (PO_IND_BOUNDS_L1)
			!------------------------------------------------
			! objval = delta * sum(abs(x))
			!------------------------------------------------
			objval = calc_obj_L1s(xopt)

		case (PO_IND_CARD)
			objval = 0.0_dp

		case (PO_IND_BOUNDS_CARD)
			objval = 0.0_dp

		case (PO_IND_LINSYS)
			!------------------------------------------------
            ! objval = ||Ax - b||^2
			!------------------------------------------------
			objval = calc_obj_linsys(xopt, linsys_mF)

		case (PO_IND_LINSYS_MA57)
			objval = calc_obj_linsys57(xopt, linsys_mF)

		case (PO_IND_BILINEAR_MULT)
			objval = 0.0_dp

		case DEFAULT
			objval = 0.0_dp
			write (*, *) "proximal operator ind=",po_ind," not recognized!"
			write (*, *) "ind must be in [1:4, 13]"
!			error stop 'Unrecognized proximal operator ind'
			stop 'Unrecognized proximal operator ind'

		end select

	end function

	subroutine cleanup_PO(po_ind)

        integer, 	intent(in)		:: po_ind

		select case (po_ind)
		case (PO_IND_LINSYS)
			call cleanup_linsys()
		case default
			! Nothing to do
		end select

	end subroutine

end module
