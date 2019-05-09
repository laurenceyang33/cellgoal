module po_bounds_card

	!--------------------------------------------------------
	! Proximal operator for:
	! arg min  rho/2*||x-n||2 subject to l <= x <= u
	!	   x
	! xopt = Projection{l,u}(x) 
	!
	! And cardinality constraint: hard thresholding.
	!
	! Laurence Yang, UCSD
	!--------------------------------------------------------
	! 09 Oct 2018: error stop --> stop to support PGI compiler
	! 08 Nov 2018: port from po_bounds_L1 to cardinality.
	!--------------------------------------------------------

	use po_bounds,	only: lb, ub
	use po_card,	only: l0max, card_threshold

	implicit none

	integer, 	parameter					:: ip=4, dp=8

	public opt_bounds_cards

	private

contains

	subroutine opt_bounds_cards(inform, objval, xopt,		&
			Drho, lenx, n)

		integer, 	intent(in)	 				:: lenx
		real(dp),	intent(in)	 				:: Drho(lenx)
		real(dp), 	intent(in), dimension(lenx)	:: n

		real(dp), 	intent(out),dimension(lenx)	:: xopt
		integer,	intent(out)					:: inform
		real(dp),	intent(out)					:: objval

		!  !$acc kernels
		xopt = min(max(card_threshold(n,l0max,lenx), lb), ub)
		!  !$acc end kernels

        objval = 0.0_dp

		inform = 0

	end subroutine 

end module
