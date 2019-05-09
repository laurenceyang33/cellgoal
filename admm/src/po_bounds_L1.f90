module po_bounds_L1

	!--------------------------------------------------------
	! Proximal operator for:
	! arg min  rho/2*||x-n||2 subject to l <= x <= u
	!	   x
	! xopt = Projection{l,u}(x) 
	!
	! And L1
	!
	! Laurence Yang, UCSD
	!--------------------------------------------------------
	! 09 Oct 2018: error stop --> stop to support PGI compiler
	!--------------------------------------------------------

	use po_L1,		only: l1weight, soft_threshold, l1weights
	use po_bounds,	only: lb, ub

	implicit none

	integer, 	parameter					:: ip=4, dp=8

	public opt_bounds_L1, opt_bounds_L1s, opt_bounds_L1s_par

	private

contains

	elemental real(dp) function threshold_project(v,d,l,u) result(xopt)
		real(dp), 	intent(in)	:: v,d,l,u
		if (abs(v) <= d) then
			xopt = 0
		else
			xopt = sign(1.0d+0,v)*(abs(v)-d)
		endif

		xopt = min(max(xopt, l), u)
	end function

	function opt_bounds_L1(rho, lenx, n, real_args, int_args,  bool_args) result(xopt)
		!  !$acc routine seq

		real(dp),	intent(in)	 				:: rho
		integer, 	intent(in)	 				:: lenx
		real(dp), 	intent(in), dimension(lenx)	:: n
		real(dp), 	optional, 	dimension(:,:)	:: real_args
		integer, 	optional, 	dimension(:,:)	:: int_args
		logical, 	optional, 	dimension(:,:)	:: bool_args
		real(dp),	            dimension(lenx)	:: l
		real(dp),	           	dimension(lenx)	:: u
		real(dp)	           	               	:: d

		real(dp),				dimension(lenx)	:: xopt
		integer				  				  	:: j

		if (present(real_args)) then
			l = real_args(:,1)
			u = real_args(:,2)
			d = real_args(1,3)
		else
!			error stop "Must have l,u in real_args(:,1:2) and d in real_args(1,3)"
			stop "Must have l,u in real_args(:,1:2) and d in real_args(1,3)"
		endif

		! xopt = min(max(soft_threshold(n,d), l), u)
		! Should project into bounds, then do soft_threshold
		! Should do soft_threshold and then project into bounds if infeas
		! xopt = soft_threshold(min(max(n,l),u), d/rho)
		xopt = min(max(soft_threshold(n,d/rho), l), u)

	end function

	subroutine opt_bounds_L1s(inform, objval, xopt,		&
			Drho, lenx, n)
		!  !$acc routine seq

		integer, 	intent(in)	 				:: lenx
		real(dp),	intent(in)	 				:: Drho(lenx)
		real(dp), 	intent(in), dimension(lenx)	:: n

		real(dp), 	intent(out),dimension(lenx)	:: xopt
		integer,	intent(out)					:: inform
		real(dp),	intent(out)					:: objval

		!  !$acc kernels
		xopt = min(max(soft_threshold(n, l1weights/Drho), lb), ub)
		!  !$acc end kernels

        objval = 0.0_dp

		inform = 0

	end subroutine 

	function soft_threshold1(v,d) result(x)
		!  !$acc routine seq

		real(dp), 	intent(in)	:: v
		real(dp),	intent(in)	:: d
		real(dp)				:: x

		integer					:: i,j

		if (abs(v) <= d) then
			x = 0
		else
			x = sign(1.0d+0,v)*(abs(v)-d)
		endif

	end function soft_threshold1

	subroutine opt_bounds_L1s_par(inform, objval, xopt,		&
			Drho, lenx, n)

		!----------------------------------------------------
		! Parallel version
		!----------------------------------------------------
		integer,	intent(out)					:: inform
		real(dp),	intent(out)					:: objval
		integer, 	intent(in)	 				:: lenx
		real(dp), 	intent(in), dimension(lenx)	:: n
		real(dp), 	intent(out),dimension(lenx)	:: xopt
		real(dp),	intent(in)	 				:: Drho(lenx)
		real(dp)								:: xj

		integer									:: i,j,k


		!  !$acc data present(n,l1weights,Drho,lb,ub)
		!  !$acc parallel loop
		do j=1,lenx
			xj = min(max(soft_threshold1(n(j), l1weights(j)/Drho(j)), lb(j)), ub(j))
			xopt(j) = xj
		enddo
		!  !$acc end data

	end subroutine

end module po_bounds_L1
