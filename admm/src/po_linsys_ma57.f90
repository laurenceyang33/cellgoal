module po_linsys_ma57

	!--------------------------------------------------------
	! Proximal operator for linearly constrained least-squares
	! Solve linear system.
	!
	! [A11, A12; * [x; = [b1;
	! [A21, A22]    z]    b2]
	!
	! Laurence Yang, UCSD
	!
	! 14 Aug 2018: first version
	! 09 Oct 2018: error stop --> stop to support PGI compiler.
	!			   repeated kind parameters for PGI.
	! 17 Oct 2018: ma57 assume symmetric so provide only upper
	! 			   (or lower) triangular part without duplication.
	! 			   If duplicated, data will be wrong.
	!--------------------------------------------------------

	use po_linsys
	use sparsekit
	use sparseops,	only: hstack, vstack, csr2csc, sort_columns

	implicit none

	integer, parameter 		:: ip = 4, sp = 4, dp = 8

	real(dp)							:: cntl(5)
	integer(ip)							:: icntl(20)

	integer(ip)							:: ne
	integer(ip),	allocatable			:: irn(:)	!ne
	integer(ip),	allocatable			:: jcn(:)	!ne
	real(dp),		allocatable			:: vAsymm(:)!ne
	real(dp),		allocatable			:: vAsymm0(:)!ne ! need for adapt rho
	integer(ip),	allocatable			:: irwork(:)	
	integer(ip),	allocatable			:: jcwork(:)	
	real(dp),		allocatable			:: vawork(:)	

	integer(ip)							:: lkeep
	integer(ip),	allocatable			:: keep(:)	!lkeep
	integer(ip),	allocatable			:: iwork(:)	!5*n
	integer(ip)							:: infos(40)
	real(dp)							:: rinfo(20)
	integer(ip)							:: nrhs
	real(dp), 		allocatable			:: rhs(:,:)

	integer(ip)							:: lfact
	integer(ip)							:: lifact
    real(dp),		allocatable			:: fact(:)	! lfact
	integer(ip),	allocatable			:: ifact(:)	! lifact

	integer(ip)							:: job
	integer(ip)							:: lwork
	real(dp),		allocatable			:: work(:)	! lwork

	! Solution
	real(dp),		allocatable	   		:: x(:)
	real(dp),		allocatable			:: xy(:)	! x, y (lagrange mults)

	public	opt_linsys57, calc_obj_linsys57, solve_linsys57, update_rho_linsys57

	private

contains

	subroutine update_rho_linsys57(Drho)

		!----------------------------------------------------
		! Update factorization but not RHS (separated) 
		! with new rho 
		!
		! INPUTS
		! 	Drho 	vector of rho (diagonal) diffs
		!----------------------------------------------------

		real(dp),		intent(in)	 				:: Drho(:)
		integer(ip)									:: stat
		integer(ip)									:: inform
		integer										:: i,j,j1, k,k1
		integer(ip)									:: lenx

		!----------------------------------------------------
		! Update A11 and b1 blocks with rho (Drho--diagonal rho)
		! b1 block gets updated with each opt_linsys call.
		!----------------------------------------------------
		!   [2/k*F'F + Drho	C'
		!    C			 	0]
		! = [diagFF + Drho  C'
		!  	 C			    0]
		!
		! where diagFF is defined by make_PO_matrices(...) in
		! linsys_helper.f90
		!----------------------------------------------------
		! Symbolic factorization unchanged.
		! Update numeric factorization:
		!----------------------------------------------------
		do k=1, ne
			i = irn(k)
			j = jcn(k)
			! Only diagonal
			if (i==j) then
				! Just upper left block
				if (i <= mA11) then
					!****************************************
					! TODO: DEBUG 
					! 17 Oct 2018 [LY]: diagFF now scaled in precond_pos
					vAsymm(k) = diagFF(j) + Drho(j)
					! Could exit but don't in case i,j not sorted
				endif
			endif
		enddo

		! stat = umfpack_di_numeric(pAc0, iAc0, ac, Symbolic, Numeric)
		!call ma57bd(nA, Annz, vA, fact, lfact, ifact, lifact, lkeep, keep,	&
		call ma57bd(nA, ne, vAsymm, fact, lfact, ifact, lifact, lkeep, keep,	&
					iwork, icntl, cntl, infos, rinfo)
		stat = infos(1)

		! write(*,*) 'update_rho: New Drho =', minval(Drho),maxval(Drho)

	end subroutine

	subroutine opt_linsys57(inform, objval, xopt, Drho, lenx, n,	SD)

        !----------------------------------------------------
		! Subroutine version.
        !----------------------------------------------------
		! Solve linear system
		! [A11, A12; * [x; = [b1(n);
		! [A21, A22]    z]    b2]
		!
		! where A11, A12, A21, A22, b1, b2 are possibly
		! functions of rho.
		!
		! Submodules can determine exact implementation of
		! RHS and factorization updates depending on 
		! n and rho
        !----------------------------------------------------

		integer,	intent(out)			    :: inform
		real(dp),	intent(out)			    :: objval
		integer,	intent(in)			    :: lenx 
		real(dp),	intent(out)  	        :: xopt(lenx)

		real(dp),	intent(in)	 		    :: Drho(:)
		real(dp), 	intent(in)              :: n(lenx)
		real(dp),	intent(in),	optional    :: SD(lenx)
		real(dp) 	           				:: iSD(lenx)

		!----------------------------------------------------
		! Linear system-specific variables: 
		! Get from module variables.
		!----------------------------------------------------

		integer(ip) 							:: nnz
		integer(ip) 							:: nrow
		integer(ip)								:: mA1

		!----------------------------------------------------
		! Start with inform = 0 and then modify it 
		! as necessary by subsequent steps
		!----------------------------------------------------
		inform = 0

		!----------------------------------------------------
		! Is this a warm-start?
		! 27 Aug 2018 [LY]: only stack A and csrcsc once
		!----------------------------------------------------
		if (.not. warm_start_linsys) then
			call prepare_linsys(inform, lenx, mA1)
		endif

		!----------------------------------------------------
		! Update RHS with Proximal Point
		!----------------------------------------------------
		!****************************************************
		! 21 Sep 2018: scaling matrix
		!****************************************************
		if (present(SD)) then
			iSD = SD
		else
			iSD = 1.0_dp
		endif

		!****************************************************
		! 21 Sep 2018 [LY]: changed update model to
		! b1 = b1fix + Drho*N
		!****************************************************
		! b(lbound(b1,1):ubound(b1,1)) = b1 + SD*Drho*SD*n
		! b(ubound(b1,1)+1:) = b2
		
		! Need to allocate rhs to nA so it contains solution
		if (.not. allocated(rhs)) then
			allocate(rhs(nA,1))
			rhs = 0
		endif

		! rhs(lbound(b1,1):ubound(b1,1), 1) = b1 + iSD*Drho*iSD*n
		! 18 Oct 2018 [LY]: scaling before rho added, so
		! 					don't scale the rho
		rhs(lbound(b1,1):ubound(b1,1), 1) = b1 + Drho*n
		rhs(ubound(b1,1)+lbound(b2,1):ubound(b1,1)+ubound(b2,1), 1) = b2

		!----------------------------------------------------
		! Solve
		!----------------------------------------------------
		!----------------------------------------------------
		! 27 Aug 2018 [LY]: only stack A and csrcsc once
		!----------------------------------------------------
		! xy = solve_linsys(b, nA)
		if (.not. allocated(xy)) then
			allocate(xy(nA))
			xy = 0
		endif

		call solve_linsys57(inform, xy, nA)
		xopt = xy(1:lenx)

		! Can warm-start next iteration
		warm_start_linsys = .true.

		!****************************************************
		! TODO: compute objval
		!****************************************************
		objval = 0.0_dp
		! objval = calc_obj_linsys(xopt, linsys_mF)


	end subroutine

	subroutine solve_linsys57(inform, sol, ncol, fact_sym, fact_num, solve)

        !----------------------------------------------------
		! 30 Aug 2018 [LY]: subroutine version.
        !----------------------------------------------------

        !----------------------------------------------------
        ! Solve linear system: Ax = b
		!
		! All sparse matrix indices in CSR and 1-based indexing.
		! CSR converted to CSC and 0-based indexing before factorization.
		! 
		! Inputs
		!	nrow	# of rows
		!	ncol	# of columns
		!	nnz		number of nonzeros
		! 	iac		dimension(nrow+1) pointers to non-zero values
		! 	jac		dimension(nnz) column indices to non-zero values
		! 	ac		dimension(nnz) non-zero values
		! 	x		non-zero values
        !----------------------------------------------------

		implicit none

		!----------------------------------------------------
		! iac and jac changed to 0-based for umfpack
		!----------------------------------------------------
		integer(ip),	intent(inout)						:: inform
		integer(ip),	intent(in)							:: ncol
		real(dp),		intent(out)                         :: sol(ncol)   ! x, y (lagrange mult)

		logical,		optional							:: fact_sym
		logical,		optional							:: fact_num
		logical,		optional							:: solve
		logical                 							:: ifact_sym = .true.
		logical                 							:: ifact_num = .true.
		logical                 							:: isolve	 = .true.

		integer												:: i, j, k

		!--------------------------------------------------------
		! Factorize?
		if (present(fact_sym)) then
			ifact_sym = fact_sym
		else
			ifact_sym = mod_fact_sym
		endif
		if (present(fact_num)) then
			ifact_num = fact_num
		else
			ifact_num = mod_fact_num
		endif

		if (present(solve)) isolve = solve

		!--------------------------------------------------------
		! Factorize: keep cached
		!--------------------------------------------------------

		!----------------------------------------------------
		! Symbolic factorization
		!----------------------------------------------------
		if (ifact_sym .and. inform>=0) then

			! Default options
			call ma57id(cntl, icntl)

			! Fill in row indices from row pointers of CSR format
			lkeep = 5*nA + Annz + max(nA,Annz) + 42 + 2*nA

			if (.not. allocated(keep)) allocate(keep(lkeep))
			if (.not. allocated(iwork)) allocate(iwork(5*nA))

			!************************************************
			! TODO: assumes symmetric (indefinite) so only
			! should provide aij once (not twice with aji too).
			! In fact, it seems that duplicating entries
			! accumulates the aij, so matrix is wrong?
			!************************************************
			if ( (.not. allocated(irn)) .or. (.not. allocated(jcn)) ) then
				allocate(irwork(Annz))
				allocate(jcwork(Annz))
				allocate(vawork(Annz))
				! From CSR to symmetric.  Just keep the upper triangular.
				! Count size and fill in values.
				ne = 0
				do i=1,mA
					do k=pA(i), pA(i+1)-1
						j = jA(k)
						if (j >= i) then
							ne = ne + 1
							irwork(ne) = i
							jcwork(ne) = j
							vawork(ne) = vA(k)
						endif
					enddo
				enddo
				! Finally, reallocate to the actual size and fill in
				allocate(irn(ne))
				allocate(jcn(ne))
				allocate(vAsymm(ne))

				irn = irwork(1:ne)
				jcn = jcwork(1:ne)
				vAsymm = vawork(1:ne)

				! Don't need the work arrays anymore
				deallocate(irwork)
				deallocate(jcwork)
				deallocate(vawork)
			endif

			! stat = umfpack_di_symbolic(mA, nA, pAc0, iAc0, ac, Symbolic)
			! call ma57ad(nA, Annz, irn, jA, lkeep, keep, iwork, icntl, infos, rinfo)
			call ma57ad(nA, ne, irn, jcn, lkeep, keep, iwork, icntl, infos, rinfo)

			stat = infos(1)
			if (stat < 0) then
				inform = stat
			endif
			mod_fact_sym = .false.
		endif

		!----------------------------------------------------
		! Numeric factorization
		!----------------------------------------------------
		if (ifact_num .and. inform>=0) then

			! stat = umfpack_di_numeric(pAc0, iAc0, ac, Symbolic, Numeric)
			lfact 	= 2*infos(9)
			lifact 	= 2*infos(10)
			!************************************************
			! Should also check that size hasn't changed
			!************************************************
			if (.not. allocated(fact)) allocate(fact(lfact))
			if (.not. allocated(ifact)) allocate(ifact(lifact))

			! call ma57bd(nA, Annz, vA, fact, lfact, ifact, lifact, lkeep, keep,	&
			call ma57bd(nA, ne, vAsymm, fact, lfact, ifact, lifact, lkeep, keep,	&
						iwork, icntl, cntl, infos, rinfo)

			stat = infos(1)
			if (stat < 0) then
				inform = stat
			endif
			mod_fact_num = .false.
		endif

		!--------------------------------------------------------
		! SOLVE
		!--------------------------------------------------------
		if (isolve .and. inform>=0) then
			sol		= 0.0

			! stat = umfpack_di_solve(UMFPACK_A, pAc0, iAc0, ac, sol, b, Numeric)
            lwork = nA
			if (.not. allocated(work)) allocate(work(lwork))

			! Upon exit, rhs contains the solution, so need nA x 1, not mA x 1
			job = 1	
			call ma57cd(job, nA, fact, lfact, ifact, lifact, 1, rhs, mA, &
				work, lwork, iwork, icntl, infos)

			sol = rhs(:,1)

			stat = infos(1)
			if (stat < 0) then
				inform = stat
			endif
		endif

	end subroutine

	function calc_obj_linsys57(xopt, mF) result(objval)

		!----------------------------------------------------
		! objval = ||Fx - e||^2
		!----------------------------------------------------

		real(dp), 		intent(in)				:: xopt(:)
		integer(ip), 	intent(in)				:: mF
		real(dp)								:: Fx(mF)
		real(dp) 	           					:: objval

		! TODO (low priority since calc objval once at end):
		! time speed of amux vs sp_blas
		call amux(mF, xopt, Fx, linsys_vF, linsys_jF, linsys_pF)

        objval = sum((Fx - linsys_e)**2)

	end function

end module
