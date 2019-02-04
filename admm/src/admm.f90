!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File admm_acc.f90: 
! 
! ADMM 
!
! Laurence Yang, UCSD
!
! 13 Aug 2018: first version. 
! 30 Aug 2018: modified from admmco that uses procedural pointer for 
! 			   PO factory. Instead using subroutine.
! 09 Oct 2018: error stop --> stop for compatibility with PGI.
! 04 Feb 2019: commented out all the ACC directives for now
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module ADMM

	implicit none

	public admms

	private

	integer,	parameter					:: ip=4, dp=8, sp=4
	real(dp),  	parameter					:: zero_rho = 1.0d-14

contains

	subroutine admms(	inform, objval, iter,				&
						Z, lenz, X, nprox, 					& 
						rho, gam, max_iter, 				&
						po_types, 							&
						abs_conv_in, print_freq_in,			&
						verbose, po_rhos, Uopt,				&
						Nopt, SD, adapt_rho, 				&
						rho_min, rho_max, 					&
						ratio_thresh, adapt_rho_freq,		&
						res_thresh, adapt_fold,				&
						prim_res_out, dual_res_out,			&
						stepsize, adapt_step,				& 
						adapt_step_freq,  					&
						step_thresh, step_fold,				&
						step_min, step_max, unit_iters,		&
						stall_tol, nCond)

		!--------------------------------------------------------
		! INPUTS:
		! 	- max_iter
		! 	- nprox		: number of proximal operators
		! 	- lenz		: size of Z
		! 	- Z	(inout)	: consensus variable
		! 	- X	(inout)	: proximal operator variable
		!	- SD		: diagonal scaling as array
        !
		! OUTPUTS:
		!	- Uopt
		!	- Nopt
		! 	
		! 	[optional]
		!
		!--------------------------------------------------------

		use pos_factory,	only: opt_po, calc_objval, cleanup_PO
		use proxop,			only: update_multipliers

		real(dp)   											:: tic,toc,t_elapsed
		integer												:: i,j,iprox,jstart,jend

		integer,	intent(in)								:: max_iter
		integer,	intent(in)								:: nprox	! Number of prox ops
		integer, 	intent(in)								:: lenz
		real(dp),                 	dimension(lenz,nprox)	:: Zmat
		real(dp),  	intent(inout),	dimension(lenz)			:: Z
		real(dp),  	              	dimension(lenz)			:: Z0
		real(dp),	intent(inout),  dimension(lenz,nprox)	:: X
		real(dp),  					dimension(lenz,nprox)	:: M
		real(dp)                	        				:: U(lenz,nprox)
		real(dp)                	          				:: N(lenz,nprox)
		real(dp),  	intent(out),	optional				:: Uopt(lenz,nprox)
		real(dp),	intent(out),	optional  				:: Nopt(lenz,nprox)
		real(dp),	intent(in),		optional				:: SD(lenz)
		real(dp)                            				:: iSD(lenz)

		real(dp),  					dimension(lenz,nprox)	:: M0
		real(dp),  					dimension(lenz,nprox)	:: N0
		real(dp),					dimension(lenz)			:: xopt
		real(dp),	intent(in)								:: gam
		real(dp),	intent(in)								:: rho
		real(dp),					dimension(lenz)			:: Drho
		real(dp),					dimension(nprox)		:: rhos
		real(dp),	optional,		dimension(nprox)		:: po_rhos
		real(dp)											:: m_change	! Message change
		real(dp)											:: n_change ! Message change
		real(dp)											:: mn_change ! Message change
		real(dp),					dimension(nprox)		:: mn_changes
		real(dp)											:: prim_res
		real(dp)											:: dual_res
		real(dp),	optional								:: nCond
		real(dp)											:: inCond

		real(dp),	intent(out),	optional				:: prim_res_out
		real(dp),	intent(out),	optional				:: dual_res_out

		integer,	intent(in),		dimension(nprox)		:: po_types

		real(dp)											:: abs_conv
		real(dp),	optional								:: abs_conv_in
		character(len=100)									:: format_iter
		character(len=100)									:: format_header
		integer             								:: print_freq
		integer,	optional								:: print_freq_in
		logical,	optional								:: verbose
		logical												:: iverbose
		integer,	optional								:: unit_iters
		integer												:: iunit_iters
		real(dp),	optional								:: stall_tol
		real(dp)											:: istall_tol

		! For adaptive rhos
		integer,	optional								:: adapt_rho
		integer             								:: iadapt_rho
		real(dp),	optional								:: rho_min, rho_max
		real(dp)            								:: irho_min, irho_max
		real(dp),	optional								:: ratio_thresh
		real(dp)            								:: iratio_thresh
		integer,	optional								:: adapt_rho_freq
		integer            									:: iadapt_rho_freq
		real(dp),	optional								:: res_thresh
		real(dp)            								:: ires_thresh
		real(dp),	optional								:: adapt_fold
		real(dp)            								:: iadapt_fold

		real(dp)											:: objvals(nprox)
		
		! For adaptive multiplier update stepsize
		real(dp),	intent(in),optional						:: stepsize
		real(dp)											:: istepsize
		integer,	intent(in),optional						:: adapt_step
		integer												:: iadapt_step
		integer,	intent(in),optional						:: adapt_step_freq
		integer												:: iadapt_step_freq
		real(dp),	intent(in),optional						:: step_thresh
		real(dp)											:: istep_thresh
		real(dp),	intent(in),optional						:: step_fold
		real(dp)											:: istep_fold
		real(dp),	intent(in),optional						:: step_min, step_max
		real(dp)											:: istep_min, istep_max 

		!----------------------------------------------------
		! Status output
		!----------------------------------------------------
		integer,	intent(out)								:: inform
		real(dp),	intent(out)								:: objval
		integer,	intent(out)								:: iter

		!--------------------------------------------------------
		! Initialize
		!--------------------------------------------------------
		if (present(abs_conv_in)) then
			abs_conv = abs_conv_in
		else
			abs_conv = 1.0d-9
		endif
		if (present(print_freq_in)) then
			print_freq = print_freq_in
		else
			print_freq = 10
		endif
		if (present(verbose)) then
			iverbose = verbose
		else
			iverbose = .true.
		endif
		if (present(po_rhos)) then
			rhos = po_rhos
		else
			rhos = rho
		endif

		if (present(SD)) then
			iSD = SD
		else
			iSD = 1.0_dp
		endif

		if (present(adapt_rho)) then
			iadapt_rho = adapt_rho
		else
			iadapt_rho = 0
		endif

		if (present(rho_min)) then
			irho_min = rho_min
		else
			irho_min = 1.0d-3
		endif

		if (present(rho_max)) then
			irho_max = rho_max
		else
			irho_max = 1.0d+3
		endif

		if (present(ratio_thresh)) then
			iratio_thresh = ratio_thresh
		else
			iratio_thresh = 2.0
		endif

		if (present(adapt_rho_freq)) then
			iadapt_rho_freq = adapt_rho_freq
		else
			iadapt_rho_freq = 100
		endif

		if (present(res_thresh)) then
			ires_thresh = res_thresh
		else
			ires_thresh = 1e-6
		endif

		if (present(adapt_fold)) then
			iadapt_fold = adapt_fold
		else
			iadapt_fold = 2.0
		endif

		! Adaptive step size
		if (present(stepsize)) then
			istepsize = stepsize
		else
			istepsize = 1.0
		endif

		if (present(adapt_step)) then
			iadapt_step = adapt_step
		else
			iadapt_step = 0
		endif

		if (present(adapt_step_freq)) then
			iadapt_step_freq = adapt_step_freq
		else
			iadapt_step_freq = max_iter
		endif

		if (present(step_thresh)) then
			istep_thresh = step_thresh
		else
			istep_thresh = max_iter
		endif

		if (present(step_fold)) then
			istep_fold = step_fold
		else
			istep_fold = 0.1
		endif

		if (present(step_max)) then
			istep_max = step_max
		else
			istep_max = 1.0
		endif

		if (present(step_min)) then
			istep_min = step_min
		else
			istep_min = 0.01
		endif

		if (present(unit_iters)) then
			iunit_iters = unit_iters
		else
			iunit_iters = 0
		endif

		if (present(stall_tol)) then
			istall_tol = stall_tol
		else
			istall_tol = 1e-3*abs_conv
		endif

		if (present(nCond)) then
			inCond = nCond
		else
			inCond = 1.0
		endif

		format_header = '(1x, A15, 8(",", A15))'
		format_iter   = '(1x, I15, 8(",", E15.4))'

		Z0 = Z
		
		objval = 0.0_dp
		!****************************************************
		! Use initial values to get initial msgs
		U		= 0.0_dp
		N		= X
		M 		= 0.0_dp ! X
		M0		= 0.0_dp ! M
		N0		= 0.0_dp ! N

		!--------------------------------------------------------
		! MAIN LOOP
		!--------------------------------------------------------

        ! Print stats header
		if (iverbose) then
			write(*, format_header) 'Iter','res_primal','res_dual','msg_change',	&
				'rho(1)', 'min rho', 'max rho', 'stepsize','time (s)'
			write(*,'(A)') '#--------------------------------------------------------------------------'
		endif
		if (iunit_iters .ne. 0) then
			write(iunit_iters, format_header) 'Iter','res_primal','res_dual','msg_change',	&
				'rho(1)', 'min rho', 'max rho', 'stepsize','time (s)'
			write(iunit_iters,'(A)') '#--------------------------------------------------------------------------'
		endif

		Drho = rho

		call cpu_time(tic)

		!----------------------------------------------------
		! Accelerator data management
		!----------------------------------------------------
		!  !$acc	data copy(Z,Z0,mn_change)					&
		!  !$acc 	copyin(	po_types,rhos,X,U,M,N,M0,N0,		&
		!  !$acc			abs_conv,lenz,iSD)

		do iter = 1, max_iter
			!------------------------------------------------
			! X update: optimize each PO
			!------------------------------------------------
			!  !$acc data present(X,po_types,lenz,N,rhos,iSD)
			!  !$acc parallel loop private(xopt,Drho)
			do iprox = 1, nprox
				Drho = rhos(iprox)
				if (rhos(iprox) >= zero_rho) then
					call opt_po(po_types(iprox), inform, objval, 		&
							xopt, Drho, lenz, N(:,iprox), rhos(iprox),	&
							iSD)

					if (inform < 0) then
						write(*,*) '#PO',iprox,'failed. Exiting.'
						exit
					endif
					X(:,iprox) = xopt
				endif
			end do
			!  !$acc end data

			!************************************************
			! If any PO unsuccessful, exit right away
			!************************************************
			if (inform < 0) then
				write(*,*) '#At least one PO failed. Exiting.'
				exit
			endif

			!------------------------------------------------
			! M update: M  = X + U
			!------------------------------------------------
			M0 = M
			M = X + U
			!************************************************
			! TODO: BENCHMARK
			do j=1, nprox
				if (rhos(j) < zero_rho) then
					M(:,j) = 0.0_dp
				endif
			enddo
			!************************************************

			!------------------------------------------------
			! Z update: Zj = (1-gamma)*Z0j + (sum_i rhoij*(gamma*Xij + Uij) / sum_i rhoij)  
			! Z = (1-gam)*Z + rhos*(gam*X + U) / sum(rhos)
			!------------------------------------------------
			! Save previous Z
			Z0 = Z
			! Get the consensus by summing over prox operators

			!  !$acc data copyin(X,U,rhos) copyout(Z)
			call update_Z_fast(Z, X, U, rhos, gam)
			!  !$acc end data

			! Just fill in the copies since all the same
			Zmat = spread(Z,2,nprox)

			!------------------------------------------------
			! U update: U  = U0 + gamma*X - Z + (1-gamma)*Z0
			!------------------------------------------------
			!  !$acc data present(U,X,Z,Z0,gam,istepsize)
			do j = 1, nprox !zcol
				call update_multipliers(lenz, U(:,j), X(:,j), Z, Z0, gam, istepsize)
			end do
			!  !$acc end data

			!------------------------------------------------
			! N update: N  = Z - U
			!------------------------------------------------
			N0 = N
			! N = Z - U
			do j=1,nprox
				!********************************************
				! TODO: BENCHMARK
				!********************************************
				if (rhos(j) >= zero_rho) then
					N(:,j) = Z - U(:,j)
				else
					N(:,j) = 0.0_dp
				endif
			enddo

			!------------------------------------------------
			! Check convergence
			!------------------------------------------------
			mn_change  = 0.0_dp
			mn_changes = 0.0_dp

			do j = 1, nprox !zcol
				if (rhos(j) >= zero_rho) then
					m_change = msg_change(M(:,j), M0(:,j))
					n_change = msg_change(N(:,j), N0(:,j))
				else
					m_change = 0.0_dp
					n_change = 0.0_dp
				endif
				mn_changes(j) = m_change + n_change
			end do
            
			! Update thet new message change
			mn_change = maxval(mn_changes)

			if ( (mn_change) <= abs_conv) then  
				call primal_residual(prim_res, X, Zmat, rhos)
				call dual_residual(dual_res, Z, Z0, sum(rhos)/nprox)
				call cpu_time(toc)
				t_elapsed = toc-tic

				if (iverbose) then
					write(*,format_iter) iter, prim_res, dual_res, mn_change, 	&
						rhos(1), minval(rhos), maxval(rhos), istepsize, t_elapsed 
					if (mn_change <= abs_conv) then
						write(*,'(a, E15.8)') '#Converged with message change: ', mn_change
					else
						write(*,'(a, E15.8)') '#Stalled with message change: ', mn_change
					endif
				endif
				if (iunit_iters .ne. 0) then
					write(iunit_iters,format_iter) iter, prim_res, dual_res, mn_change, 	&
						rhos(1), minval(rhos), maxval(rhos), istepsize, t_elapsed 
					if (mn_change <= abs_conv) then
						write(iunit_iters,'(a, E15.8)') '#Converged with message change: ', mn_change
					else
						write(iunit_iters,'(a, E15.8)') '#Stalled with message change: ', mn_change
					endif
				endif

				!============================================
				exit
				!============================================
			endif

			!------------------------------------------------
			! Print stats 
			!------------------------------------------------
			if (iter==1 .or. modulo(iter, print_freq) == 0) then
				call primal_residual(prim_res, X, Zmat, rhos)
				call dual_residual(dual_res, Z, Z0, sum(rhos)/nprox)
				call cpu_time(toc)
				t_elapsed = toc-tic
				if (iverbose) then
					write(*,format_iter) iter, prim_res, dual_res, mn_change, 	&
						rhos(1), minval(rhos), maxval(rhos), istepsize, t_elapsed 
!						maxval(abs(X)), maxval(abs(Z)), maxval(abs(U)), t_elapsed 
				endif
				if (iunit_iters .ne. 0) then
					write(iunit_iters,format_iter) iter, prim_res, dual_res, mn_change, 	&
						rhos(1), minval(rhos), maxval(rhos), istepsize, t_elapsed 
				endif
			endif

			!------------------------------------------------
			! Adaptive rho
			!------------------------------------------------
			if (iadapt_rho  > 0) then
				if (modulo(iter, iadapt_rho_freq) == 0) then
					call primal_residual(prim_res, X, Zmat, rhos)
					call dual_residual(dual_res, Z, Z0, sum(rhos)/nprox)
					call update_rhos(rhos, iadapt_rho, po_types, 	&
						X, Z, U, irho_min, irho_max,				&
						prim_res, dual_res, iratio_thresh,			&
						ires_thresh, iadapt_fold, inCond)
				endif
			endif

			!------------------------------------------------
			! Adaptive step-size alpha
			!------------------------------------------------
			if (iadapt_step > 0) then
				if (modulo(iter, iadapt_step_freq) == 0) then
					if (mn_change <= istep_thresh) then
						call update_stepsize(istepsize, istep_fold, istep_min, istep_max)
					endif
				endif
			endif

		!--------------------------------------------------------
		! END MAIN LOOP
		end do
		! !$acc end data
		!--------------------------------------------------------

		!--------------------------------------------
		! Only need to check objval at the end
		!--------------------------------------------
		do iprox=1, nprox
			objvals(iprox) = calc_objval(po_types(iprox), Z)
		enddo
		objval = sum(objvals)
		if (iverbose) write(*,*) '#objvals:',objvals
		if (iunit_iters .ne. 0) write(iunit_iters,*) '#objvals:',objvals 

		!--------------------------------------------
		! Populate optional outputs
		!--------------------------------------------
		if (present(Uopt)) Uopt = U
		if (present(Nopt)) Nopt = N

		if (present(prim_res_out)) prim_res_out = prim_res
		if (present(dual_res_out)) dual_res_out = dual_res

		!--------------------------------------------
		! Clean up: free C pointers, etc.
		!--------------------------------------------
		do iprox=1, nprox
			call cleanup_PO(po_types(iprox))
		enddo

	end subroutine admms

	subroutine primal_residual(res, x, z, rhos)

		real(dp), 	intent(inout) 							:: res
		real(dp), 	intent(in), 	dimension(:,:) 			:: x
		real(dp), 	intent(in), 	dimension(:,:) 			:: z
		real(dp),	intent(in),		dimension(:)			:: rhos

		!----------------------------------------------------
		! If rhoj = 0, exclude from residual calc
		!----------------------------------------------------
		res = maxval(abs(x-z)*spread(merge(0.0_dp, 1.0_dp, rhos==0),1,size(x,1)))

	end subroutine primal_residual

	subroutine dual_residual(res, z, z0, rho)

		real(dp),	intent(inout)  							:: res
		real(dp),	intent(in),    	dimension(:)			:: z
		real(dp),	intent(in), 	dimension(:)			:: z0
		real(dp), 	intent(in)								:: rho

		res = maxval(abs(rho*(z-z0)))

	end subroutine dual_residual

	real(dp) function msg_change(msg,msg0) result(dmsg)

		real(dp),	intent(in)						:: msg(:)
		real(dp),	intent(in)						:: msg0(:)

		!  !$acc kernels
		dmsg = maxval(abs(msg-msg0))
		!  !$acc end kernels

	end function

	subroutine update_Z_fast(z, X, U, rhos, g)

		!----------------------------------------------------
		! Z = (1-gam)*Z + rhos*(gam*X + U) / sum(rhos)
		!----------------------------------------------------

		integer										:: i,j
		integer										:: nx, nprox
		real(dp),	intent(in)						:: g	! gamma
		real(dp),	intent(in)		              	:: X(:,:)
		real(dp),	intent(in)		              	:: U(:,:)
		real(dp),	intent(in)		                :: rhos(:)
		real(dp),	intent(inout)	                :: z(:)
		real(dp)	                              	:: rho_norm(size(X,1),size(X,2))
		real(dp)									:: denom, zi

		!----------------------------------------------------
		! Serial version:
		! rho_norm = spread(rhos/sum(rhos),1,size(z))
		! z = (1-g)*z + sum(rho_norm * (g*X + U), 2) 
		!----------------------------------------------------

		!----------------------------------------------------
		! ACC version
		denom = sum(rhos)

		!  !$acc data present(rhos,X,U)
		!  !$acc parallel loop private(zi)
		do i=1, size(X,1)
			zi = 0
			!  !$acc loop reduction(+:zi)
			do j=1, size(X,2)
				zi = zi + rhos(j) * (g*X(i,j) + U(i,j))
			enddo
			z(i) = (1-g)*z(i) + zi / denom
		enddo
		!  !$acc end data

		!----------------------------------------------------
			
	end subroutine

	subroutine update_rhos(	rhos_new, adapt_rule, po_types, X, Z, U, 	&
							rho_min, rho_max, prim_res, dual_res, 		&
							ratio_thresh, res_thresh, adapt_fold, nCond)

		!------------------------------------------------
		! Update rho based on adapt_rule
		! 0 : no adaptive rule
		! 1 : residual balancing, increase or decrease by fixed ratio
		! 2 : residual balancing, computed ratio
		!------------------------------------------------

		use adaptive_rhos,	only: res_bal_scaled, res_bal_scaled1, res_bal1
		use proxop_twa,		only: update_PO_rho

		integer,	intent(in)						:: adapt_rule
		integer,	intent(in)						:: po_types(:)
		real(dp),	intent(inout)					:: rhos_new(:)
		real(dp)                 					:: rhos0(size(rhos_new))
		real(dp),	intent(in)						:: X(:,:), U(:,:)
		real(dp),	intent(in)						:: z(:)

		real(dp),	intent(in)						:: rho_min, rho_max
		logical										:: rho_changed
		real(dp),	intent(in)						:: prim_res, dual_res
		real(dp),	intent(in)						:: ratio_thresh
		real(dp),	intent(in)						:: res_thresh
		real(dp),	intent(in)						:: adapt_fold
		real(dp),	intent(in)						:: nCond

		integer										:: i,j,k

		if (adapt_rule > 0) then
			if ((prim_res < res_thresh) .or. (dual_res < res_thresh)) then
				rhos0 = rhos_new
                !------------------------------------------------------------
				! Loop through proxops
				! Each ProxOp has potentially different adaptive rule
                !------------------------------------------------------------
				do i=1, size(rhos_new)
					if (rhos_new(i) > 0) then
						select case(adapt_rule)
						case (1)
							rhos_new(i) = res_bal1(rhos_new(i), 			&
								prim_res, dual_res, ratio_thresh,			&
								rho_min, rho_max, adapt_fold)
						case (2)
							rhos_new(i) =  res_bal_scaled1(rhos_new(i),		&
								X(:,i), z, U(:,i), prim_res, dual_res, 		&
								ratio_thresh, nCond, po_types(i), 			&
								rho_min, rho_max)
						case default
							write(*,*) '#Received adapt_rule:',adapt_rule
							stop 'Unknown adapt_rule'
						end select

						rho_changed = abs(rhos_new(i) - rhos0(i))/rhos0(i) > 0.1

						if (rho_changed) then
							!----------------------------------------
							! Need to update LHS of linsys
							!----------------------------------------
							call update_PO_rho(spread(rhos_new(i),1,size(z)), po_types(i))
						endif
					endif
				enddo
                !------------------------------------------------------------
			endif
		endif

	end subroutine

	subroutine update_stepsize(stepsize, step_fold, step_min, step_max)

		!----------------------------------------------------
		! Update step-size for multiplier updates
		! 0 : no adaptive rule
		! 1 : multiply by fixed factor
		! 2 : line search
		! ...
		!----------------------------------------------------

		real(dp),	intent(inout)					:: stepsize
		real(dp),	intent(in)						:: step_fold
		real(dp),	intent(in)						:: step_min
		real(dp),	intent(in)						:: step_max

		stepsize = max(step_min, min(step_max, stepsize*step_fold))

	end subroutine
	
	! subroutine linesearch_stepsize(stepsize, step_min, step_max, Z, Z0)

	! 	!----------------------------------------------------
	! 	! Line search adaptive step size.
	!	! Now, wouldn't that be nice...?
	! 	!----------------------------------------------------

	! 	real(dp),	intent(inout)					:: stepsize
	! 	real(dp),	intent(in)						:: step_fold
	! 	real(dp),	intent(in)						:: step_min
	! 	real(dp),	intent(in)						:: step_max

	! end subroutine

end module
