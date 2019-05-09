!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File accfits57.f90: 
!
! Flux fit problem for any GEM using fine-grained ADMM.
!
! Laurence Yang, UCSD
!
! 12 Oct 2018: first version. 
! 16 Oct 2018: option to use MA57 linear solver
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program accfits57

	use admm_acc,	only: admms
	use po_constants

	use mm_helper,	only: csr_from_mm, array_from_mm, matrix_from_mm
	use mm_io,		only: mm_file_write
	use sparsekit
	use sparseops
	use po_linsys

	use linsys_helper,	only: make_PO_matrices
	use po_bounds,		only: lb, ub
	use po_L1,			only: l1weights
	use po_card,		only: l0max
	use po_bilinear_multi,	only: ds, csenses, nxs, xstarts, xends, ystarts, yends
	use precondition,	only: SD, SE, unscale_ruiz
	use precond_pos,	only: scale_po, unscale_po

	! For testing
	use proxop,			only: update_multipliers
	use pos_factory,	only: calc_objval

	implicit none

	!--------------------------------------------------------
	integer,	parameter				:: ip = 4, dp = 8, sp=4
	integer								:: i,j,max_iter,print_freq
	logical            					:: verbose=.true.
	character(*), parameter				:: version = '0.1'

	!--------------------------------------------------------
	! Command arguments
	!--------------------------------------------------------
	character(100)						:: specs_file 	! = 'options.spec'
	integer(ip)							:: unit_specs
	integer(ip)							:: ios
	character(80)						:: pname
	real(dp)							:: pval
    integer								:: nargs
	character(250)              		:: arg

	!--------------------------------------------------------
	! Output file
	!--------------------------------------------------------
	character(150)						:: out_file
	integer(ip)							:: unit_out
	integer(ip)							:: ival(1)
	real(sp)							:: rval(1)
	complex(sp)							:: cval(1)
	integer(ip)							:: indx(1), jndx(1)
	character(7)						:: field
	character(14)						:: id
	character(10)						:: rep
	character(19)						:: symm
	character(6)						:: type

	integer								:: inform
	real(dp)							:: objval
	integer								:: iter

	!--------------------------------------------------------
	! Stats
	!--------------------------------------------------------
	real(sp)   							:: tic,toc,t_elapsed

	!--------------------------------------------------------
	! Problem data
	!--------------------------------------------------------
	integer              				:: lenx  ! = 429
	integer,	parameter				:: nprox = 3

	integer,	parameter				:: IND_LIN  = 1
	integer,	parameter				:: IND_BND_REG  = 2
	integer,	parameter				:: IND_BI 	= 3

	real(dp)							:: rho
	real(dp)							:: gam
	integer,	dimension(nprox)   		:: po_types
	real(dp),	dimension(nprox)   		:: po_rhos

	real(dp)							:: rho_bi
	real(dp)							:: rho_linsys
	real(dp)							:: rho_bnd_reg

	real(dp)							:: rho_min
	real(dp)							:: rho_max
	integer								:: adapt_rho_freq
	real(dp)							:: res_thresh
	real(dp)							:: adapt_fold

	real(dp)							:: stepsize
	integer								:: adapt_step
	real(dp)							:: step_min, step_max
	integer								:: adapt_step_freq
	real(dp)							:: step_thresh
	real(dp)							:: step_fold

	real(dp)							:: nCond

	!--------------------------------------------------------
	! Preconditioner
	!--------------------------------------------------------
	logical								:: precond

	!--------------------------------------------------------
	! PO Bilinear
	!--------------------------------------------------------
	integer              				:: nbilin  != 72
	real(dp),	allocatable				:: bi_data(:,:)

	!--------------------------------------------------------
	! PO Linsys variables: KKT system
	!--------------------------------------------------------
	integer                             :: mF     
	integer                             :: nF     
	integer                             :: nnzF       
	integer	,	allocatable				:: pF(:)
	integer	,	allocatable				:: jF(:)
	real(dp),	allocatable				:: vF(:)

	integer              				:: mC
	integer              				:: nC
	integer              				:: nnzC
	integer	,	allocatable				:: pC(:)      
	integer	,	allocatable				:: jC(:)      
	real(dp),	allocatable				:: vC(:)      

	! RHS
	real(dp),	allocatable				:: rhs_e(:)
	real(dp),	allocatable				:: rhs_d(:)
	real(dp),	allocatable             :: Drho(:)	! diagonal rho

	! shift
	real(dp)							:: shift

	!--------------------------------------------------------
	! After unscaling
	!--------------------------------------------------------
	real(dp)							:: objvals(nprox)

	!--------------------------------------------------------
	! PO L1
	!--------------------------------------------------------
	real(dp)                            :: delta
	integer(ip)							:: reg_norm	! L0 or L1

	!--------------------------------------------------------
	! Proximal point, solution
	!--------------------------------------------------------
	real(dp),	allocatable		        :: n(:)		! proximal point
	integer(ip)							:: ierr

	!--------------------------------------------------------
	! Distributed Z,X,U,N,M 
	!--------------------------------------------------------
	real(dp),	allocatable             :: zpar(:)! [*]
	real(dp),	allocatable             :: xpar(:,:)! [*]
	real(dp),	allocatable             :: upar(:,:)! [*]
	real(dp),	allocatable             :: mpar(:,:)! [*]
	real(dp),	allocatable             :: npar(:,:)! [*]

	! For testing
	real(dp),	allocatable             :: zpar0(:)! [*]

	!--------------------------------------------------------
	! Adaptive rho
	!--------------------------------------------------------
    integer		   	 					:: adapt_rho

	!--------------------------------------------------------
	! Data Files
	!--------------------------------------------------------
	character(150)           			:: filebase ! = 'e_coli_core_'
	! Linsys  1
	character(150)           			:: file_F ! = 'e_coli_core_F.mtx'
	character(150)           			:: file_C ! = 'e_coli_core_C.mtx'
	character(150)           			:: file_b ! = 'e_coli_core_b.mtx'
	character(150)           			:: file_d ! = 'e_coli_core_d.mtx'
	! Bounds  1           
	character(150)           			:: file_lb!  = 'e_coli_core_lb.mtx'
	character(150)           			:: file_ub!  = 'e_coli_core_ub.mtx'
	! Bilinear1
	character(150)           			:: file_bi!  = 'e_coli_core_bi.mtx'

	! Iters file
	character(150)						:: iters_file
	integer		   	 					:: unit_iters

	!--------------------------------------------------------
	! ADMM parameters
	!--------------------------------------------------------
	real(dp)							:: abs_conv 

	real(dp)							:: stall_tol 
	! optional outputs
	real(dp)							:: prim_res, dual_res
	
	!--------------------------------------------------------
	! Writing
	!--------------------------------------------------------
	character(150)						:: format_args

	!========================================================
	! Get command-line arguments: e.g., filebase
	!========================================================

	format_args = '(A, *(E15.6))'

	! Default options file
	specs_file = 'options.spec'

    nargs = command_argument_count()

	do i = 1, nargs
		call get_command_argument(i, arg)
		select case (i)
		!----------------------------------------------------
		! Basefile or options
		!----------------------------------------------------
		case (1)
			select case(arg)
				case ('-h', '--help')
					call print_help()
					stop
				case ('-v', '--version')
					print '(2a)', 'gem_fit version ', version
					stop
				case default
					print '(a,a)', '#Fitting for model filebase: ', trim(arg)
					filebase = trim(arg)
			end select
		!----------------------------------------------------
		! Output file
		!----------------------------------------------------
		case (2)
			select case(arg)
			case default
				print '(a,a)', '#Will save to output file: ', trim(arg)
				out_file = trim(arg)
			end select
		case (3)
			select case(arg)
			case default
				print '(a,a)', '#Using specs file: ', trim(arg)
				specs_file = trim(arg)
			end select
		end select
	enddo

	!--------------------------------------------------------
	! If didn't exit by now and still don't have 2 args, 
	! need to stop since outputfile not specified
	!--------------------------------------------------------
	if (nargs < 2) then
		call print_help()
		stop
	endif

	!========================================================
	! READ options from specs file
	!========================================================
	abs_conv = 1.0d-9
	stall_tol = 1.0d-12
	delta = 1e-3
	rho = 1.0
	gam = 1.5
	max_iter = 100000
	print_freq = 1000
	shift = 0.0_dp
	rho_bi = -1
	rho_linsys= -1
	rho_bnd_reg = -1
	precond = .true.
	rho_min = 1.0d-3
	rho_max = 1.0d+3
	adapt_rho = 1
	adapt_rho_freq = 1000
	res_thresh = 1.0d-6
	adapt_fold = 2.0
	reg_norm = 1
    ! Adaptive step size for multiplier updates
	stepsize = 1.0
	adapt_step = 0
	step_fold = 0.1
	step_min = 0.01
	step_max = 1.0
	adapt_step_freq = 10000
	step_thresh = 1.0d-4
	nCond = 1.0_dp

	open(newunit = unit_specs, file=trim(specs_file), status='old')
	do
		read(unit_specs, *, iostat=ios) pname, pval
		if (ios /= 0) then
			write(*,'(A)') '#End of file or not name - value format'
			exit
		else
			write(*,'(A,A,*(E15.4))') '#Read parameter:', trim(pname),pval
			select case (trim(pname))
				case ('delta')
					delta = pval
					write(*,format_args) '#Set delta =',pval !delta
				case ('rho')
					rho = pval
					write(*,format_args) '#Set rho =',pval !rho
				case ('gamma')
					gam = pval
					write(*,format_args) '#Set gam =',pval !gam
				case ('max_iter')
					max_iter = pval
					write(*,format_args) '#Set max_iter =',pval !max_iter
				case ('print_freq')
					print_freq = pval
					write(*,format_args) '#Set print_freq =',pval !print_freq
				case ('abs_conv')
					abs_conv = pval
					write(*,format_args) '#Set abs_conv =',pval !abs_conv
				case ('stall_tol')
					stall_tol = pval
					write(*,format_args) '#Set stall_tol =',pval !stall_tol
				case ('shift')
					shift	 = pval
					write(*,format_args) '#Set shift =',pval !shift
				case ('rho_bi')
					rho_bi	 = pval
					write(*,format_args) '#Set rho_bi =',pval !rho_bi
				case ('rho_linsys')
					rho_linsys	 = pval
					write(*,format_args) '#Set rho_linsys =',pval !rho_linsys
				case ('rho_bnd_l1')
					! TODO: Deprecation warning
					rho_bnd_reg	 = pval
					write(*,format_args) '#******* rho_bnd_l1 DEPRECATED. Use rho_bnd_reg with reg_norm=1'
					write(*,format_args) '#Set rho_bnd_reg =',pval !rho_bnd_reg
				case ('rho_bnd_reg')
					rho_bnd_reg	 = pval
					write(*,format_args) '#Set rho_bnd_reg =',pval !rho_bnd_reg
				case ('precond')
					if (pval==0) precond = .false. 
					write(*,format_args) '#Set precond =',pval !precond
				case ('rho_min')
					rho_min	 = pval
					write(*,format_args) '#Set rho_min =',pval !rho_min
				case ('rho_max')
					rho_max	 = pval
					write(*,format_args) '#Set rho_max =',pval !rho_max
				case ('adapt_rho')
					adapt_rho	 = int(pval)
					write(*,format_args) '#Set adapt_rho =',pval ! adapt_rho
				case ('adapt_rho_freq')
					adapt_rho_freq	 = pval
					write(*,format_args) '#Set adapt_rho_freq =',pval ! adapt_rho_freq
				case ('res_thresh')
					res_thresh	 = pval
					write(*,format_args) '#Set res_thresh =',pval ! res_thresh
				case ('adapt_fold')
					adapt_fold	 = pval
					write(*,format_args) '#Set adapt_fold =',pval ! adapt_fold
				case ('reg_norm')
					reg_norm	 = int(pval)
					write(*,format_args) '#Set reg_norm =',pval ! reg_norm
				case ('stepsize')
					stepsize	 = pval
					write(*,format_args) '#Set stepsize =',pval ! stepsize
				case ('adapt_step')
					adapt_step	 = int(pval)
					write(*,format_args) '#Set adapt_step =',pval ! adapt_step
				case ('step_min')
					step_min	 = pval
					write(*,format_args) '#Set step_min =',pval !step_min
				case ('step_max')
					step_max	 = pval
					write(*,format_args) '#Set step_max =',pval !step_max
				case ('adapt_step_freq')
					adapt_step_freq	 = pval
					write(*,format_args) '#Set adapt_step_freq =',pval ! adapt_step_freq
				case ('step_thresh')
					step_thresh	 = pval
					write(*,format_args) '#Set step_thresh =',pval ! step_thresh
				case ('step_fold')
					step_fold	 = pval
					write(*,format_args) '#Set step_fold =',pval ! step_fold

				case ('nCond')
					nCond 		= pval
					write(*,format_args) '#Set nCond =',pval
			end select
		endif
	enddo
	close(unit_specs)

    ! Were PO-specific rhos not set by spec file?
	if (rho_bi < 0) rho_bi = rho
	if (rho_linsys < 0) rho_linsys = rho
	if (rho_bnd_reg < 0) rho_bnd_reg = rho

	ierr = 0

	!--------------------------------------------------------
	! Specify proximal operators
	!--------------------------------------------------------
	po_types(IND_LIN) = PO_IND_LINSYS_MA57        ! linsys
	select case (reg_norm)
	case (0)
		po_types(IND_BND_REG) = PO_IND_BOUNDS_CARD        ! Bounds  
	case (1)
		po_types(IND_BND_REG) = PO_IND_BOUNDS_L1        ! Bounds  
	case default
		print '(a,a,/)', 'Unknown regularization norm: ', reg_norm
	end select

	po_types(IND_BI)  = PO_IND_BILINEAR_MULT        ! Bilinear

	if (verbose) write(*,'(A)') '#Setting PO-dependent rhos:'
	po_rhos = rho
	po_rhos(IND_BI) = rho_bi	!0.1*rho
	po_rhos(IND_LIN) = rho_linsys	!10*rho
	po_rhos(IND_BND_REG) = rho_bnd_reg	!0.1*rho
	if (verbose) write(*,'(A, *(E15.6))') '#po_rhos =',po_rhos

	!--------------------------------------------------------
	! Provide problem data
	!--------------------------------------------------------

	!--------------------------------------------------------
	! Read matrices, arrays for linsys
	!--------------------------------------------------------

	write(*,'(A)') '#Reading data files...'

	file_F  = trim(filebase)//'_F.mtx'
	file_C  = trim(filebase)//'_C.mtx'
	file_b  = trim(filebase)//'_b.mtx'
	file_d  = trim(filebase)//'_d.mtx'
	file_lb = trim(filebase)//'_lb.mtx'
	file_ub = trim(filebase)//'_ub.mtx'
	file_bi = trim(filebase)//'_bi.mtx'

	write(*,'(A)') '#Reading F'
	call csr_from_mm(file_F, pF, jF, vF, nnzF, mF, nF)
	write(*,'(A)') '#Reading C'
	call csr_from_mm(file_C, pC, jC, vC, nnzC, mC, nC)

	write(*,'(A)') '#Reading b'
	call array_from_mm(file_b, rhs_e) 
	write(*,'(A)') '#Reading d'
	call array_from_mm(file_d, rhs_d) 

	write(*,'(A)') '#Reading lb'
	call array_from_mm(file_lb, lb) 
	write(*,'(A)') '#Reading ub'
	call array_from_mm(file_ub, ub) 

	write(*,'(A)') '#Reading bi'
	call matrix_from_mm(file_bi, bi_data) 
    !************************************************************
	! TODO: DEBUG
    !************************************************************
	write(*,'(A,*(I15))') '#bi_data shape:',shape(bi_data)
	write(*,format_args) '#bi_data:',bi_data

	!--------------------------------------------------------
	! Get sizes and allocate arrays
	!--------------------------------------------------------
	write(*,'(A)') '#Allocating arrays'

	lenx = size(lb,1)
	if (lenx /= size(ub,1)) then
!		error stop 'size of lb and ub must match!'
		stop 'size of lb and ub must match!'
	endif

	write(*,'(A)') '#Allocating n'
	allocate(n(lenx))

	write(*,'(A)') '#Allocating zpar,xpar,upar,mpar,npar'
	allocate(zpar(lenx))
	allocate(xpar(lenx,nprox))
	allocate(upar(lenx,nprox))
	allocate(mpar(lenx,nprox))
	allocate(npar(lenx,nprox))

	!--------------------------------------------------------
	! Initialize variable data
	!--------------------------------------------------------
	call random_number(zpar)
	call random_number(xpar)

	upar = 0.0_dp
	mpar = 0.0_dp
	npar = 0.0_dp

	! Initialize proximal point, solution 
	n = 0.0_dp

	!--------------------------------------------------------
	! Parse bilinear data
	!--------------------------------------------------------
	allocate(ds(size(bi_data,1)))

    ds  = bi_data(:,1)
    nxs = bi_data(:,2)

    xstarts = int(bi_data(:,3))
    xends   = int(bi_data(:,4))
    ystarts = int(bi_data(:,5))
    yends   = int(bi_data(:,6))

	!********************************************************
	! TODO: DEBUG
    ! write(*,*) '#ds:',ds
    ! write(*,*) '#nxs:',nxs
    ! write(*,*) '#xstarts:',xstarts
    ! write(*,*) '#xends:',xends
    ! write(*,*) '#ystarts:',ystarts
    ! write(*,*) '#yends:',yends
    ! write(*,*) '#shift:',shift
	!********************************************************

	!--------------------------------------------------------
	! Linsys
	!--------------------------------------------------------
	allocate(Drho(nF))
	Drho = rho_linsys

	! Construct problem:
	! TODO: Issue: if shift != 0, then nnzF is potentially 
	! wrong--should change to reflect all nonzero diag
	call make_PO_matrices(mF, nF, nnzF, pF, jF, vF,	&
						  mC, nC, nnzC, pC, jC, vC,	&
						  rhs_e, rhs_d, Drho, npar(:,IND_LIN), &
						  shift, nCond)

	write(*,'(A)') "#Finished make_PO_matrices"
	!--------------------------------------------------------
	! Bilinear
	!--------------------------------------------------------
	! x'y >= d
	csenses = [(.true., i=1,size(ds))]
	!--------------------------------------------------------

	!--------------------------------------------------------
	! L1
	select case (reg_norm)
	case (0) 
        l0max = int(delta)
	case (1)
		allocate(l1weights(nF))
		l1weights = 0.0_dp
		! L1 only on y
		do i=1, size(ds)
			l1weights(ystarts(i):yends(i)) = delta
		enddo
	case default
		print '(a,a)', '#Unknown regularization norm: ', reg_norm
	end select

	!l1weight = delta / dble(lenx)
	!write(*,*) "Normalizing L1 weight by problem dimension to", l1weight
	!--------------------------------------------------------

	!--------------------------------------------------------
	! Precondition all PO data
	!--------------------------------------------------------
	if (precond) then
		if (verbose) write(*,'(A)') '#Preconditioning all POs'

		do i=1, size(po_types)
			call scale_po(po_types(i), [(po_rhos(i),j=1,lenx)])
		enddo
	endif

	!--------------------------------------------------------
	! QC the scaled data
	!--------------------------------------------------------
	if (verbose) then
		if (precond) then
			write(*,format_args) '#min(SD) =', minval(SD)
			write(*,format_args) '#max(SD) =', maxval(SD)
			write(*,format_args) '#min(SE) =', minval(SE)
			write(*,format_args) '#max(SE) =', maxval(SE)
		endif

		write(*,format_args) '#max(abs(lb)) =', maxval(abs(lb))
		write(*,format_args) '#max(abs(ub)) =', maxval(abs(ub))

		write(*,format_args) '#min(lb) =', minval(lb)
		write(*,format_args) '#max(lb) =', maxval(lb)
		write(*,format_args) '#min(ub) =', minval(ub)
		write(*,format_args) '#max(ub) =', maxval(ub)

		write(*,format_args) '#Bilinear'
		write(*,format_args) '#min(ds) =', minval(ds)
		write(*,format_args) '#max(ds) =', maxval(ds)
		write(*,format_args) '#min(ds) =', minval(ds)
		write(*,format_args) '#max(ds) =', maxval(ds)

		write(*,format_args) '#Linsys'
		write(*,format_args) '#min(b1) =', minval(b1)
		write(*,format_args) '#max(b1) =', maxval(b1)
		write(*,format_args) '#min(b2) =', minval(b2)
		write(*,format_args) '#max(b2) =', maxval(b2)
		!----------------------------------------------------
		! Confirm that shifted zero block untouched
		! by preconditioner
		!----------------------------------------------------
		write(*,format_args) '#min(vA22) =', minval(vA22)
		write(*,format_args) '#max(vA22) =', maxval(vA22)


		if (reg_norm .eq. 1) then
			write(*,'(A)') '#L1'
			write(*,format_args) '#min(l1weights) =', minval(l1weights)
			write(*,format_args) '#max(l1weights) =', maxval(l1weights)
		endif
	endif

	! Check for obvious infeasibilities
	if ( any(lb > ub) ) then
!		error stop "lb > ub"
		stop "lb > ub"
	endif

	!--------------------------------------------------------
	! ADMM part
	!--------------------------------------------------------

	if (verbose) write(*,'(A)') '#Running ADMM'

	! Write iters to file
	iters_file = trim(filebase)//'.iters'
	open(newunit = unit_iters, file=trim(iters_file), status='replace')

	call cpu_time(tic)

    call admms( inform, objval, iter,	&
		zpar, lenx, &
		xpar, nprox, & 
		rho, gam, max_iter, &
		po_types, &
		verbose=verbose, &
		print_freq_in = print_freq, abs_conv_in=abs_conv,		&
		po_rhos = po_rhos, Uopt=upar, Nopt=npar, SD=SD,			&
		adapt_rho = adapt_rho, adapt_rho_freq=adapt_rho_freq, 	&
		rho_min=rho_min, rho_max=rho_max, res_thresh=res_thresh,&
		adapt_fold=adapt_fold, prim_res_out=prim_res, dual_res_out=dual_res,&
		stepsize = stepsize,									&
		adapt_step=adapt_step, adapt_step_freq=adapt_step_freq,	&
		step_min=step_min, step_max=step_max, step_fold=step_fold,	&
		step_thresh=step_thresh, unit_iters=unit_iters,			&
		stall_tol=stall_tol, nCond=nCond) !iters_file=iters_file)

	close(unit_iters)

	!--------------------------------------------------------
	! Unscale solution back
	!--------------------------------------------------------
	if (precond) then
		call unscale_ruiz(zpar, SD)
		do j=1, nprox
			call unscale_ruiz(xpar(:,j), SD)
			call unscale_ruiz(upar(:,j), SD)
			call unscale_ruiz(npar(:,j), SD)

			! Also unscale bounds, l1weights, bilinear d, etc.
			call unscale_po(po_types(j))
		enddo
	endif

	!--------------------------------------------------------
	! Time
	!--------------------------------------------------------
	call cpu_time(toc)
	t_elapsed = toc - tic

	!--------------------------------------------------------
	! Recompute objective using unscaled sol
	!--------------------------------------------------------
	if (precond) then
		do i=1, nprox
			objvals(i) = calc_objval(po_types(i), zpar)
		enddo
		objval = sum(objvals)
	endif


	if (verbose) then
		write(*,'(A)') '#--------------------------------------'
		write(*,'(A,I15)') '#Finished with inform =',inform 
		write(*,format_args) '#Objval =',objval

		write(*,format_args) '#Z(1:5) =', zpar(1:5)
		write(*,format_args) '#U(1:5,:) =', upar(1:5,:)

		write(*,format_args) '#X(1:5,:) =', xpar(1:5,:)
		write(*,format_args) '#N(1:5,:) =', npar(1:5,:)
	endif

	!--------------------------------------------------------
	! Save solution
	!--------------------------------------------------------
	call write_stats()
	call write_solution()

	!--------------------------------------------------------
	! Final solution stats
	!--------------------------------------------------------
	if (verbose) then
		call show_solution_quality(zpar, ds, xstarts, xends,	&
			ystarts, yends, mC, vC, jC, pC, rhs_d, lb, ub)
	endif

	!--------------------------------------------------------
	! Some diagnostics
	!--------------------------------------------------------
	if (verbose) then
		write(*,'(A)') '#--------------------------------------'
		write(*,format_args) '#Primal res: max |X-Z| =', maxval(abs(xpar-spread(zpar,2,nprox)))
		write(*,format_args) '#For linsys: max |X-Z| =', maxval(abs(xpar(:,IND_LIN)-zpar))
		write(*,format_args) '#For bnds L1: max |X-Z| =', maxval(abs(xpar(:,IND_BND_REG)-zpar))
		write(*,format_args) '#For bilinear: max |X-Z| =', maxval(abs(xpar(:,IND_BI)-zpar))
		write(*,format_args) '#max |Zt| =', maxval(abs(zpar))
		write(*,format_args) '#max |Ut| =', maxval(abs(upar))
		write(*,format_args) '#max |Xt| =', maxval(abs(xpar))
		write(*,format_args) '#max |Nt| =', maxval(abs(npar))

		write(*,format_args) '#max |Nt+1|', maxval(abs(spread(zpar,2,nprox) - upar))
		allocate(zpar0(size(zpar)))
		zpar0 = zpar
		write(*,'(A,*(I15))') '#shape(xpar+upar) =', shape(xpar+upar)
		write(*,'(A,*(I15))') '#shape(spread(po_rhos,1,lenx)) =',shape(spread(po_rhos,1,lenx))
		zpar = sum(spread(po_rhos,1,lenx)*(xpar+upar),2)
		write(*,format_args) '#max |Zt+1|', maxval(abs(zpar))
		call update_multipliers(lenx, upar, xpar, zpar, zpar0, gam)
		write(*,format_args) '#max |Ut+1|', maxval(abs(upar))
		deallocate(zpar0)
	endif

	!--------------------------------------------------------
	! 20 Sep 2018 [LY]: be nice to compiler and deallocate
	! everything to hush valgrind
	!--------------------------------------------------------
	deallocate(bi_data)
	deallocate(pF)
	deallocate(jF)
	deallocate(vF)
	deallocate(pC)
	deallocate(jC)
	deallocate(vC)
	deallocate(rhs_e)
	deallocate(rhs_d)
	deallocate(Drho)
	deallocate(n)
	deallocate(zpar)
	deallocate(xpar)
	deallocate(upar)
	deallocate(mpar)
	deallocate(npar)
	if (allocated(l1weights)) deallocate(l1weights)


contains

    subroutine print_help()
		print '(a)', 'usage: gem_fit [filebase] [outfile] [specsfile]'
		print '(a)', ''
		print '(a)', 'filebase: base for the following files (must be present):'
		print '(a)', '	[filebase]_F.mtx'
		print '(a)', '	[filebase]_C.mtx'
		print '(a)', '	[filebase]_b.mtx'
		print '(a)', '	[filebase]_d.mtx'
		print '(a)', '	[filebase]_lb.mtx'
		print '(a)', '	[filebase]_ub.mtx'
		print '(a)', '	[filebase]_bi.mtx'
		print '(a)', ''
		print '(a)', 'outfile: solution saved to this file. '
		print '(a)', ''
		print '(a)', '(optional) specsfile: additional specs file.'
		print '(a)', ''
		print '(a)', 'cmdline options:'
		print '(a)', '	-h, --help        print usage information and exit'
		print '(a)', '	-v, --version     print version information and exit'
	end subroutine

	subroutine show_solution_quality(zopt, bi_d,			&
									 xinds0, xinds1,		&
									 yinds0, yinds1,		&
									 mM, vM, jM, pM, rhs,	&
									 xl, xu)

		!----------------------------------------------------
		! Show solution quality: constraint violations
		!----------------------------------------------------
		real(dp),		intent(in)			:: zopt(:)
		real(dp),		intent(in)			:: bi_d(:)						
		integer(ip),	intent(in)			:: xinds0(:)
		integer(ip),	intent(in)			:: xinds1(:)
		integer(ip),	intent(in)			:: yinds0(:)
		integer(ip),	intent(in)			:: yinds1(:)

		real(dp)							:: xy, cons_viol
		integer(ip),	intent(in)			:: mM
		real(dp),		intent(in)			:: vM(:)
		integer(ip),	intent(in)			:: jM(:)
		integer(ip),	intent(in)			:: pM(:)
		real(dp),		intent(in)			:: rhs(:)
		real(dp),		intent(in)			:: xl(:)
		real(dp),		intent(in)			:: xu(:)

        real(dp)							:: Mx(mM)

		integer(ip)							:: nprint
		integer(ip)							:: i,j,k

		!--------------------------------------------------------
		! Final constraint viol
		!--------------------------------------------------------
		do i=1, size(bi_d)
			xy = dot_product( zopt(xinds0(i):xinds1(i)),	&
							  zopt(yinds0(i):yinds1(i)))  
			cons_viol = bi_d(i) - xy

			write(*,'(A,I15)') '#Cond:',i
			write(*,'(A)') '#--------------------------------------'
			write(*,format_args) '#Final bilinear x*y =', xy
			write(*,format_args) '#Final bilinear d =', bi_d(i)
			write(*,format_args) '#final bilinear cons viol =', cons_viol
		enddo
		
		write(*,'(A)') '#--------------------------------------'
		write(*,format_args) '#Final max bound viol, lb:', maxval(xl - zopt)
		write(*,format_args) '#Final max bound viol, ub:', maxval(zopt - xu)
		
		nprint = 0
		if (maxval(xl - zopt) > abs_conv) then
			if (rho_bnd_reg > 0.0_dp) then
				write(*,'(A)') '#--------------------------------------'
				write(*,'(A)') '#Infeas lbs:'
				write(*,'(A15, A15, A15)') '#ind', 'xl','zopt'
				do i=1, lenx
					if ( (xl(i)-zopt(i)) > abs_conv) then
						write(*,'(A, I15, E15.4, E15.4)') '#', i, xl(i), zopt(i)
						nprint = nprint + 1
						if (nprint >= 10) then
							exit
						endif
					endif
				enddo
				write(*,'(A)') '#Only showed first 10 infeasible'
			endif
		endif

		Mx = 0.0
		call amux(mM, zopt, Mx, vM, jM, pM)
		write(*,'(A)') '#--------------------------------------'
		write(*,format_args) '#Final ||Cx - d||inf =', maxval(abs(Mx-rhs))
		write(*,format_args) '#Final ||Cx - d||2 =', sum((Mx-rhs)**2)

	end subroutine

	subroutine write_solution()

		!----------------------------------------------------
		! Write soution to mtx
		!----------------------------------------------------

		if (verbose) write(*,'(*(A))') '#Writing solution Z to file:', trim(out_file)

		id = '%%MatrixMarket'
		type = 'matrix'
		rep = 'array'
		field = 'double'
		symm = 'general'

		open(newunit = unit_out, file=trim(out_file), status='replace')

		call mm_file_write ( unit_out, id, type, rep, field, symm,				&
							size(zpar,1), 1, product(shape(zpar)), 	&
							indx, jndx, ival, rval, zpar, cval)

		close(unit_out)

		if (verbose) write(*,'(A)') "#Done"

	end subroutine

	subroutine write_stats()

		!----------------------------------------------------
		! Write stats
		!----------------------------------------------------

		character(150)           			:: stats_file
		integer(ip)							:: unit_stats
		character(150)						:: format_stats

		format_stats = "(*(a,a,e15.8))"

		stats_file = trim(filebase)//'.stats'

		if (verbose) write(*,'(*(A))') '#Writing stats to file:', trim(stats_file)

		open(newunit = unit_stats, file=trim(stats_file), status='replace')

		!----------------------------------------------------
		! inform
		! obval
		! iter
		! time
		!----------------------------------------------------
		write(unit_stats, format_stats) 'inform', ',', dble(inform)
		write(unit_stats, format_stats) 'objval', ',', objval
		write(unit_stats, format_stats) 'iter', ',', dble(iter)
		write(unit_stats, format_stats) 'time', ',', t_elapsed
		write(unit_stats, format_stats) 'prim_res', ',', prim_res
		write(unit_stats, format_stats) 'dual_res', ',', dual_res

		close(unit_stats)

		if (verbose) write(*,'(A)') '#Done'

	end subroutine

end program
