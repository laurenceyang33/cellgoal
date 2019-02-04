!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File fgfitsg.f90: 
!
! Flux fit problem for any model using fine-grained ADMM.
!
! Laurence Yang, UCSD
!
! 05 Nov 2018: first version. 
! 20 Nov 2018: each PO group (PORG) optimizes its fine-graining.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program fgfitsg

	use admm_fg_groups,	only: admms
	use po_constants

	use mm_helper,	only: csr_from_mm, array_from_mm, matrix_from_mm
	use mm_io,		only: mm_file_write
	use sparsekit
	use sparseops
	use pog_linsys,	only: pA,jA,vA,nnzA,mA,nA, rhs_b,	&
						  pF,jF,vF,nnzF,mF,nF, rhs_e

	use po_bounds,		only: lb, ub
	use po_L1,			only: l1weights
	use po_bilinear_multi,	only: ds, csenses, nxs, xstarts, xends, ystarts, yends

	use pogs_factory,	only: create_po_graph
	use precondition,	only: SD, SE, unscale_ruiz

	implicit none

	!--------------------------------------------------------
	integer,	parameter				:: ip = 4, dp = 8, sp=4
	integer								:: i,j,max_iter,print_freq
	logical            					:: verbose=.false.
	character(*), parameter				:: version = '0.1'

	!--------------------------------------------------------
	! Command arguments
	!--------------------------------------------------------
	character(20)						:: specs_file 	! = 'options.spec'
	integer(ip)							:: unit_specs
	integer(ip)							:: ios
	character(50)						:: pname
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
	integer,	parameter				:: npogs = 3
	integer								:: nprox

	integer,	parameter				:: IND_LIN  = 1
	integer,	parameter				:: IND_BND_L1  = 2
	integer,	parameter				:: IND_BI 	= 3

	real(dp)							:: rho
	real(dp)							:: gam
	integer,	allocatable				:: po_types(:)
	real(dp),	allocatable				:: po_rhos(:)
	integer,	allocatable				:: po_sizes(:)
	integer,	allocatable				:: po_starts(:)

	real(dp)							:: rho_bi
	real(dp)							:: rho_linsys
	real(dp)							:: rho_bnd_l1

	real(dp)							:: rho_min
	real(dp)							:: rho_max
	integer								:: adapt_rho_freq
	real(dp)							:: res_thresh
	real(dp)							:: adapt_fold

	real(dp)							:: stepsize

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
	real(dp),	allocatable             :: Drho(:)	! diagonal rho

	! shift
	real(dp)							:: shift

	!--------------------------------------------------------
	! Graph filter in CSR
	!--------------------------------------------------------
	integer,	allocatable				:: pgraph(:)
	integer,	allocatable				:: jgraph(:)

	!--------------------------------------------------------
	! After unscaling
	!--------------------------------------------------------
	real(dp),	allocatable				:: objvals(:)

	!--------------------------------------------------------
	! PO L1
	!--------------------------------------------------------
	real(dp)                            :: delta

	!--------------------------------------------------------
	! Proximal point, solution
	!--------------------------------------------------------
	real(dp),	allocatable		        :: n(:)		! proximal point
	integer(ip)							:: ierr

	!--------------------------------------------------------
	! Distributed Z,X,U,N,M 
	!--------------------------------------------------------
	real(dp),	allocatable             :: zpar(:)
	real(dp),	allocatable             :: xpar(:,:)
	real(dp),	allocatable             :: upar(:,:)
	real(dp),	allocatable             :: mpar(:,:)
	real(dp),	allocatable             :: npar(:,:)

	! For testing
	real(dp),	allocatable             :: zpar0(:)

	!--------------------------------------------------------
	! Adaptive rho
	!--------------------------------------------------------
    integer		   	 					:: adapt_rho

	!--------------------------------------------------------
	! Data Files
	!--------------------------------------------------------
	character(150)           			:: filebase
	! Linsys  1
	character(150)           			:: file_F 
	character(150)           			:: file_C 
	character(150)           			:: file_b 
	character(150)           			:: file_d 
	! Bounds  1           
	character(150)           			:: file_lb
	character(150)           			:: file_ub
	! Bilinear1
	character(150)           			:: file_bi

	! Iters file
	character(150)						:: iters_file
	integer		   	 					:: unit_iters
	character(150)						:: iters_file_num

	integer								:: istart

	!--------------------------------------------------------
	! ADMM parameters
	!--------------------------------------------------------
	real(dp)							:: abs_conv 

	! optional outputs
	real(dp)							:: prim_res, dual_res
	

	!========================================================
	! Get command-line arguments: e.g., filebase
	!========================================================

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
					print '(2a)', '# gem_fit version ', version
					stop
				case default
					print '(a,a,/)', '# Fitting for model filebase: ', arg
					filebase = trim(arg)
			end select
		!----------------------------------------------------
		! Output file
		!----------------------------------------------------
		case (2)
			select case(arg)
			case default
				print '(a,a,/)', '# Will save to output file: ', arg
				out_file = trim(arg)
			end select
		case (3)
			select case(arg)
			case default
				print '(a,a,/)', '# Using specs file: ', arg
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
	delta = 1e-3
	rho = 1.0
	gam = 1.5
	max_iter = 100000
	print_freq = 1000
	shift = 0.0_dp
	rho_bi = -1
	rho_linsys= -1
	rho_bnd_l1 = -1

	precond = .false. !.true.

	rho_min = 1.0d-3
	rho_max = 1.0d+3
	adapt_rho = 1
	adapt_rho_freq = 1000
	res_thresh = 1.0d-6
	adapt_fold = 2.0
	stepsize = 1.0

	iters_file = trim(filebase)//'.iters'

	open(newunit = unit_specs, file=trim(specs_file), status='old')
	do
		read(unit_specs, *, iostat=ios) pname, pval
		if (ios /= 0) then
			if (verbose) write(*,*) '# End of file or not name - value format'
			exit
		else
			if (verbose) write(*,*) '# Read parameter:', pname,pval
			select case (trim(pname))
				case ('delta')
					delta = pval
					if (verbose) write(*,*) '# Set delta =',delta
				case ('rho')
					rho = pval
					if (verbose) write(*,*) '# Set rho =',rho
				case ('gamma')
					gam = pval
					if (verbose) write(*,*) '# Set gam =',gam
				case ('max_iter')
					max_iter = pval
					if (verbose) write(*,*) '# Set max_iter =',max_iter
				case ('print_freq')
					print_freq = pval
					if (verbose) write(*,*) '# Set print_freq =',print_freq
				case ('abs_conv')
					abs_conv = pval
					if (verbose) write(*,*) '# Set abs_conv =',abs_conv
				case ('shift')
					shift	 = pval
					if (verbose) write(*,*) '# Set shift =',shift
				case ('rho_bi')
					rho_bi	 = pval
					if (verbose) write(*,*) '# Set rho_bi =',rho_bi
				case ('rho_linsys')
					rho_linsys	 = pval
					if (verbose) write(*,*) '# Set rho_linsys =',rho_linsys
				case ('rho_bnd_l1')
					rho_bnd_l1	 = pval
					if (verbose) write(*,*) '# Set rho_bnd_l1 =',rho_bnd_l1
				case ('precond')
					if (pval==0) precond = .false. 
					if (verbose) write(*,*) '# Set precond =',precond
				case ('rho_min')
					rho_min	 = pval
					if (verbose) write(*,*) '# Set rho_min =',rho_min
				case ('rho_max')
					rho_max	 = pval
					if (verbose) write(*,*) '# Set rho_max =',rho_max
				case ('adapt_rho')
					adapt_rho	 = int(pval)
					if (verbose) write(*,*) '# Set adapt_rho =', adapt_rho
				case ('adapt_rho_freq')
					adapt_rho_freq	 = pval
					if (verbose) write(*,*) '# Set adapt_rho_freq =', adapt_rho_freq
				case ('res_thresh')
					res_thresh	 = pval
					if (verbose) write(*,*) '# Set res_thresh =', res_thresh
				case ('adapt_fold')
					adapt_fold	 = pval
					if (verbose) write(*,*) '# Set adapt_fold =', adapt_fold
				case ('stepsize')
					stepsize	 = pval
					if (verbose) write(*,*) '# Set stepsize =', stepsize
				case ('iters_file_num')
					write(iters_file_num, '(I5.5)') int(pval)
					iters_file = trim(filebase)//trim(iters_file_num)//'.iters'
					if (verbose) write(*,*) '# iters_file_num =', iters_file
				case ('verbose')
					if (pval==0) then
						verbose = .false.
					else
						verbose = .true.
					endif
					if (verbose) write(*,*) '# Set verbose =', verbose
			end select
		endif
	enddo
	close(unit_specs)

    ! Were PO-specific rhos not set by spec file?
	if (rho_bi < 0) rho_bi = rho
	if (rho_linsys < 0) rho_linsys = rho
	if (rho_bnd_l1 < 0) rho_bnd_l1 = rho

	ierr = 0

	!--------------------------------------------------------
	! Read matrices, arrays for linsys
	!--------------------------------------------------------

	if (verbose) write(*,*) '# Reading data files...'

	file_F  = trim(filebase)//'_F.mtx'
	file_C  = trim(filebase)//'_C.mtx'
	file_b  = trim(filebase)//'_b.mtx'
	file_d  = trim(filebase)//'_d.mtx'
	file_lb = trim(filebase)//'_lb.mtx'
	file_ub = trim(filebase)//'_ub.mtx'
	file_bi = trim(filebase)//'_bi.mtx'

	if (verbose) write(*,*) '# Reading F'
	call csr_from_mm(file_F, pF, jF, vF, nnzF, mF, nF)
	if (verbose) write(*,*) '# Reading C'
	call csr_from_mm(file_C, pA, jA, vA, nnzA, mA, nA)

	if (verbose) write(*,*) '# Reading b'
	call array_from_mm(file_b, rhs_e) 
	if (verbose) write(*,*) '# Reading d'
	call array_from_mm(file_d, rhs_b) 

	if (verbose) write(*,*) '# Reading lb'
	call array_from_mm(file_lb, lb) 
	if (verbose) write(*,*) '# Reading ub'
	call array_from_mm(file_ub, ub) 

	if (verbose) write(*,*) '# Reading bi'
	call matrix_from_mm(file_bi, bi_data) 
    !************************************************************
	! TODO: DEBUG
    !************************************************************
	if (verbose) write(*,*) '# bi_data shape:',shape(bi_data)
	if (verbose) write(*,*) '# bi_data:',bi_data

	!--------------------------------------------------------
	! Specify proximal operators
	!--------------------------------------------------------
	allocate(po_types(npogs))
	allocate(po_rhos(npogs))
    allocate(po_sizes(npogs))
    allocate(po_starts(npogs))

	po_types(IND_LIN) = POG_IND_LINSYS
	po_types(IND_BND_L1) = POG_IND_BOUNDS_L1        ! Bounds + L1
	po_types(IND_BI)  = POG_IND_BILINEAR_MULT_CPU        ! Bilinear

	if (verbose) write(*,*) '# size of pA:', size(pA)
	if (verbose) write(*,*) '# size of pF:', size(pF)
	if (verbose) write(*,*) '# mF:', mF
	if (verbose) write(*,*) '# nF:', nF

	po_sizes(IND_LIN) = mA
	po_sizes(IND_BND_L1) = 1
	po_sizes(IND_BI) = size(ds)  ! nConds

	istart = 0
	do i=1, npogs
		po_starts(i) = istart
		istart = istart + po_sizes(i)
	enddo

	nprox = sum(po_sizes)

	if (verbose) write(*,*) '# Allocated po_types:',size(po_types)

	if (verbose) write(*,*) '# Setting PO-dependent rhos:'
	po_rhos = rho
	!--------------------------------------------------------
	! Get sizes and allocate arrays
	!--------------------------------------------------------
	if (verbose) write(*,*) '# Allocating arrays'

	lenx = size(lb,1)
	if (lenx /= size(ub,1)) then
		stop 'size of lb and ub must match!'
	endif

	if (verbose) write(*,*) '# Allocating n'
	allocate(n(lenx))

	if (verbose) write(*,*) '# Allocating zpar,xpar,upar,mpar,npar'
	allocate(zpar(lenx))
	!--------------------------------------------------------
	! Too many proximal operators to be stored as dense arrays
	!--------------------------------------------------------
	allocate(xpar(lenx,nprox))
	allocate(upar(lenx,nprox))
	allocate(mpar(lenx,nprox))
	allocate(npar(lenx,nprox))

	if (verbose) write(*,*) '# lenx:',lenx,'nprox:',nprox
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

	if (verbose) write(*,*) 'Done initializing upar,mpar,npar,n'
	!--------------------------------------------------------
	! Parse bilinear data
	!--------------------------------------------------------
    ds  = bi_data(:,1)
    nxs = bi_data(:,2)

    xstarts = int(bi_data(:,3))
    xends   = int(bi_data(:,4))
    ystarts = int(bi_data(:,5))
    yends   = int(bi_data(:,6))

	if (verbose) then
		write(*,*) '# ds:',ds
		write(*,*) '# nxs:',nxs
		write(*,*) '# xstarts:',xstarts
		write(*,*) '# xends:',xends
		write(*,*) '# ystarts:',ystarts
		write(*,*) '# yends:',yends
		write(*,*) '# shift:',shift
	endif
	!********************************************************

	!--------------------------------------------------------
	! Linsys
	!--------------------------------------------------------
	!--------------------------------------------------------
	! Bilinear
	!--------------------------------------------------------
	! x'y >= d
	csenses = [(.true., i=1,size(ds))]
	!--------------------------------------------------------

	!--------------------------------------------------------
	! L1
	allocate(l1weights(nF))
	l1weights = 0.0_dp
	! L1 only on y
	do i=1, size(ds)
		l1weights(ystarts(i):yends(i)) = delta
	enddo

	! Check for obvious infeasibilities
	if ( any(lb > ub) ) then
		stop "lb > ub"
	endif
	!--------------------------------------------------------
	! Create the PO graph
	!--------------------------------------------------------
	allocate(pgraph(lenx+1))
	if (verbose) write(*,*) 'create_po_graph...'
    call create_po_graph(po_types, po_starts, po_sizes, lenx, npogs, nprox, pgraph, jgraph)
	if (verbose) write(*,*) '... done.'

	if (verbose) write(*,*) '# lenx:', lenx
	if (verbose) write(*,*) '# size of po_types:', size(po_types)
	if (verbose) write(*,*) '# size of pgraph:', size(pgraph)
	if (verbose) write(*,*) '# size of jgraph:', size(jgraph)
	if (verbose) write(*,*) '# max(jgraph):', maxval(jgraph)

	!--------------------------------------------------------
	! ADMM part
	!--------------------------------------------------------
	if (verbose) write(*,*) '# po_types:', po_types

	! Write iters to file
	open(newunit = unit_iters, file=trim(iters_file), status='replace')

	if (verbose) write(*,*) '# Running ADMM'

	call cpu_time(tic)

    call admms( inform, objval, iter,		&
		zpar, lenx, xpar, nprox, npogs, 	& 
		rho, gam, max_iter, 				&
		po_types, po_sizes, po_starts,		& 
		pgraph, jgraph,						&
		verbose=verbose, 					&
		print_freq_in = print_freq, abs_conv_in=abs_conv,		&
		po_rhos = po_rhos, Uopt=upar, Nopt=npar, SD=SD,			&
		adapt_rho = adapt_rho, adapt_rho_freq=adapt_rho_freq, 	&
		rho_min=rho_min, rho_max=rho_max, res_thresh=res_thresh,&
		adapt_fold=adapt_fold, prim_res_out=prim_res, dual_res_out=dual_res,&
		stepsize=stepsize, unit_iters=unit_iters)

	!--------------------------------------------------------
	! Time
	!--------------------------------------------------------
	call cpu_time(toc)
	t_elapsed = toc - tic

	!--------------------------------------------------------
	! Recompute objective using unscaled sol
	!--------------------------------------------------------
	allocate(objvals(nprox))
	objvals = 0

	if (verbose) then
		write(*,*) '# --------------------------------------'
		write(*,*) '# Finished with inform =',inform 
		write(*,*) '# Objval =',objval

		write(*,*) '# Z(1:5) =', zpar(1:5)
		write(*,*) '# U(1:5,1:5) =', upar(1:5,1:5)

		write(*,*) '# X(1:5,1:5) =', xpar(1:5,1:5)
		write(*,*) '# N(1:5,1:5) =', npar(1:5,1:5)
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
			ystarts, yends, mA, vA, jA, pA, rhs_b, lb, ub)
	endif

	!--------------------------------------------------------
	! Some diagnostics
	!--------------------------------------------------------
	write(*,*) '# --------------------------------------'
	write(*,*) '# Primal res: max |X-Z| =', maxval(abs(xpar-spread(zpar,2,nprox)))
	write(*,*) '# max |Zt| =', maxval(abs(zpar))
	write(*,*) '# max |Ut| =', maxval(abs(upar))
	write(*,*) '# max |Xt| =', maxval(abs(xpar))
	write(*,*) '# max |Nt| =', maxval(abs(npar))

	write(*,*) '# max |Nt+1|', maxval(abs(spread(zpar,2,nprox) - upar))
	allocate(zpar0(size(zpar)))
	zpar0 = zpar
	write(*,*) '# shape(xpar+upar) =', shape(xpar+upar)
	write(*,*) '# shape(spread(po_rhos,1,lenx)) =',shape(spread(po_rhos,1,lenx))
	zpar = sum(spread(po_rhos,1,lenx)*(xpar+upar),2)
	write(*,*) '# max |Zt+1|', maxval(abs(zpar))

	!--------------------------------------------------------
	! Deallocate
	!--------------------------------------------------------
	if (allocated(bi_data)) deallocate(bi_data)
	deallocate(pF)
	deallocate(jF)
	deallocate(vF)
	deallocate(pA)
	deallocate(jA)
	deallocate(vA)
	deallocate(rhs_e)
	deallocate(rhs_b)
	if (allocated(Drho)) deallocate(Drho)
	deallocate(n)
	!
	deallocate(zpar)
	deallocate(xpar)
	deallocate(upar)
	deallocate(mpar)
	deallocate(npar)
	!
	deallocate(pgraph)
	deallocate(jgraph)
	if (allocated(l1weights)) deallocate(l1weights)

	deallocate(po_types)
	deallocate(po_rhos)
	deallocate(po_starts)
	deallocate(po_sizes)
	if (allocated(objvals)) deallocate(objvals)
	deallocate(zpar0)

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

			if (verbose) then
				write(*,*) '# Cond:',i
				write(*,*) '# --------------------------------------'
				write(*,*) '# Final bilinear x*y =', xy
				write(*,*) '# Final bilinear d =', bi_d(i)
				write(*,*) '# final bilinear cons viol =', cons_viol
			endif
		enddo
		
		if (verbose) then
			write(*,*) '# --------------------------------------'
			write(*,*) '# Final max bound viol, lb:', maxval(xl - zopt)
			write(*,*) '# Final max bound viol, ub:', maxval(zopt - xu)
		endif
		
		nprint = 0
		if (maxval(xl - zopt) > abs_conv) then
			if (rho_bnd_l1 > 0.0_dp) then
				if (verbose) then
					write(*,*) '# --------------------------------------'
					write(*,*) '# Infeas lbs:'
					write(*,'(A15, A15, A15)') 'ind', 'xl','zopt'
				endif
				do i=1, lenx
					if ( (xl(i)-zopt(i)) > abs_conv) then
						if (verbose) write(*,'(I15, E15.4, E15.4)') i, xl(i), zopt(i)
						nprint = nprint + 1
						if (nprint >= 10) then
							exit
						endif
					endif
				enddo
				if (verbose) write(*,*) '# Only showed first 10 infeasible'
			endif
		endif

		Mx = 0.0
		call amux(mM, zopt, Mx, vM, jM, pM)
        if (verbose) then
			write(*,*) '# --------------------------------------'
			write(*,*) '# Final ||Ax - b||inf =', maxval(abs(Mx-rhs))
			write(*,*) '# Final ||Ax - b||2 =', sum((Mx-rhs)**2)
		endif

	end subroutine

	subroutine write_solution()

		!----------------------------------------------------
		! Write soution to mtx
		!----------------------------------------------------

		if (verbose) write(*,*) '# Writing solution Z to file:', out_file

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

		if (verbose) write(*,*) "Done"

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

		if (verbose) write(*,*) '# Writing stats to file:', stats_file

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

		if (verbose) write(*,*) "Done"

	end subroutine

end program
