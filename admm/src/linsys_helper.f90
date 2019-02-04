module linsys_helper

	!----------------------------------------------------
	! Helper for creating matrices of linear constraints PO:
	!
	! min 1/2||Fx-e||2 + rho/2||x-n||2
	!  x
	! s.t. Cx = d
	!
	! Laurence Yang, UCSD
	!
	! 19 Aug 2018: first version
	! 11 Sep 2018: normalized objective so magnitude (relative to
	! 			   proximal term) stays invariant with problem dimension.
	! 09 Oct 2018: error stop --> stop for compatibility with PGI.
	!----------------------------------------------------

	use po_linsys,	only: 	opt_linsys, &
							mA11,  nA11,  A11nnz,  pA11,  jA11,  vA11, &
							mA12,  nA12,  pA12,  jA12,  vA12,  A12nnz, &
							mA21,  nA21,  pA21,  jA21,  vA21,  A21nnz, &
							mA22,  nA22,  pA22,  jA22,  vA22,  A22nnz, &
							mA, nA,  pA,  jA,  vA,  Annz, &
							pAc0, iAc0, &
							b1, b2, stat,	&
							linsys_mF, linsys_nF,				&
							linsys_pF, linsys_jF, linsys_vF,	&
							linsys_e, 							&
							linsys_pFF, linsys_jFF, linsys_vFF,	&
							linsys_pFt, linsys_jFt, linsys_vFt, &
							diagFF

	use sparsekit
	use sparseops

	implicit none

	integer, parameter 		:: ip = 4, sp = 4, dp = 8

	public make_PO_matrices

	private

contains

	subroutine make_PO_matrices(mF, nF, nnzF, pF, jF, vF,	&
								mC, nC, nnzC, pC, jC, vC,	& 
								e, d, Drho, n, shift, nCond)

		!----------------------------------------------------
		! Make matrix blocks and rhs for PO given F, C, d for
		!
		! min ||Fx-e||2 + rho/2||x-n||2
		!  x
		! s.t. Cx = d
        !
		! INPUTS
		! 	mF, nF		row and col dimensions of F
		! 	mC, nC		row and col dimensions of C
		!	pF, jF, vF	F matrix in CSR	
		!	pC, jC, vC	C matrix in CSR	
		! 	e			real array
		! 	d			real array
		! 	Drho		diagonal(nF) matrix of rho
		! 	n			real array, proximal point
        !
		! OUTPUTS
		! 	None.
		!
		! 	Modifies A11, A12, A21, A22, b1, b2
		! 	in the po_linsys module.
		!----------------------------------------------------

		integer(ip),	intent(in)					:: mF 
		integer(ip),	intent(in)					:: nF 
		integer(ip),	intent(in)					:: mC 
		integer(ip),	intent(in)					:: nC 
		integer(ip),	intent(in)					:: nnzF 
		integer(ip),	intent(in)					:: nnzC 

		integer(ip),	intent(in)					:: pF(mF+1)
		integer(ip),	intent(in)					:: jF(nnzF)
		real(dp),		intent(in)					:: vF(nnzF)

		! transpose(F) 
		integer(ip)									:: mFt, nFt, nnzFt
		integer(ip)									:: pFt(nF+1)
		integer(ip)									:: jFt(nnzF)
		real(dp)									:: vFt(nnzF)

		integer(ip),	intent(in)					:: pC(mC+1)
		integer(ip),	intent(in)					:: jC(nnzC)
		real(dp),		intent(in)					:: vC(nnzC)

		! Fully stacked A: CSC
		integer(ip),   	allocatable                 :: pAc(:)	! 1-based indexing
		integer(ip),   	allocatable                 :: iAc(:)	! 1-based indexing
		real(dp),	   	allocatable              	:: ac(:)

		real(dp),		intent(in)					:: e(mF)
		real(dp),		intent(in)					:: d(mC)
		real(dp),		intent(in)					:: Drho(nF)
		real(dp),		intent(in)					:: n(nF)

		real(dp),		intent(in),		optional	:: shift
		real(dp)                                    :: ishift

		real(dp),		intent(in),		optional	:: nCond
		real(dp)									:: inCond

		! Internal
		integer(ip)									:: ierr
		! integer(ip)									:: iw(2*nF)
		integer(ip)									:: iw(nF)
		! integer(ip),	allocatable					:: iw(:)
		integer(ip)									:: ndegr(nF)

		! CSR rep of diagonal Drho
		integer(ip)									:: pDrho(nF+1)	! nrow+1
		integer(ip)									:: jDrho(nF) 
		real(dp)									:: vDrho(nF) 

		real(dp)									:: Fe(nF)
!		real(dp)									:: RhoN(nF)

		integer(ip)									:: mFF, nFF
		integer(ip) 	          					:: pFF(nF+1)
		integer(ip),	allocatable 	          	:: jFF(:) ! nzmax !nnzF
		real(dp),		allocatable 	          	:: vFF(:) ! nzmax !nnzF
		integer(ip)									:: nzmax
		integer(ip)									:: nnzFF

		! Temporary arrays before finalizing FF
		integer(ip),	allocatable 	          	:: pFF0(:) ! nF+1
		integer(ip),	allocatable 	          	:: jFF0(:) ! nzmax !nnzF
		real(dp),		allocatable 	          	:: vFF0(:) ! nzmax !nnzF

		integer(ip)									:: ndiag
		integer(ip)									:: idiag
		integer(ip)									:: ioff(nF)

		integer(ip)									:: mFFc, nFFc, nnzFFc

		integer(ip)									:: i,j,k,l


		!----------------------------------------------------
		! Create matrix blocks:
		! [2*F'*F + Drho,  C';  * [x; = [2*F'*e + Drho*n;
		!  C          ,  0 ]     z]                d]
		!----------------------------------------------------

		if (present(shift)) then
			ishift = shift
		else
			ishift = 0.0_dp
		endif

		if (present(nCond)) then
			inCond = nCond
		else
			inCond = 1.0_dp
		endif

		!----------------------------------------------------
		! A11
		!----------------------------------------------------
        ! F'*F
		! AMUB performs the matrix product C = A * B.

		nzmax = nnzF

		! Initialize
		pFt = 0
		jFt = 0
		vFt = 0.0

		call sp_transpose(mF, nF , nnzF, pF, jF, vF,	&
						  mFt,nFt,nnzFt,pFt,jFt,vFt)

		!----------------------------------------------------
		! Need this for later: adaptive rho, precond
		!----------------------------------------------------
		linsys_pFt = pFt
		linsys_jFt = jFt
		linsys_vFt = vFt


		! F'*F
		! 18 Oct 2018 [LY]: nzmax should be enough but if not,
		! catch with ierr > 0 by amub()
        ! Look ahead to get right size
		call amubdg(nF, nF, nF, jFt, pFt, jF, pF, ndegr, nnzFF, iw)
		nzmax = nnzFF

		if (allocated(vFF0)) deallocate(vFF0)
		if (allocated(jFF0)) deallocate(jFF0)
		allocate(vFF0(nzmax))
		allocate(jFF0(nzmax))

		if (allocated(pFF0)) deallocate(pFF0)
		allocate(pFF0(nF+1))

		mFF = nF
		nFF = nF

		! Initialize solution arrays
		pFF0 = 0
		jFF0 = 0
		vFF0 = 0.0

		call amub (	nF, nF, 1, 				&
					vFt, jFt, pFt, 			&	! F'
					vF, jF, pF, 			&	! F
					vFF0, jFF0, pFF0, nzmax, 	&
				   	iw, ierr )

		if (ierr > 0) then
			write(*,*) 'ERROR: amub() for F*F failed with ierr=',ierr
			write(*,*) 'This means nnz of A*B exceeds nzmax provided: nzmax=',nzmax
			stop
		endif
	
		vFF0 = 2.0_dp/inCond * vFF0

		! If shift != 0, then might change nnz
		if (ishift == 0) then
			pFF = pFF0
			jFF = jFF0
			vFF = vFF0
		else
            call addDiag(mFF, pFF0, jFF0, vFF0, [(-ishift, i=1, mFF)], &
						 nnzFF, pFF, jFF, vFF)
		endif

		!----------------------------------------------------
        ! A11: 2*F'F + Drho - s*I
		!----------------------------------------------------
		! APLDIA adds a diagonal matrix to a general sparse matrix:  B = A + Diag.
		mA11 = nF

		! Get nnz correctly
        ! Initialize arrays
		pDrho = 0
		jDrho = 0
		vDrho = 0.0
		iw = 0

        call coocsr(nF, nF, Drho, [(i,i=1,nF)], [(i,i=1,nF)], vDrho, jDrho, pDrho)

		! SIZES
		nA11 = nF

		! Allocation
		if (allocated(pA11)) deallocate(pA11)
		allocate(pA11(mA11+1))

		!----------------------------------------------------
		! Need this for later: residual balancing, etc.
		!----------------------------------------------------
		linsys_vFF = vFF
		linsys_jFF = jFF
		linsys_pFF = pFF

		! Need diagonal part of vFF before adding rho
		if (allocated(diagFF)) deallocate(diagFF)
		allocate(diagFF(nF))
		diagFF = 0.0_dp

		do i=1, nF
			do k=pFF(i), pFF(i+1)-1
				j = jFF(k)
				if (i==j) then
					diagFF(i) = vFF(k)
				endif
			enddo
		enddo

		!----------------------------------------------------
		! Add diagonal rho
		!----------------------------------------------------
		call addDiag(mFF, pFF, jFF, vFF, Drho, A11nnz, pA11, jA11, vA11)

		!----------------------------------------------------
		! A12: C'
		!----------------------------------------------------

		!****************************************************
		! ASSERTION
		if (nC /= nF) then
			write(*,*) 'nC =',nC, 'nF=',nF
			stop  'nC /= nF'
		endif
		!****************************************************

		mA12 = nC
        A12nnz = nnzC

		if (allocated(vA12)) deallocate(vA12)
		if (allocated(pA12)) deallocate(pA12)
		if (allocated(jA12)) deallocate(jA12)
		allocate(vA12(A12nnz))
		allocate(jA12(A12nnz))
		allocate(pA12(mA12+1))

		! Initialize arrays
		pA12 = 0
		jA12 = 0
		vA12 = 0.0

		call sp_transpose(mC, nC, nnzC, pC, jC, vC,	&
						  mA12, nA12, A12nnz, pA12, jA12, vA12)

		!----------------------------------------------------
		! A21
		!----------------------------------------------------
		! Since Fortran 2003, allocatable arrays automatically
		! allocated upon assignment, including reallocating to
		! correct shape if mismatch.
		! But maybe don't count on this...
        mA21 = mC
		nA21 = nC
		A21nnz = nnzC

		if (allocated(vA21)) deallocate(vA21)
		if (allocated(pA21)) deallocate(pA21)
		if (allocated(jA21)) deallocate(jA21)
		allocate(vA21(A21nnz))
		allocate(jA21(A21nnz))
		allocate(pA21(mA21+1))

		pA21 = pC
		jA21 = jC
		vA21 = vC

		!----------------------------------------------------
		! A22
		!----------------------------------------------------
		if (ishift == 0) then
			A22nnz = 0
		else
			!------------------------------------------------
			! -shift*I
			!------------------------------------------------
			if (allocated(vA22)) deallocate(vA22)
			if (allocated(pA22)) deallocate(pA22)
			if (allocated(jA22)) deallocate(jA22)
			allocate(vA22(A22nnz))
			allocate(jA22(A22nnz))
			allocate(pA22(mA22+1))
            A22nnz 	= mC
			mA22	= mC
			nA22	= mC
			pA22	= [(i, i=1,mA22+1)]
			jA22	= [(i, i=1,A22nnz)]
			vA22	= [(-ishift, i=1,A22nnz)]

            !************************************************
			! DEBUG
            !write(*,'(*(A,I15))') '#nF =',nF, ', nC =',nC
            !write(*,'(*(A,I15))') '#mC =',mC
            !write(*,'(*(A,I15))') '#nA11 =',nA11, ', nA12 =',nA12
            !write(*,'(*(A,I15))') '#nA21 =',nA21, ', nA22 =',nA22
            !write(*,'(*(A,I15))') '#nA1 = nA11 + nA12 =', nA11+nA12
            !write(*,'(*(A,I15))') '#nA2 = nA21 + nA22 =', nA21+nA22
			if ((nA11+nA12) /= (nA21+nA22)) then
				stop 'columns of A1 and A2 not matching'
			endif
            !************************************************

		endif

		!====================================================
		! RHS
		!====================================================
		! AMUX multiplies a CSR matrix A times a vector.
		! Initialize
		Fe = 0.0

		call amux(nF, e, Fe, vFt, jFt, pFt)

		! Drho*n
		! AMUXD multiplies a DIA matrix times a vector.
		ndiag = nF
		idiag = nF
		ioff = 0
        ! Initialize array

		!****************************************************
		! 25 Aug 2018 [LY]: Fortran 2003 has automatic allocation
		! so don't need to allocate(b1(mA1)) first. Just assigning
		! it works.
		! 11 Sep 2018 [LY]: normalize by dimension so magnitude
		! of objective relative to proximal penalty invariant
		! to data 
		! 21 Sep 2018 [LY]: changed update model to
		! b1 = b1fix + Drho*N
		!****************************************************
		b1 = 2.0_dp/inCond*Fe
		b2 = d

		!----------------------------------------------------
		! 03 Sep 2018 [LY]: save F and e for objval calc
		!----------------------------------------------------
        linsys_mF = mF
        linsys_nF = nF
		linsys_vF = vF
		linsys_jF = jF
		linsys_pF = pF
		linsys_e  = e

		!----------------------------------------------------
		! Deallocate
		!----------------------------------------------------
		deallocate(jFF0)
		deallocate(vFF0)
		deallocate(pFF0)

	end subroutine

end module
