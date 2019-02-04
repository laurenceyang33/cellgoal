module po_constants

	!--------------------------------------------------------
	! Proximal operator constants
	!--------------------------------------------------------

	implicit none

	integer,		parameter			:: ip=4, dp=8
	integer(ip),	parameter			:: PO_IND_BOUNDS 		= 1
	integer(ip),	parameter			:: PO_IND_BILINEAR 		= 2
	integer(ip),	parameter			:: PO_IND_L1 			= 3
	integer(ip),	parameter			:: PO_IND_LINSYS 		= 4
	integer(ip),	parameter			:: PO_IND_CARD 			= 5
	integer(ip),	parameter			:: PO_IND_BOUNDS_L1 	= 13
	integer(ip),	parameter			:: PO_IND_BOUNDS_CARD 	= 15
	integer(ip),	parameter			:: PO_IND_LINSYS_MA57 	= 54

	!--------------------------------------------------------
	! Multi-condition
	!--------------------------------------------------------
	integer(ip),	parameter			:: PO_IND_BILINEAR_MULT	= 102
	integer(ip),	parameter			:: PO_IND_LINSYS_MULT 	= 104

	!--------------------------------------------------------
	! fine-grained
	!--------------------------------------------------------
	integer(ip),	parameter			:: PO_IND_LINSYS_FG 	= 1004

	!--------------------------------------------------------
	! Sets of PO indices
	!--------------------------------------------------------
	integer(ip),	parameter			:: LINSYS_PO_INDS(4) = [	&
								PO_IND_LINSYS, PO_IND_LINSYS_MA57,	&
								PO_IND_LINSYS_MULT, PO_IND_LINSYS_FG]

	private ip, dp

	public

end module
