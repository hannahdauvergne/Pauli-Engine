"""
	solve_eigs(H; nE, vin)
	solve_eigs(H)
	
	Diagonalise H. It uses ARPACK for sparse matrices and LinearAlgebra for dense matrices. For sparse matrices it computes the low part of the spectrum.

	Inputs:

	H:			Hamiltonian matrix to diagonalise.
	nE:			Number of smallest real eigenvalues to obtain (only for sparse matrix).
	v0:			Optional starting vector for diagonalisation. It is useful to use a previously calculated similar solution.

	Returns:

	eigsE:		Vector containing eigenvalues.
	eigsψ:		Matrix where each row is an eigenvector.

"""
function solve_eigs(H::SparseMatrixCSC; nE::Int=1, tol::Real=0, v0::AbstractVector{<:Number}=Float64[])
	# Extract basis' size
	D = size(H,1)
    # Solve eigenproblem for nE smallest real eigenvalues
	eigsE, eigsψ, nconv = eigs(H, nev=nE, ncv=2*min(D,max(20,2*nE+1)), which=:SR, maxiter=10*D, tol=tol, v0=v0)
	# Check convergence
	nconv >= nE || println("WARNING: Diagonalisation did not converge.")
	# Return
	return eigsE, eigsψ
end

function solve_eigs(H::Matrix{<:Number})
    # Solve eigenproblem
	F = eigen(H)
	eigsE = F.values
	eigsψ = F.vectors
	# Return
	return eigsE, eigsψ
end

