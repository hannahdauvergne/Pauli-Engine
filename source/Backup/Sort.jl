#====================================
    Sort.jl
    Routine to sort eigensolutions.
====================================#

"""
	sort_eigs!(eigsψ, eigsE, eigsψ0, δovlap)

	Sort the eigenvectors and eigenergies eigsψ and eigsE so that each solution has the largest overlap with eigsψ0.

	Inputs:

	eigψ:          	Target eigenvectors, given as a matrix where each column is an eigenvector.
	eigE:			Target eigenenergies.
	eigψ:			Eigenvectors to compare.
	δovlap:			Threshold to define a finite overlap.

"""

function sort_eigs!(
	eigsψ :: Matrix{<:Number}, 
	eigsE :: Vector{<:Number}, 
	eigsψ0 :: AbstractMatrix{<:Number}, 
	δovlap :: Real = 1E-2
	)
	
	# Check inputs
	size(eigsψ) == size(eigsψ0) || error("Target and reference eigenvectors matrices do not have the same size.")
	size(eigsψ,2) == length(eigsE) || error("Eigenvectors and eigenenergies arrays do not contain the same number of solutions.")
	# Number of eigenvectors
	nE = size(eigsψ,2)
	# Variables for overlaps
	lovlap0, lovlap,  = 0.0, 0.0
	novlap = 0
	# Loop over eigsψ0
	for n = 1:(nE-1)
		# Loop over eigsψ
		for n1 = n:nE
			lovlap0 = abs(dot(eigsψ0[:,n],eigsψ[:,n1]))
			if lovlap0 > lovlap
				lovlap = lovlap0
				novlap = n1
			end
		end
		# If overlap too small, assign larger eigensolution
		if lovlap < δovlap
			@views novlap = argmax(real(eigsE[n:nE]))+(n-1)
		end
		# Sort
		if n != novlap
			eigsψ[:,n], eigsψ[:,novlap] = eigsψ[:,novlap], eigsψ[:,n]
			eigsE[n], eigsE[novlap] = eigsE[novlap], eigsE[n]
		end
		# Reset lovlap
		lovlap = 0.0
	end
	# Return
	nothing
end
