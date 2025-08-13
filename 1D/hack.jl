function create_periodic_fd_matrix(stencil::Vector{T}, N::Int) where T<:Number
    # --- Input Validation ---
    if iseven(length(stencil))
        error("Stencil length must be odd to have a clear center.")
    end

    # --- Initialization ---
    # For efficient sparse matrix construction, we first create vectors for the
    # row indices (I), column indices (J), and the non-zero values (V).
    I = Int[]
    J = Int[]
    V = T[]

    # The center of the stencil corresponds to the main diagonal (offset = 0).
    center_idx = (length(stencil) + 1) ÷ 2

    # --- Matrix Population ---
    # Iterate through each row of the conceptual dense matrix.
    for row in 1:N
        # For each row, apply the stencil.
        for (stencil_idx, stencil_val) in enumerate(stencil)
            # We only need to store non-zero values.
            if stencil_val != 0
                # Calculate the column offset from the main diagonal.
                # A negative offset is a sub-diagonal, positive is a super-diagonal.
                offset = stencil_idx - center_idx

                # Calculate the column index for the current row and offset.
                # mod1(x, N) performs modulo arithmetic in the range 1 to N,
                # which elegantly handles the periodic (wrap-around) boundaries.
                col = mod1(row + offset, N)

                # Store the coordinates and value.
                push!(I, row)
                push!(J, col)
                push!(V, stencil_val)
            end
        end
    end

    # --- Sparse Matrix Assembly ---
    # Create the sparse matrix from the I, J, and V vectors.
    return sparse(I, J, V, N, N)
end