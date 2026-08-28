reciprocity(matrix) = (matrix + transpose(matrix)) / 2
reciprocity!(matrix) = (matrix .= (matrix .+ transpose(matrix)) ./ 2; matrix)

function ideal_transposition!(matrix::AbstractMatrix)
    n = checksquare(matrix)
    coefficients = similar(diag(matrix))
    @inbounds for offset in 0:(n - 1)
        total = zero(eltype(matrix))
        for row in 1:n
            total += matrix[row, 1 + mod(row - 1 + offset, n)]
        end
        coefficients[offset + 1] = total / n
    end
    @inbounds for row in 1:n, column in 1:n

        matrix[row, column] = coefficients[mod1(column - row + 1, n)]
    end
    return matrix
end

function offdiagonal_ratio(matrix::AbstractMatrix)
    n = checksquare(matrix)
    values = nominal(matrix)
    matrix_norm = sqrt(sum(abs2, values))
    offdiagonal_norm = sqrt(sum(
        abs2(values[row, column])
        for column in 1:n for row in 1:n if row != column;
        init = zero(real(matrix_norm))
    ))
    return offdiagonal_norm / max(
        matrix_norm,
        eps(float(real(matrix_norm)))
    )
end

function suppress_offdiagonal_residue!(matrix::AbstractMatrix, tolerance::Real)
    ratio = offdiagonal_ratio(matrix)
    ratio <= tolerance || return ratio
    n = checksquare(matrix)
    for column in 1:n, row in 1:n
        row == column || (matrix[row, column] = zero(eltype(matrix)))
    end
    return ratio
end
