reciprocity(matrix) = (matrix + transpose(matrix)) / 2

function reciprocity!(matrix::AbstractMatrix)
    n = checksquare(matrix)
    @inbounds for column in 1:n, row in 1:(column - 1)

        value = (matrix[row, column] + matrix[column, row]) / 2
        matrix[row, column] = value
        matrix[column, row] = value
    end
    return matrix
end

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
    first_value = nominal(matrix[1, 1])
    matrix_norm_squared = zero(real(abs2(first_value)))
    offdiagonal_norm_squared = zero(matrix_norm_squared)
    @inbounds for column in 1:n, row in 1:n

        squared = abs2(nominal(matrix[row, column]))
        matrix_norm_squared += squared
        row == column || (offdiagonal_norm_squared += squared)
    end
    matrix_norm = sqrt(matrix_norm_squared)
    offdiagonal_norm = sqrt(offdiagonal_norm_squared)
    return offdiagonal_norm / max(
        matrix_norm,
        eps(float(real(matrix_norm)))
    )
end

function suppress!(matrix::AbstractMatrix, tolerance::Real)
    ratio = offdiagonal_ratio(matrix)
    ratio <= tolerance || return ratio
    n = checksquare(matrix)
    for column in 1:n, row in 1:n

        row == column || (matrix[row, column] = zero(eltype(matrix)))
    end
    return ratio
end
