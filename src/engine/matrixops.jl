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
