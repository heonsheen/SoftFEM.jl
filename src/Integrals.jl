using LinearAlgebra

# Line integral of a quadratic function X^T A X from P0 to P1
function quadratic_line_int(
    A::Matrix{Float64}, # 2x2 matrix
    P0::Vector{Float64}, P1::Vector{Float64})
    x0 = P0[1]; y0 = P0[2]
    x1 = P1[1]; y1 = P1[2]
    a11 = A[1,1]; a12 = A[1,2]; a21 = A[2,1]; a22 = A[2,2]
    norm(P1 - P0) * (a11/3.0 * (x1^2 + x1*x0 + x0^2) +
                     a22/3.0 * (y1^2 + y1*y0 + y0^2) +
                     (a12 + a21) * (x0*y1/6.0 + y0*x1/6.0 + x1*y1/3.0 + x0*y0/3.0))
end

# Line integral of a linear function B X from P0 to P1
function linear_line_int(
    B::Matrix{Float64}, # 1x2 matrix
    P0::Vector{Float64}, P1::Vector{Float64})
    x0 = P0[1]; y0 = P0[2]
    x1 = P1[1]; y1 = P1[2]
    B11 = B[1,1]; B12 = B[1,2]
    norm(P1 - P0) * (B11*x0 + B12*y0 + 0.5 * B11*(x1-x0) + 0.5 * B12*(y1-y0))
end

# Line integral of a constant c from P0 to P1
function const_line_int(
    c::Float64,
    P0::Vector{Float64}, P1::Vector{Float64})
    norm(P1 - P0) * c
end