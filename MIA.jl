using LinearAlgebra
using SciPy

function Prod_Diff(m, N)
    P = zeros(2 * N + 1, 2 * N + 1)
    P[1, 1] = 1.0

    for i in 2:2 * N + 1
        P[i-1, 2] = (-1)^(i) * m[i - 1]
    end

    for j in 4:2 * N + 2
        for i in 2:2 * N + 4 - j
            P[i-1, j - 1] = (P[1, j - 2] * P[i, j - 3]) - (P[1, j - 3] * P[i, j - 2])
        end
    end

    zeta = zeros(2 * N)
    for i in 3:2 * N + 1
        if P[1, i - 1] * P[1, i - 2] > 0.0
            zeta[i - 1] = P[1, i] / (P[1, i - 1] * P[1, i - 2])
        else
            zeta[i - 1] = 0.0
        end
    end

    a = zeros(N)
    for i in 2:N+1
        a[i - 1] = zeta[2 * (i-1) ] + zeta[2 * (i-1) - 1]
    end
    

    b = zeros(N - 1)
    for i in 2:N 
        b[i - 1] = zeta[2 * (i-1)+1] * zeta[2 * (i-1)]
    end

    nodes, vectors = SciPy.linalg.eigh_tridiagonal(a, -sqrt.(b))
    #nodes, vectors = eigen(SymTridiagonal(a, b))
    weights = zeros(N)
    for i in 1:N
        weights[i] = vectors[1, i]^2
    end
    weights *= m[1]
    # println("ffp: ", P)
    # println("ffzeta: ", zeta)
    # println("ffa: ", a)
    # println("ffb: ", b)
    # println("ffnodes: ", nodes)
    # println("ffweights: ", weights)

    return weights, nodes
end

function Realisable(m, N)
    matrix0 = zeros(N, N)
    matrix1 = zeros(N, N)

    for i in 1:N
        matrix0[i, :] .= m[i:N + i - 1]
        matrix1[i, :] .= m[i + 1:N + i]
    end

    det0 = det(matrix0)
    det1 = det(matrix1)
    if det0 < 0 || det1 < 0
        #println("Moment set is unrealisable")
        return false
    else
        # println("Moment set is realisable")
        return true
    end
end

# Rest of your code...
