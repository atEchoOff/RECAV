rows = rowvals(Q_skew)
vals = nonzeros(Q_skew)

function matrix_iterate_1(Q_skew, Q_skew_nzranges, rows, vals)
    accum = 0
    index = 0
    i = 0
    j = 0
    nij = 0
    skew_internal = 0
    for j in axes(Q_skew, 2)
        index = 0
        for skew_internal in nzrange(Q_skew, j)
            index += 1
            i = rows[skew_internal]
            nij = vals[skew_internal]

            accum += nij + index + i + j # make sure we use everything
        end
    end

    return accum
end

Q_skew_nzranges = UnitRange{Int64}[]

for j in axes(Q_skew, 2)
    push!(Q_skew_nzranges, nzrange(Q_skew, j))
end

function matrix_iterate_2(Q_skew, Q_skew_nzranges, rows, vals)
    accum = 0
    index = 0
    for j in axes(Q_skew, 2)
        index = 0
        for skew_internal in Q_skew_nzranges[j]
            index += 1
            i = rows[skew_internal]
            nij = vals[skew_internal]

            accum += nij + index + i + j # make sure we use everything
        end
    end

    return accum
end

function matrix_iterate_3(Q_skew, Q_skew_nzranges)
    accum = 0
    index = 0
    i = 0
    j = 0
    nij = 0
    skew_internal = 0
    for j in axes(Q_skew, 2)
        index = 0
        for i in Q_skew_nzranges[j]
            index += 1
            nij = Q_skew[i, j]

            accum += nij + index + i + j # make sure we use everything
        end
    end

    return accum
end

matrix_iterate_1(Q_skew, Q_skew_nzranges, rows, vals)
matrix_iterate_2(Q_skew, Q_skew_nzranges, rows, vals)
matrix_iterate_3(Q_skew, Q_skew_nzranges)

# @benchmark matrix_iterate_1(Q_skew, Q_skew_nzranges, rows, vals)
@benchmark matrix_iterate_2(Q_skew, Q_skew_nzranges, rows, vals)
# @benchmark matrix_iterate_3(Q_skew)