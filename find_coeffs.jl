for type in subtypes(SourceOfCoefficients)
    operator = nothing
    try
        operator = derivative_operator(type(), derivative_order=1, accuracy_order=2,
                               xmin=xmin, xmax=xmax, N=num_nodes)
    catch
        try
            operator = derivative_operator(type(:central), derivative_order=1, accuracy_order=2,
                               xmin=xmin, xmax=xmax, N=num_nodes)
        catch
            continue
        end
    end

    try
        operator = derivative_operator(type(), derivative_order=1, accuracy_order=4,
                               xmin=xmin, xmax=xmax, N=num_nodes)
    catch
        try
            operator = derivative_operator(type(:central), derivative_order=1, accuracy_order=4,
                               xmin=xmin, xmax=xmax, N=num_nodes)
        catch
            continue
        end
    end

    println(type)
end