for type in subtypes(SourceOfCoefficients)
    print(type, " ")
    for i in 2:10
        operator = nothing
        try
            operator = derivative_operator(type(), derivative_order=1, accuracy_order=i,
                                xmin=xmin, xmax=xmax, N=num_nodes)
            print(i, " ")
        catch
            try
                operator = derivative_operator(type(:central), derivative_order=1, accuracy_order=i,
                                xmin=xmin, xmax=xmax, N=num_nodes)
                print(i, "(central) ")
            catch
            end
        end
    end

    println()
end