# function getLastOneIndex(matrix::Array{Int, 2}, vector::Array{Int, 1}, vector2::Array{Int, 1})
#     m, n = size(matrix)
#     result = Vector{Int}(undef, m)
    
#     for i in 1:m
#         idx = findlast(==(1), matrix[i, :])
#         result[i] = idx === nothing ? max(0, vector[m] - 24) : max(0, vector2[m] - (24 + 1 - idx))
#     end
    
#     return result
# end

# # Sample usage:
# matrix = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
#           0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
#           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
# vector = [0;0;15]
# vector2 = [8;168;168]

# vector = getLastOneIndex(matrix, vector, vector2)
# println(vector)

# using JuMP

# import Gurobi
# model = Model(Gurobi.Optimizer)

# set_silent(model)

# @variable(model, x >= 0)

# @constraint(model, c1, x >= 2)
# @constraint(model, c2, x <= 1)

# optimize!(model)
# compute_conflict!(model)


# iis_model, _ = copy_conflict(model)
# print(iis_model)

original_vector = [i for i in 1:25]  # Replace this with your 25x1 vector
m = 10  # Replace this with your value for m

# Create the new m x 25 matrix
new_matrix = repeat(original_vector', m, 1)

# Expand each row to make it a m x 300 matrix
expanded_matrix = repeat(new_matrix, inner = (1, 12))
 

println(expanded_matrix)