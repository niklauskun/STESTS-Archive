using STESTS, Flux, JLD2, Statistics

input_size = 60;
dense_size = 60;
output_size = 50;
activation_fn = relu;
model_filenames =
    ["models/WEST_1.jld2", "models/WEST_2.jld2", "models/WEST_3.jld2"]
num_segments = 5
segment_length = 50 ÷ num_segments

function load_models(
    model_filenames,
    input_size,
    dense_size,
    output_size,
    activation_fn,
)
    models = Dict{String,Chain}()

    for filename in model_filenames
        jldopen(filename, "r") do file
            # Define the model architecture and create a model
            model = Chain(
                Dense(input_size, dense_size, activation_fn),
                Dense(dense_size, dense_size, activation_fn),
                Dense(dense_size, output_size),
            )

            # Load the model state from the JLD2 file
            model_state = file["model_state"]
            Flux.loadmodel!(model, model_state)

            # Store the loaded model in the dictionary with filename as key
            return models[filename] = model
        end
    end

    return models
end

models = load_models(
    model_filenames,
    input_size,
    dense_size,
    output_size,
    activation_fn,
)

# models = Dict{String,Any}()  # Use a dictionary if you want to reference models by name

# for filename in model_filenames
#     jldopen(filename, "r") do file
#         # Assume that the model state is saved under the "model" key in the JLD2 file
#         model = Chain(
#             Dense(input_size, dense_size, activation_fn),   # activation function inside layer
#             Dense(dense_size, dense_size, activation_fn),   # activation function inside layer
#             Dense(dense_size, output_size),
#         )
#         model_state = file["model_state"]
#         Flux.loadmodel!(model, model_state)
#         return models[filename] = model  # If using a dictionary
#     end
# end

function use_models(models)
    # Here you can access and use the models
    for (name, model) in models
        # Do something with each model
        v = model(zeros(Float32, 60, 1))
        segment_averages = [
            mean(v[(i-1)*segment_length+1:i*segment_length]) for
            i in 1:num_segments
        ]
        cb = segment_averages .* 0.9
        db = segment_averages ./ 0.9 .+ 10
        println("Using model stored in $name")
        println("charge bid $cb")
        println("discharge bid $db")
        # For example, make a prediction
        # prediction = model(some_input_data)
    end
end

use_models(models)
