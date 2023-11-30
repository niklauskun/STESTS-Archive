# generate energy storage bid using deep learning

using Flux, JLD2, Random

function loadbidmodels(
    model_filenames;
    input_size = 60,
    dense_size = 60,
    output_size = 50,
    activation_fn = relu,
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

function assign_models_to_storages(models, num_storages)
    storage_to_model_map = []

    model_keys = collect(keys(models))  # Get the keys of the models dictionary to choose from

    for storage_id in 1:num_storages
        # Select a model for this storage in a round-robin fashion
        selected_model_index = ((storage_id - 1) % length(model_keys)) + 1
        selected_model_key = model_keys[selected_model_index]
        push!(storage_to_model_map, (storage_id, models[selected_model_key]))
    end

    return storage_to_model_map
end