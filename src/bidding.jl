# generate energy storage bid using deep learning

using Flux, JLD2

function loadbidmodel(
    model_path;
    input_size = 60,
    dense_size = 60,
    output_size = 50,
    activation_fn = relu,
)
    model = Chain(
        Dense(input_size, dense_size, activation_fn),   # activation function inside layer
        Dense(dense_size, dense_size, activation_fn),   # activation function inside layer
        Dense(dense_size, output_size),
    )

    Flux.loadmodel!(model, model_state)
    return model
end