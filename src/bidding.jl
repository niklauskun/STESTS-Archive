# generate energy storage bid using deep learning

using Flux, JLD2, Random

function loadbidmodels(
    model_filenames;
    input_size = 60,
    dense_size = 60,
    output_size = 5,
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

function assign_models_to_storages(params, models, num_storages)
    storage_to_model_map = []
    storage_to_index_map = []

    model_keys = collect(keys(models))  # Get the keys of the models dictionary to choose from

    for storage_id in 1:num_storages
        if params.EStrategic[storage_id] == 1
            # Select a model for this storage in a round-robin fashion
            selected_model_index = ((storage_id - 1) % length(model_keys)) + 1
            selected_model_key = model_keys[selected_model_index]
            push!(
                storage_to_model_map,
                (storage_id, models[selected_model_key]),
            )
            push!(storage_to_index_map, (storage_id, selected_model_index)) # Track storage_id and model index
        end
    end
    df = DataFrame(StorageID = first.(storage_to_index_map), SelectedModelIndex = last.(storage_to_index_map))
    CSV.write("storage_to_index_map.csv", df)

    return storage_to_model_map
end

function update_battery_storage!(
    params::STESTS.ModelParams,
    ratio::Float64,
    output_folder::String,
)
    # Assume additional fields in Params for simplicity of explanation
    new_storagemap = copy(params.storagemap)
    new_Eeta = []
    new_EPC = []
    new_EPD = []
    new_ESOC = []
    new_ESOCini = []
    new_EStrategic = []

    # Track indices to remove
    remove_indices = []

    for region in 1:size(params.storagemap, 2)
        # Find battery storage indices in this region
        battery_indices = findall(
            i -> params.Eeta[i] == 0.9 && params.storagemap[i, region] == 1,
            1:length(params.Eeta),
        )

        if isempty(battery_indices)
            continue
        end

        # Aggregate values
        total_EPC = sum(params.EPC[battery_indices])
        total_EPD = sum(params.EPD[battery_indices])
        total_ESOC = sum(params.ESOC[battery_indices])
        total_ESOCini = sum(params.ESOCini[battery_indices])

        # Prepare data for two new storages
        append!(new_Eeta, [0.9, 0.9])
        append!(new_EPC, [total_EPC * ratio, total_EPC * (1 - ratio)])
        append!(new_EPD, [total_EPD * ratio, total_EPD * (1 - ratio)])
        append!(new_ESOC, [total_ESOC * ratio, total_ESOC * (1 - ratio)])
        append!(
            new_ESOCini,
            [total_ESOCini * ratio, total_ESOCini * (1 - ratio)],
        )
        append!(new_EStrategic, [1, 0])

        # Update storagemap for new entries
        # Assuming new_row_strategic is a 1xN vector
        new_row_strategic = zeros(Int64, size(new_storagemap, 2))
        new_row_strategic[region] = 1

        # To add two identical rows to the matrix, ensure they are treated as separate rows
        new_storagemap = vcat(
            new_storagemap,
            reshape(new_row_strategic, 1, :),
            reshape(new_row_strategic, 1, :),
        )

        append!(remove_indices, battery_indices)
    end

    # Save the new entries to csv
    df = DataFrame(
        Eeta = new_Eeta,
        EPC = new_EPC,
        EPD = new_EPD,
        ESOC = new_ESOC,
        ESOCini = new_ESOCini,
        EStrategic = new_EStrategic,  # Assuming you also have this array
    )

    CSV.write(joinpath(output_folder * "/Strategic", "ADDED_ES.csv"), df)

    # Now, remove the aggregated entries and update params
    # Append new entries
    params.Eeta = vcat(params.Eeta, new_Eeta)
    params.EPC = vcat(params.EPC, new_EPC)
    params.EPD = vcat(params.EPD, new_EPD)
    params.ESOC = vcat(params.ESOC, new_ESOC)
    params.ESOCini = vcat(params.ESOCini, new_ESOCini)
    params.EStrategic = vcat(params.EStrategic, new_EStrategic)

    # Create a mask that is true for indices that should be kept
    remove_indices = sort(unique(remove_indices))
    mask = trues(length(params.Eeta))
    mask[remove_indices] .= false
    params.Eeta = params.Eeta[mask]
    params.EPC = params.EPC[mask]
    params.EPD = params.EPD[mask]
    params.ESOC = params.ESOC[mask]
    params.ESOCini = params.ESOCini[mask]
    params.EStrategic = params.EStrategic[mask]
    params.storagemap = new_storagemap[mask, :]
    return params
end

function CalcValueNoUnc(d, c, eta, vi, iC, iD)
    """
    ####### Completely translated to Julia, gave exact numbers to Python #######

    Title: Calculate Risk-Neutral value function using deterministic price
    Inputs:
        d - price right now
        c - marginal discharge cost
        eta - efficiency
        vi - input value function for the next time period, which equals to
        v_t(e) where e is sampled from 0 to 1 at the granularity e
    Outputs:
        vo - value function for the current time period sampled at ed
    """
    # add a large number of upper and lower v, where the first point is
    # v_t(0-) = +infty, and the second point is v_t(0), the second largest one is
    # v_t(1), and the largest one is v_t(1+) = -infty
    lNum = fill(1e5)
    v_foo = [lNum; vi; -lNum]

    # calculate soc after charge vC = v_t(e+P*eta)
    vC = v_foo[iC]

    # calculate soc after discharge vC = v_t(e-P/eta)
    vD = v_foo[iD]

    # calculate CDF and PDF
    FtEC = (vi .* eta .> d)
    FtCC = (vC .* eta .> d)
    FtED = ((vi ./ eta .+ c) .* ((vi ./ eta .+ c) .> 0) .> d)
    FtDD = ((vD ./ eta .+ c) .* ((vD ./ eta .+ c) .> 0) .> d)

    # calculate terms
    Term1 = vC .* FtCC
    Term2 = d .* (vC .* eta .<= d) .* (vi .* eta .> d) ./ eta
    Term3 = vi .* (FtED - FtEC)
    Term4 =
        d .* (((vi ./ eta .+ c) .* ((vi ./ eta .+ c) .> 0)) .<= d) .*
        (((vD ./ eta .+ c) .* ((vD ./ eta .+ c) .> 0)) .> d) .* eta
    Term5 = -c .* eta .* (FtDD - FtED)
    Term6 = vD .* (1 .- FtDD)

    # output new value function sampled at ed
    vo = Term1 + Term2 + Term3 + Term4 + Term5 + Term6
    return vo
end

function generate_value_function(
    T,
    P,
    eta,
    num_segment,
    RTP;
    c = 10.0,
    ed = 0.001,
    Ne = 1001,
    ef = 0.5,
)
    """
    ####### Completely translated to Julia, gave exact numbers to Python #######
    Generate value function v and downsampled value function vAvg
    """

    # start_time = time()

    # Set final SoC level
    vEnd = zeros(Ne)
    vEnd[1:floor(Int64, ef * Ne)] .= 1e2 # Use 100 as the penalty for final discharge level

    # Define the risk-neutral value function and populate the final column.
    # v[1, 1] is the marginal value of 0% SoC at the beginning of day 1, v[Ne, T] is the marginal value of 100% SoC at the beginning of the last operating day
    v = zeros(Ne, T + 1) # initialize the value function series
    v[:, end] = vEnd  # v.shape == (1001, 210241)

    # Process indices: discretize vt by modeling it as a vector v_{t,j} in which each element is associated with equally spaced SoC samples
    es = collect(0:ed:1)

    # the number of samples is J = 1 + E/ed
    Ne = length(es)

    # Calculate SoC after charge vC = v_t(e+P*eta)
    eC = es .+ P * eta  # [0.0375, 0.0385, 0.0395, ..., 1.0355, 1.0365, 1.0375]
    iC = max.(0, min.(Ne + 1, ceil.(Int, eC / ed))) .+ 1

    # Calculate SoC after discharge vC = v_t(e-P/eta)
    eD = es .- P / eta
    iD = max.(0, min.(Ne + 1, floor.(Int, eD / ed))) .+ 1

    # Populate value function
    for t in T:-1:1 # start from the last day and move backwards
        vi = v[:, t+1] # input value function of next time stamp
        vo = CalcValueNoUnc(RTP[t], c, eta, vi, iC, iD)
        v[:, t] = vo # record the result
    end

    # end_time = time()
    # println("Time:", end_time - start_time)

    # Downsample
    vAvg = reshape(
        reduce(
            vcat,
            [
                mean(
                    v[
                        (i-1)*Int(1 / ed / num_segment)+1:i*Int(
                            1 / ed / num_segment,
                        ),
                        :,
                    ],
                    dims = 1,
                ) for i in 1:num_segment
            ],
        ),
        num_segment,
        :,
    )

    return vAvg
end