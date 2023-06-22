using JLD2
using DataFrames

"""
    read_jld2(filename::String)

Read system parameters and save into individual variables from a .jld2 file. 
Variables could be directly used in the unit commitment and economic dispatch construction.

# Example
```julia
read_jld2("./data/WECC240.jld2")
```
"""
function read_jld2(filename::String)
    @info "Reading data from $filename..."
    data = load(filename)
    @info "Done"
    return (
        UCL = data["UCL"],
        genmap = data["genmap"],
        GPmax = data["GPmax"],
        GPmin = data["GPmin"],
        GMC = data["GMC"],
        GSMC = data["GSMC"],
        GINCPmax = data["GINCPmax"],
        transmap = data["transmap"],
        TX = data["TX"],
        TFmax = data["TFmax"],
        GNLC = data["GNLC"],
        GRU = data["GRU"],
        GRD = data["GRD"],
        GSUC = data["GSUC"],
        GUT = data["GUT"],
        GDT = data["GDT"],
        GPIni = data["GPIni"],
        hydromap = data["hydromap"],
        HAvail = data["HAvail"],
        renewablemap = data["renewablemap"],
        RAvail = data["RAvail"],
        storagemap = data["storagemap"],
        EPC = data["EPC"],
        EPD = data["EPD"],
        Eeta = data["Eeta"],
        ESOC = data["ESOC"],
        ESOCini = data["ESOCini"],
        EDL = data["EDL"],
        EDHAvail = data["EDHAvail"],
        EDRAvail = data["EDRAvail"],
    )
end
