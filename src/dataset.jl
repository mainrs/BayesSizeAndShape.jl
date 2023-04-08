function dataset(dataset_name::AbstractString)
    
    basename = joinpath(@__DIR__, "..", "data")
 
    rdaname = joinpath(basename, string(dataset_name, ".rda"))
    #println(rdaname)
    if isfile(rdaname)
        
        dataset = load(rdaname)[dataset_name].data
        

        return (x = dataset[1], no = dataset[2], time = dataset[3])
    end

    #csvname = joinpath(basename, string(dataset_name, ".csv.gz"))
    #if isfile(csvname)
    #    return open(csvname,"r") do io
    #        uncompressed = IOBuffer(read(GzipDecompressorStream(io)))
    #        DataFrame(CSV.File(uncompressed, delim=',', quotechar='\"', missingstring="NA",
    #                  types=get(Dataset_typedetect_rows, (package_name, dataset_name), nothing)) )
    #    end
    #end
    #error("Unable to locate dataset file $rdaname or $csvname")
end

"""
    dataset_desciption(dataset_name::AbstractString)

It gives a description of the elements of the dataset (dataset_name) - the datasets are taken form the R package shape
"""
function dataset_desciption(dataset_name::AbstractString)

        if dataset_name == "rats"
            print("
            Description:

                Rat skulls data, from X rays. 8 landmarks in 2 dimensions, 18 individuals observed at 7, 14, 21, 30, 40, 60, 90, 150 days

                
            Format:    

                x: An array of landmark configurations 144 x 2 x 2
                no: Individual rat number (note rats 3, 13, 20 missing due to incomplete data)
                time: observed time in days
                
            ")
        else

            println("dataset_name not available")

        end

end

"""
    dataset_names()

It return the names of the available datasets - the datasets are taken form the R package shape
"""
function dataset_names()

    return ["rats"]

end