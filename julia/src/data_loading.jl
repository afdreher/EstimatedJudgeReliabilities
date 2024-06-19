using Arrow                 # Compressed loading (lz4)
using CSV                   # CSV loading
using DataFrames            # Data frames (table storage)
using DataFramesMeta
using OrderedCollections     # OrderedDict

function load_preferences(preference_file::String, compressed::Bool = false)
    if lowercase(last(preference_file, 3)) == "csv"
        # Try to read the preference file. This should be updated to take the winner / loser format, but that's not for now
        preferences = CSV.read(preference_file, DataFrame;
                         header=["attribute", "hit_id", "user_id", "font_A_name", "font_B_name", "user_choice"])
    elseif lowercase(last(preference_file, 3)) == "lz4"
        preferences = DataFrame(Arrow.Table(preference_file))
        if !("winner_id" in names(preferences))
            DataFramesMeta.rename!(preferences, [:winner] .=>  [:winner_id])
        end
        preferences.user_choice = map(x -> x == 1 ? "more" : "less", preferences.winner_id)

        @eachrow! preferences begin
            @newcol :winner::Vector{String}
            @newcol :loser::Vector{String}

            if :winner_id == 1
                :winner = string(:font_A_name)
                :loser = string(:font_B_name)
            else
                :loser = string(:font_A_name)
                :winner = string(:font_B_name)
            end
        end

        #print(preferences)
    end
    if compressed && (ncol(preferences) == 6)
        # Merge to something better...
        df = DataFrame(map(x -> OrderedDict(:attribute => x.attribute,
                                            :user_id => x.user_id,
                                            :winner => x.user_choice == "more" ? x.font_A_name : x. font_B_name, 
                                            :loser => x.user_choice == "more" ? x.font_B_name : x. font_A_name), 
                        eachrow(preferences)))
        preferences = df
    end

    return preferences
end