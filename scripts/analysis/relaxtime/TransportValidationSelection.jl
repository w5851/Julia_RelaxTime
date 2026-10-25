module TransportValidationSelection

export requested_validation_rows

"""Select explicit validation points without accepting missing or repeated IDs."""
function requested_validation_rows(source,point_ids)
    ids=String.(collect(point_ids))
    isempty(ids) && throw(ArgumentError("point_ids must not be empty"))
    length(unique(ids))==length(ids) || throw(ArgumentError("point_ids must be unique"))
    selected=filter(row->String(row.point_id) in ids,source)
    found=String.(getproperty.(selected,:point_id))
    length(found)==length(ids) && Set(found)==Set(ids) ||
        throw(ArgumentError("requested validation points are missing or duplicated in the input"))
    return selected
end

end
