function cors_middleware(handler)
    return function (req::HTTP.Request)
        if req.method == "OPTIONS"
            return HTTP.Response(
                200,
                [
                    "Access-Control-Allow-Origin" => "*",
                    "Access-Control-Allow-Methods" => "POST, GET, OPTIONS",
                    "Access-Control-Allow-Headers" => "Content-Type",
                    "Access-Control-Max-Age" => "86400",
                ],
            )
        end

        resp = handler(req)
        HTTP.setheader(resp, "Access-Control-Allow-Origin" => "*")
        HTTP.setheader(resp, "Access-Control-Allow-Methods" => "POST, GET, OPTIONS")
        HTTP.setheader(resp, "Access-Control-Allow-Headers" => "Content-Type")
        return resp
    end
end

function serve_static_file(path::String, repo_root::String)
    if contains(path, "..") || contains(path, "\\")
        return HTTP.Response(403, "Forbidden")
    end

    clean_path = startswith(path, "/") ? path[2:end] : path
    file_path = if path == "/" || path == ""
        joinpath(repo_root, "web", "index.html")
    else
        joinpath(repo_root, "web", clean_path)
    end

    if !isfile(file_path)
        return HTTP.Response(404, "Not Found: $path")
    end

    content_type = if endswith(file_path, ".html")
        "text/html; charset=utf-8"
    elseif endswith(file_path, ".css")
        "text/css; charset=utf-8"
    elseif endswith(file_path, ".js")
        "application/javascript; charset=utf-8"
    elseif endswith(file_path, ".json")
        "application/json; charset=utf-8"
    else
        "application/octet-stream"
    end

    try
        content = read(file_path)
        return HTTP.Response(200, ["Content-Type" => content_type], body = content)
    catch e
        @error "Error serving file" file = file_path exception = e
        return HTTP.Response(500, "Internal Server Error")
    end
end
