function cors_middleware(handler; allow_origins::AbstractString=get(ENV, "JRT_CORS_ALLOW_ORIGINS", "*"))
    origins = Set(strip.(split(allow_origins, ','; keepempty=false)))
    wildcard = "*" in origins
    return function (req::HTTP.Request)
        origin = HTTP.header(req, "Origin", "")
        if !wildcard && !isempty(origin) && !(origin in origins)
            return HTTP.Response(403, ["Vary" => "Origin"], "Origin not allowed")
        end

        resp = req.method == "OPTIONS" ? HTTP.Response(200) : handler(req)
        filter!(header -> lowercase(first(header)) != "access-control-allow-origin", resp.headers)
        if wildcard
            HTTP.setheader(resp, "Access-Control-Allow-Origin" => "*")
        else
            isempty(origin) || HTTP.setheader(resp, "Access-Control-Allow-Origin" => origin)
            vary = HTTP.header(resp, "Vary", "")
            if !any(value -> lowercase(strip(value)) in ("origin", "*"), split(vary, ','))
                HTTP.setheader(resp, "Vary" => (isempty(vary) ? "Origin" : "$vary, Origin"))
            end
        end
        HTTP.setheader(resp, "Access-Control-Allow-Methods" => "POST, GET, OPTIONS")
        HTTP.setheader(resp, "Access-Control-Allow-Headers" => "Content-Type, Idempotency-Key")
        req.method == "OPTIONS" && HTTP.setheader(resp, "Access-Control-Max-Age" => "86400")
        return resp
    end
end

function _within_static_root(path::String, root::String)
    relative = relpath(path, root)
    return !isabspath(relative) && !(".." in splitpath(relative))
end

function serve_static_file(path::String, repo_root::String)
    if any(char -> char in ('\\', ':', '\0'), path)
        return HTTP.Response(403, "Forbidden")
    end

    clean_path = startswith(path, "/") ? path[2:end] : path
    if startswith(clean_path, "/") || isabspath(clean_path) || ".." in split(clean_path, '/')
        return HTTP.Response(403, "Forbidden")
    end
    web_root = abspath(joinpath(repo_root, "web"))
    file_path = abspath(joinpath(web_root, isempty(clean_path) ? "index.html" : clean_path))
    _within_static_root(file_path, web_root) || return HTTP.Response(403, "Forbidden")

    if !isfile(file_path)
        return HTTP.Response(404, "Not Found: $path")
    end
    _within_static_root(realpath(file_path), realpath(web_root)) || return HTTP.Response(403, "Forbidden")

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
