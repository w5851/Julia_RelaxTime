function route_request(req::HTTP.Request, repo_root::String)
    path = String(HTTP.URIs.unescapeuri(req.target))

    if path == "/health"
        return HTTP.Response(200, ["Content-Type" => "text/plain"], "OK")
    elseif path == "/compute"
        return handle_compute(req)
    elseif path == "/api/modules" && req.method == "GET"
        return handle_modules_list()
    elseif path == "/api/modules/pnjl-gap/run"
        return handle_pnjl_single_point(req)
    else
        static_path = String(split(path, '?')[1])
        return serve_static_file(static_path, repo_root)
    end
end

function build_app(repo_root::String)
    return cors_middleware(req -> route_request(req, repo_root))
end
