function route_request(req::HTTP.Request, repo_root::String)
    path = String(split(String(HTTP.URIs.unescapeuri(req.target)), '?')[1])
    job_prefix = "/api/modules/pnjl-scan/jobs/"

    if path == "/health"
        return HTTP.Response(200, ["Content-Type" => "text/plain"], "OK")
    elseif path == "/compute"
        return handle_compute(req)
    elseif path == "/api/modules" && req.method == "GET"
        return handle_modules_list()
    elseif path == "/api/modules/pnjl-gap/run"
        return handle_pnjl_single_point(req)
    elseif path == "/api/modules/pnjl-scan/jobs" && req.method == "POST"
        return handle_pnjl_scan_job_create(req)
    elseif startswith(path, job_prefix) && req.method == "GET"
        suffix = path[(length(job_prefix) + 1):end]
        if endswith(suffix, "/result")
            job_id = suffix[1:end-length("/result")]
            isempty(job_id) && return HTTP.Response(400, ["Content-Type" => "text/plain"], "Bad Request")
            return handle_pnjl_scan_job_result(job_id)
        else
            job_id = suffix
            isempty(job_id) && return HTTP.Response(400, ["Content-Type" => "text/plain"], "Bad Request")
            return handle_pnjl_scan_job_status(job_id)
        end
    elseif startswith(path, job_prefix) && req.method == "POST"
        suffix = path[(length(job_prefix) + 1):end]
        if endswith(suffix, "/cancel")
            job_id = suffix[1:end-length("/cancel")]
            isempty(job_id) && return HTTP.Response(400, ["Content-Type" => "text/plain"], "Bad Request")
            return handle_pnjl_scan_job_cancel(job_id)
        end
        return HTTP.Response(404, ["Content-Type" => "text/plain"], "Not Found")
    else
        return serve_static_file(path, repo_root)
    end
end

function build_app(repo_root::String)
    return cors_middleware(req -> route_request(req, repo_root))
end
