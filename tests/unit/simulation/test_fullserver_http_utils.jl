using Test
using HTTP

module FullServerHttpTest
using HTTP
include(joinpath(@__DIR__, "..", "..", "..", "src", "simulation", "fullserver", "http_utils.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "simulation", "fullserver", "routing.jl"))
end

@testset "FullServer HTTP behavior" begin
    @testset "configured CORS policy" begin
        calls = Ref(0)
        handler = req -> begin
            calls[] += 1
            HTTP.Response(200, ["Vary" => "Accept-Encoding"], "OK")
        end
        app = withenv("JRT_CORS_ALLOW_ORIGINS" => "https://frontend.test, https://other.test") do
            FullServerHttpTest.cors_middleware(handler)
        end
        allowed = app(HTTP.Request("GET", "/health", ["Origin" => "https://frontend.test"]))
        @test allowed.status == 200
        @test HTTP.header(allowed, "Access-Control-Allow-Origin") == "https://frontend.test"
        @test HTTP.header(allowed, "Vary") == "Accept-Encoding, Origin"

        preflight = app(HTTP.Request("OPTIONS", "/health", ["Origin" => "https://other.test"]))
        @test preflight.status == 200
        @test HTTP.header(preflight, "Access-Control-Allow-Origin") == "https://other.test"
        @test occursin("Idempotency-Key", HTTP.header(preflight, "Access-Control-Allow-Headers"))
        @test calls[] == 1

        for method in ("GET", "POST", "OPTIONS")
            blocked = app(HTTP.Request(method, "/health", ["Origin" => "https://unlisted.test"]))
            @test blocked.status == 403
            @test isempty(HTTP.header(blocked, "Access-Control-Allow-Origin"))
            @test HTTP.header(blocked, "Vary") == "Origin"
        end
        @test calls[] == 1

        no_origin = app(HTTP.Request("GET", "/health"))
        @test no_origin.status == 200
        @test isempty(HTTP.header(no_origin, "Access-Control-Allow-Origin"))
        @test calls[] == 2

        local_app = FullServerHttpTest.cors_middleware(handler; allow_origins="*")
        @test HTTP.header(local_app(HTTP.Request("GET", "/health")), "Access-Control-Allow-Origin") == "*"
    end

    @testset "static root and URL decoding" begin
        mktempdir() do root
            web = mkpath(joinpath(root, "web"))
            write(joinpath(web, "index.html"), "index")
            write(joinpath(web, "a b.js"), "export const value = 1;")
            write(joinpath(web, "version..js"), "export const version = 2;")
            outside = joinpath(root, "private.txt")
            write(outside, "outside web root")

            route(target) = FullServerHttpTest.route_request(HTTP.Request("GET", target), root)
            @test String(route("/").body) == "index"
            @test route("/a%20b.js?cache=1").status == 200
            @test HTTP.header(route("/a%20b.js"), "Content-Type") == "application/javascript; charset=utf-8"
            @test route("/version..js").status == 200
            @test route("/index.html%3Fignored").status != 200
            @test route("/missing.js").status == 404
            @test route("/%2e%2e/private.txt").status == 403
            @test route("/" * replace(outside, '\\' => '/')).status == 403
            @test route("//index.html").status == 403
            @test route("/C%3A/private.txt").status == 403
            @test route("/x%5C..%5Cprivate.txt").status == 403
            @test route("/index.html%00").status == 403

            if !Sys.iswindows()
                symlink(outside, joinpath(web, "linked.txt"))
                @test route("/linked.txt").status == 403
            end
        end
    end
end
