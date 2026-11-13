# Phase-2 Split Deployment Runbook

## 1. Scope

- Target: frontend static assets + Julia backend API split deployment baseline.
- Profiles: `localhost`, `staging`, `remote`.
- Goal: keep localhost default path stable while allowing profile-based runtime switching.

## 2. Runtime Policy Matrix

| Profile | API Base URL (frontend default) | CORS Allow Origins (backend) | Max Running | Max Pending | Job Timeout (s) |
|---|---|---|---:|---:|---:|
| localhost | `http://localhost:8080` | `*` | 2 | 32 | 0 |
| staging | `https://staging.jrt.local` | `https://staging.jrt.local` | 2 | 64 | 180 |
| remote | `https://api.jrt.example.com` | `https://api.jrt.example.com` | 4 | 128 | 300 |

## 3. Startup

### 3.1 Localhost (default)

```bash
julia --project=. scripts/server/server_full.jl
```

- No extra env required.
- `ServerLauncher` fills missing env with localhost defaults.
- `PNJL_SCAN_MAX_RUNNING` and `PNJL_SCAN_MAX_PENDING` override positive queue limits when the app is built. Restart the server to apply changed limits; active jobs prevent reconfiguration. The startup banner reports effective values.

### 3.2 Staging profile

```powershell
$env:JRT_DEPLOY_PROFILE = 'staging'
julia --project=. scripts/server/server_full.jl
```

### 3.3 Remote profile

```powershell
$env:JRT_DEPLOY_PROFILE = 'remote'
julia --project=. scripts/server/server_full.jl
```

## 4. Frontend Switching

Frontend profile priority:

1. `globalThis.__JRT_DEPLOY_PROFILE__`
2. URL query `?deploy_profile=...`
3. `localStorage.JRT_DEPLOY_PROFILE`
4. without an explicit profile, HTTP(S) pages use their own origin and port; other contexts fall back to `http://localhost:8080`

Manual override (stronger than profile URL mapping):

- `globalThis.__JRT_API_BASE_URL__`
- URL query `?api_base_url=...`
- `localStorage.JRT_API_BASE_URL`

For split deployment, set `JRT_CORS_ALLOW_ORIGINS` to the frontend origin(s), separated by commas. The middleware reads it when the app is built. `*` permits all origins; an explicit list rejects requests with an unlisted `Origin` with HTTP 403. Requests without `Origin` remain available. Allowed responses include `Vary: Origin`; preflight supports `Content-Type` and `Idempotency-Key`.

## 5. Health Checks

- Backend health endpoint:

```bash
curl http://localhost:8080/health
```

- Modules endpoint:

```bash
curl http://localhost:8080/api/modules
```

- Split-config contract tests:

```bash
julia --project=. -e 'include("tests/integration/relaxtime/test_frontend_backend_config_contract.jl")'
node web/js/runtime_config.contract.test.mjs
```

## 6. Rollback

If staging/remote profile causes issue:

1. reset profile to localhost
   - `$env:JRT_DEPLOY_PROFILE = 'localhost'`
   - clear any explicit backend CORS override: `Remove-Item Env:JRT_CORS_ALLOW_ORIGINS -ErrorAction SilentlyContinue`
2. clear frontend profile overrides
   - remove query params `deploy_profile` / `api_base_url`
   - clear `localStorage` keys `JRT_DEPLOY_PROFILE`, `JRT_API_BASE_URL`
3. restart backend with default command
   - `julia --project=. scripts/server/server_full.jl`

## 7. Troubleshooting

- Symptom: frontend still hits old API base URL
  - Check `localStorage` override first.
  - Check URL query override.
- Symptom: queue saturates quickly
  - Confirm current `JRT_DEPLOY_PROFILE` and corresponding `PNJL_SCAN_MAX_PENDING`.
- Symptom: timeout appears in status
  - Check `PNJL_SCAN_JOB_TIMEOUT_SECONDS` and long task characteristics.
