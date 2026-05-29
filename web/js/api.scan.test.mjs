import assert from 'node:assert/strict';
import {
    API,
    build_api_url,
    normalize_api_error,
    summarize_scan_request_metrics,
    validate_pnjl_gap_request,
    validate_scan_request,
    validate_script_task_request,
    validate_transport_point_request,
} from './api.js';

function test_validate_scan_request() {
    const ok = validate_scan_request({
        kind: 'tmu',
        params: {
            T_values: [150.0, 160.0],
            mu_values: [0.0, 100.0],
            xi: 0.0,
        },
    });
    assert.equal(ok.valid, true);

    const bad_kind = validate_scan_request({ kind: 'bad', params: { T_values: [150.0], mu_values: [0.0], xi: 0.0 } });
    assert.equal(bad_kind.valid, false);

    const bad_xi_strategy = validate_scan_request({
        kind: 'tmu',
        params: {
            T_values: [150.0],
            mu_values: [0.0],
            xi: 0.0,
            xi_values: [0.0, 0.1],
        },
    });
    assert.equal(bad_xi_strategy.valid, false);

    const bad_trho = validate_scan_request({ kind: 'trho', params: { T_values: [150.0], mu_values: [0.0], xi: 0.0 } });
    assert.equal(bad_trho.valid, false);

    const bad_xi_range = validate_scan_request({
        kind: 'tmu',
        params: {
            T_values: [150.0],
            mu_values: [0.0],
            xi: 1.5,
        },
    });
    assert.equal(bad_xi_range.valid, false);

    const bad_temperature_range = validate_scan_request({
        kind: 'tmu',
        params: {
            T_values: [-10.0],
            mu_values: [0.0],
            xi: 0.0,
        },
    });
    assert.equal(bad_temperature_range.valid, false);

    const bad_empty_numeric = validate_scan_request({
        kind: 'tmu',
        params: {
            T_values: [''],
            mu_values: [0.0],
            xi: 0.0,
        },
    });
    assert.equal(bad_empty_numeric.valid, false);
}

function test_validate_point_requests() {
    const gapOk = validate_pnjl_gap_request({
        params: {
            T_mev: 150.0,
            mu_mev: 0.0,
            xi: 0.0,
        },
    });
    assert.equal(gapOk.valid, true);

    const gapFixedRhoOk = validate_pnjl_gap_request({
        params: {
            T_mev: 150.0,
            rho_target: 0.1,
            xi: 0.0,
        },
    });
    assert.equal(gapFixedRhoOk.valid, true);

    const gapBad = validate_pnjl_gap_request({
        params: {
            T_mev: -1.0,
            mu_mev: 0.0,
            xi: 0.0,
        },
    });
    assert.equal(gapBad.valid, false);

    const transportOk = validate_transport_point_request({
        params: {
            T_mev: 150.0,
            mu_mev: 0.0,
            xi: 0.0,
            tau: 1.0,
            transport: {
                p_nodes: 8,
                p_max: 3.5,
            },
        },
    });
    assert.equal(transportOk.valid, true);

    const transportBad = validate_transport_point_request({
        params: {
            T_mev: 150.0,
            mu_mev: 0.0,
            xi: 0.0,
            tau: -1.0,
        },
    });
    assert.equal(transportBad.valid, false);
}

function test_validate_script_task_request() {
    const task = {
        id: 'run-gap-transport-scan',
        default_preset: 'smoke',
        presets: {
            smoke: { heavy: false },
            canonical: { heavy: true },
            custom: { heavy: true },
        },
    };

    const smoke = validate_script_task_request({
        task_id: 'run-gap-transport-scan',
        preset: 'smoke',
    }, task);
    assert.equal(smoke.valid, true);

    const heavyWithoutConfirm = validate_script_task_request({
        task_id: 'run-gap-transport-scan',
        preset: 'canonical',
    }, task);
    assert.equal(heavyWithoutConfirm.valid, false);

    const customOk = validate_script_task_request({
        task_id: 'run-gap-transport-scan',
        preset: 'custom',
        confirm_heavy: true,
        custom_args: ['--output', 'data/outputs/results/demo.csv'],
    }, task);
    assert.equal(customOk.valid, true);

    const customMissingArgs = validate_script_task_request({
        task_id: 'run-gap-transport-scan',
        preset: 'custom',
        confirm_heavy: true,
        custom_args: [],
    }, task);
    assert.equal(customMissingArgs.valid, false);
}

function test_build_api_url() {
    globalThis.__JRT_API_BASE_URL__ = 'http://127.0.0.1:9000';
    assert.equal(build_api_url('/health'), 'http://127.0.0.1:9000/health');
    delete globalThis.__JRT_API_BASE_URL__;
}

function test_normalize_api_error() {
    const err = normalize_api_error({
        message: 'queue is full',
        status: 429,
        payload: {
            code: 'QUEUE_FULL',
            diagnostics: { max_pending: 32 },
        },
    });
    assert.equal(err.code, 'QUEUE_FULL');
    assert.equal(err.status, 429);
    assert.deepEqual(err.diagnostics, { max_pending: 32 });
}

async function test_format_error_for_timeout_and_offline() {
    globalThis.__JRT_API_BASE_URL__ = 'http://127.0.0.1:9000';
    const originalFetch = globalThis.fetch;

    let abortAttempts = 0;
    globalThis.fetch = async () => {
        abortAttempts += 1;
        const err = new Error('The operation was aborted');
        err.name = 'AbortError';
        throw err;
    };

    let timeoutMessage = '';
    try {
        await API.getJobStatus('job-timeout');
    } catch (error) {
        timeoutMessage = API.formatError(error);
    }
    assert.equal(abortAttempts, 2);
    assert.equal(timeoutMessage, '请求超时，请检查服务器是否正常运行');

    let offlineAttempts = 0;
    globalThis.fetch = async () => {
        offlineAttempts += 1;
        throw new TypeError('fetch failed');
    };

    let offlineMessage = '';
    try {
        await API.getJobStatus('job-offline');
    } catch (error) {
        offlineMessage = API.formatError(error);
    }
    assert.equal(offlineAttempts, 2);
    assert.equal(offlineMessage, '无法连接到服务器，请确保Julia服务器正在运行 (julia server.jl)');

    globalThis.fetch = originalFetch;
    delete globalThis.__JRT_API_BASE_URL__;
}

async function test_scan_request_400_ratio_metrics() {
    const originalFetch = globalThis.fetch;
    let fetchCalls = 0;

    globalThis.fetch = async () => {
        fetchCalls += 1;
        return {
            ok: true,
            status: 202,
            json: async () => ({ status: 'accepted', job_id: `job-${fetchCalls}` }),
        };
    };

    const requests = [
        {
            kind: 'tmu',
            params: {
                T_values: [150.0],
                mu_values: [0.0],
                xi: 0.0,
            },
        },
        {
            kind: 'tmu',
            params: {
                T_values: [-1.0],
                mu_values: [0.0],
                xi: 0.0,
            },
        },
        {
            kind: 'tmu',
            params: {
                T_values: [160.0],
                mu_values: [20.0],
                xi: 0.2,
            },
        },
        {
            kind: 'tmu',
            params: {
                T_values: [''],
                mu_values: [0.0],
                xi: 0.0,
            },
        },
        {
            kind: 'tmu',
            params: {
                T_values: [170.0],
                mu_values: [50.0],
                xi: 1.2,
            },
        },
    ];

    const metrics = await summarize_scan_request_metrics(requests, async (payload) => {
        await API.createScanJob(payload);
    });

    assert.equal(metrics.total_requests, 5);
    assert.equal(metrics.client_blocked, 3);
    assert.equal(metrics.backend_requests_after, 2);
    assert.equal(metrics.backend_400_baseline, 3);
    assert.equal(metrics.backend_400_after, 0);
    assert.equal(metrics.backend_400_ratio_baseline, 0.6);
    assert.equal(metrics.backend_400_ratio_after, 0);
    assert.equal(metrics.backend_400_reduction, 3);

    assert.equal(fetchCalls, 2);

    globalThis.fetch = originalFetch;
}

async function test_cancel_job_request_contract() {
    const originalFetch = globalThis.fetch;
    globalThis.__JRT_API_BASE_URL__ = 'http://127.0.0.1:9000';

    let capturedUrl = '';
    let capturedMethod = '';
    globalThis.fetch = async (url, options = {}) => {
        capturedUrl = String(url);
        capturedMethod = String(options.method || 'GET');
        return {
            ok: true,
            status: 200,
            json: async () => ({
                status: 'ok',
                job_id: 'job-cancel-1',
                job_status: 'cancelled',
            }),
        };
    };

    const payload = await API.cancelJob('job-cancel-1');
    assert.equal(capturedMethod, 'POST');
    assert.equal(capturedUrl, 'http://127.0.0.1:9000/api/modules/pnjl-scan/jobs/job-cancel-1/cancel');
    assert.equal(payload.job_status, 'cancelled');

    globalThis.fetch = originalFetch;
    delete globalThis.__JRT_API_BASE_URL__;
}

async function test_point_service_request_contracts() {
    const originalFetch = globalThis.fetch;
    globalThis.__JRT_API_BASE_URL__ = 'http://127.0.0.1:9000';

    const calls = [];
    globalThis.fetch = async (url, options = {}) => {
        calls.push({
            url: String(url),
            method: String(options.method || 'GET'),
            body: JSON.parse(String(options.body || '{}')),
        });
        return {
            ok: true,
            status: 200,
            json: async () => ({ status: 'ok', result: {} }),
        };
    };

    await API.runPnjlGap({
        params: {
            T_mev: 150.0,
            mu_mev: 0.0,
            xi: 0.0,
        },
    });
    await API.runTransportPoint({
        params: {
            T_mev: 150.0,
            mu_mev: 0.0,
            xi: 0.0,
            tau: 1.0,
            p_num: 8,
            t_num: 4,
            transport: {
                p_nodes: 8,
                p_max: 3.5,
            },
        },
    });

    assert.equal(calls[0].method, 'POST');
    assert.equal(calls[0].url, 'http://127.0.0.1:9000/api/modules/pnjl-gap/run');
    assert.equal(calls[0].body.params.T_mev, 150.0);
    assert.equal(calls[1].method, 'POST');
    assert.equal(calls[1].url, 'http://127.0.0.1:9000/api/modules/transport-point/run');
    assert.equal(calls[1].body.params.transport.p_nodes, 8);

    assert.equal(API.buildPnjlGapUrl(), 'http://127.0.0.1:9000/api/modules/pnjl-gap/run');
    assert.equal(API.buildTransportPointUrl(), 'http://127.0.0.1:9000/api/modules/transport-point/run');

    globalThis.fetch = originalFetch;
    delete globalThis.__JRT_API_BASE_URL__;
}

async function test_script_task_request_contracts() {
    const originalFetch = globalThis.fetch;
    globalThis.__JRT_API_BASE_URL__ = 'http://127.0.0.1:9000';

    const calls = [];
    globalThis.fetch = async (url, options = {}) => {
        calls.push({
            url: String(url),
            method: String(options.method || 'GET'),
            body: options.body ? JSON.parse(String(options.body)) : null,
        });
        const last = calls[calls.length - 1];
        if (last.url.endsWith('/api/modules/script-tasks')) {
            return {
                ok: true,
                status: 200,
                json: async () => ({
                    status: 'ok',
                    tasks: [{
                        id: 'run-gap-transport-scan',
                        default_preset: 'smoke',
                        presets: { smoke: { heavy: false } },
                    }],
                }),
            };
        }
        if (last.url.endsWith('/api/modules/script-tasks/jobs')) {
            return {
                ok: true,
                status: 202,
                json: async () => ({ status: 'accepted', job_id: 'script-job-1', kind: 'script-task' }),
            };
        }
        if (last.url.endsWith('/cancel')) {
            return {
                ok: true,
                status: 200,
                json: async () => ({ status: 'ok', job_id: 'script-job-1', job_status: 'cancelled' }),
            };
        }
        if (last.url.endsWith('/result')) {
            return {
                ok: true,
                status: 200,
                json: async () => ({ status: 'ok', job_id: 'script-job-1', job_status: 'succeeded', result: {} }),
            };
        }
        return {
            ok: true,
            status: 200,
            json: async () => ({ status: 'ok', job_id: 'script-job-1', job_status: 'running' }),
        };
    };

    const catalog = await API.getScriptTaskCatalog();
    await API.createScriptTaskJob({ task_id: 'run-gap-transport-scan', preset: 'smoke' }, catalog.tasks[0]);
    await API.getScriptTaskJobStatus('script-job-1');
    await API.getScriptTaskJobResult('script-job-1');
    await API.cancelScriptTaskJob('script-job-1');

    assert.equal(calls[0].method, 'GET');
    assert.equal(calls[0].url, 'http://127.0.0.1:9000/api/modules/script-tasks');
    assert.equal(calls[1].method, 'POST');
    assert.equal(calls[1].url, 'http://127.0.0.1:9000/api/modules/script-tasks/jobs');
    assert.equal(calls[1].body.task_id, 'run-gap-transport-scan');
    assert.equal(calls[2].url, 'http://127.0.0.1:9000/api/modules/script-tasks/jobs/script-job-1');
    assert.equal(calls[3].url, 'http://127.0.0.1:9000/api/modules/script-tasks/jobs/script-job-1/result');
    assert.equal(calls[4].url, 'http://127.0.0.1:9000/api/modules/script-tasks/jobs/script-job-1/cancel');
    assert.equal(API.buildScriptTaskCatalogUrl(), 'http://127.0.0.1:9000/api/modules/script-tasks');
    assert.equal(API.buildScriptTaskCreateUrl(), 'http://127.0.0.1:9000/api/modules/script-tasks/jobs');

    globalThis.fetch = originalFetch;
    delete globalThis.__JRT_API_BASE_URL__;
}

async function run() {
    test_validate_scan_request();
    test_validate_point_requests();
    test_validate_script_task_request();
    test_build_api_url();
    test_normalize_api_error();
    await test_format_error_for_timeout_and_offline();
    await test_scan_request_400_ratio_metrics();
    await test_cancel_job_request_contract();
    await test_point_service_request_contracts();
    await test_script_task_request_contracts();
    console.log('api.scan.test.mjs: PASS');
}

run().catch((error) => {
    console.error('api.scan.test.mjs: FAIL', error);
    process.exit(1);
});
