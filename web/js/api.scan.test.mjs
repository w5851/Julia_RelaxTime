import assert from 'node:assert/strict';
import { build_api_url, normalize_api_error, validate_scan_request } from './api.js';

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

    const { API } = await import('./api.js');
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

async function run() {
    test_validate_scan_request();
    test_build_api_url();
    test_normalize_api_error();
    await test_format_error_for_timeout_and_offline();
    console.log('api.scan.test.mjs: PASS');
}

run().catch((error) => {
    console.error('api.scan.test.mjs: FAIL', error);
    process.exit(1);
});
