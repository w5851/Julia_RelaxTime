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

function run() {
    test_validate_scan_request();
    test_build_api_url();
    test_normalize_api_error();
    console.log('api.scan.test.mjs: PASS');
}

run();
