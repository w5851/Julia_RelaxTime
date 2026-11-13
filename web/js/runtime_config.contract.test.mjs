import assert from 'node:assert/strict';
import {
    DEFAULT_API_BASE_URL,
    DEFAULT_DEPLOY_PROFILE,
    get_api_base_url,
    get_deploy_profile,
    resolve_profile_api_base_url,
} from './runtime_config.js';

function test_default_localhost_profile() {
    delete globalThis.__JRT_API_BASE_URL__;
    delete globalThis.__JRT_DEPLOY_PROFILE__;

    assert.equal(DEFAULT_DEPLOY_PROFILE, 'localhost');
    assert.equal(get_deploy_profile(), 'localhost');
    assert.equal(get_api_base_url(), DEFAULT_API_BASE_URL);
}

function test_profile_mapping() {
    assert.equal(resolve_profile_api_base_url('localhost'), 'http://localhost:8080');
    assert.equal(resolve_profile_api_base_url('staging'), 'https://staging.jrt.local');
    assert.equal(resolve_profile_api_base_url('remote'), 'https://api.jrt.example.com');
}

function test_global_profile_override() {
    globalThis.__JRT_DEPLOY_PROFILE__ = 'staging';
    delete globalThis.__JRT_API_BASE_URL__;
    assert.equal(get_deploy_profile(), 'staging');
    assert.equal(get_api_base_url(), 'https://staging.jrt.local');
    delete globalThis.__JRT_DEPLOY_PROFILE__;
}

function test_browser_origin_and_override_priority() {
    const values = new Map();
    globalThis.window = {
        location: { protocol: 'http:', origin: 'http://127.0.0.1:8081', search: '' },
        localStorage: { getItem: key => values.get(key) ?? null },
    };
    try {
        assert.equal(get_api_base_url(), 'http://127.0.0.1:8081');
        globalThis.__JRT_DEPLOY_PROFILE__ = 'remote';
        globalThis.__JRT_API_BASE_URL__ = 'http://127.0.0.1:9123/';
        assert.equal(get_api_base_url(), 'http://127.0.0.1:9123');
        delete globalThis.__JRT_API_BASE_URL__;
        window.location.search = '?api_base_url=http%3A%2F%2Flocalhost%3A9001';
        assert.equal(get_api_base_url(), 'http://localhost:9001');
        window.location.search = '';
        values.set('JRT_API_BASE_URL', 'http://localhost:9002');
        assert.equal(get_api_base_url(), 'http://localhost:9002');
        values.clear();
        delete globalThis.__JRT_DEPLOY_PROFILE__;
        window.localStorage.getItem = () => { throw new Error('storage disabled'); };
        assert.equal(get_api_base_url(), 'http://127.0.0.1:8081');
        window.location.protocol = 'file:';
        assert.equal(get_api_base_url(), DEFAULT_API_BASE_URL);
    } finally {
        delete globalThis.window;
        delete globalThis.__JRT_API_BASE_URL__;
        delete globalThis.__JRT_DEPLOY_PROFILE__;
    }
}

function run() {
    test_default_localhost_profile();
    test_profile_mapping();
    test_global_profile_override();
    test_browser_origin_and_override_priority();
    console.log('runtime_config.contract.test.mjs: PASS');
}

run();
