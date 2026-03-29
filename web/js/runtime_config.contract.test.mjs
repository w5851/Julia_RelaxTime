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

function run() {
    test_default_localhost_profile();
    test_profile_mapping();
    test_global_profile_override();
    console.log('runtime_config.contract.test.mjs: PASS');
}

run();
