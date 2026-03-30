import assert from 'node:assert/strict';
import {
    get_default_templates,
    load_task_history,
    load_task_templates,
    save_task_history,
    upsert_task_history_entry,
} from './task_store.js';

function createMemoryStorage() {
    const data = new Map();
    return {
        getItem(key) {
            return data.has(key) ? data.get(key) : null;
        },
        setItem(key, value) {
            data.set(key, String(value));
        },
        clear() {
            data.clear();
        },
    };
}

function setupStorage() {
    const storage = createMemoryStorage();
    globalThis.window = {
        localStorage: storage,
    };
    return storage;
}

function test_default_templates() {
    const templates = get_default_templates();
    assert.ok(Array.isArray(templates));
    assert.ok(templates.length >= 2);
    assert.equal(templates[0].id, 'tmu-smoke');
}

function test_history_roundtrip() {
    const storage = setupStorage();
    storage.clear();

    assert.deepEqual(load_task_history(), []);
    save_task_history([{ job_id: 'job-1', job_status: 'queued', kind: 'tmu' }]);
    const history = load_task_history();
    assert.equal(history.length, 1);
    assert.equal(history[0].job_id, 'job-1');
}

function test_upsert_history() {
    const storage = setupStorage();
    storage.clear();

    upsert_task_history_entry({ job_id: 'job-1', job_status: 'queued', kind: 'tmu' });
    upsert_task_history_entry({ job_id: 'job-2', job_status: 'running', kind: 'trho' });
    upsert_task_history_entry({ job_id: 'job-1', job_status: 'succeeded', kind: 'tmu' });

    const history = load_task_history();
    assert.equal(history.length, 2);
    assert.equal(history[0].job_id, 'job-1');
    assert.equal(history[0].job_status, 'succeeded');
}

function test_template_load_fallback() {
    const storage = setupStorage();
    storage.clear();

    const templates = load_task_templates();
    assert.ok(Array.isArray(templates));
    assert.ok(templates.length >= 1);
}

function run() {
    test_default_templates();
    test_history_roundtrip();
    test_upsert_history();
    test_template_load_fallback();
    console.log('task_store.test.mjs: PASS');
}

run();
