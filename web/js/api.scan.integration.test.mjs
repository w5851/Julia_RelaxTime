import assert from 'node:assert/strict';
import test from 'node:test';
import { API } from './api.js';

async function createRecordedJob(payload, acceptedIds) {
    const created = await API.createScanJob(payload);
    if (typeof created.job_id === 'string') {
        acceptedIds.push(created.job_id);
    }
    assert.equal(created.status, 'accepted');
    assert.equal(typeof created.job_id, 'string');
    return created;
}

function buildValidPayload() {
    return {
        kind: 'tmu',
        params: {
            T_values: [150.0, 160.0, 170.0],
            mu_values: [0.0, 100.0, 200.0],
            xi_values: [0.0, 0.1],
            max_retries: 0,
            p_num: 16,
            t_num: 8,
        },
    };
}

async function testResultAvailability(acceptedIds) {
    const created = await createRecordedJob(buildValidPayload(), acceptedIds);
    let result;
    try {
        result = await API.getJobResult(created.job_id);
    } catch (error) {
        assert.equal(error.code, 'JOB_NOT_SUCCEEDED');
        assert.equal(error.status, 409);
        return;
    }
    // A fast server may finish before the result request arrives.
    assert.equal(result.status, 'ok');
    assert.equal(result.job_status, 'succeeded');
    assert.equal(result.job_id, created.job_id);
}

async function testInvalidOutputPath(acceptedIds) {
    const badPayload = {
        kind: 'tmu',
        params: {
            T_values: [150.0],
            mu_values: [0.0],
            xi: 0.0,
            max_retries: 0,
            output_path: '../outside_outputs_forbidden.csv',
        },
    };

    await assert.rejects(() => createRecordedJob(badPayload, acceptedIds), {
        code: 'INVALID_REQUEST',
        status: 400,
    });
}

async function testQueueFull(acceptedIds) {
    const heavyPayload = {
        kind: 'tmu',
        params: {
            T_values: Array.from({ length: 24 }, (_, i) => 120 + i * 5),
            mu_values: Array.from({ length: 24 }, (_, i) => i * 30),
            xi_values: [0.0, 0.05, 0.1, 0.15],
            max_retries: 0,
            p_num: 24,
            t_num: 12,
        },
    };

    let queueFullHit = false;
    for (let i = 0; i < 60; i += 1) {
        try {
            await createRecordedJob(heavyPayload, acceptedIds);
        } catch (error) {
            if (error.code === 'QUEUE_FULL' && error.status === 429) {
                queueFullHit = true;
                break;
            }
            throw error;
        }
    }

    assert.equal(queueFullHit, true, 'queue-full case should be reproducible');
}

async function cancelCreatedJobs(acceptedIds) {
    const results = await Promise.allSettled(acceptedIds.map(async (jobId) => {
        try {
            await API.cancelJob(jobId);
        } catch (error) {
            if (error.status !== 409 || error.code !== 'JOB_NOT_CANCELLABLE') {
                throw error;
            }
        }
    }));
    const failures = results.filter((result) => result.status === 'rejected');
    assert.equal(failures.length, 0, `failed to cancel test jobs: ${failures.map((result) => result.reason.message).join('; ')}`);
}

test('live scan API', { skip: process.env.JRT_RUN_LIVE_API_TESTS !== '1' }, async (t) => {
    assert.ok(process.env.JRT_TEST_API_BASE_URL?.trim(), 'set JRT_TEST_API_BASE_URL to a dedicated test server');
    globalThis.__JRT_API_BASE_URL__ = process.env.JRT_TEST_API_BASE_URL;
    assert.equal(await API.checkHealth(), true, 'server health should be OK before integration tests');

    const acceptedIds = [];
    try {
        await t.test('invalid output path is rejected before registration', () => testInvalidOutputPath(acceptedIds));
        await t.test('result availability', () => testResultAvailability(acceptedIds));
        await t.test('queue capacity', { skip: process.env.JRT_TEST_QUEUE_FULL !== '1' }, () => testQueueFull(acceptedIds));
    } finally {
        await cancelCreatedJobs(acceptedIds);
    }
});
