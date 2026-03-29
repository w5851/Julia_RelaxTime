import assert from 'node:assert/strict';
import { API } from './api.js';

async function sleep(ms) {
    return new Promise((resolve) => setTimeout(resolve, ms));
}

async function waitForTerminal(jobId, timeoutMs = 120000) {
    const start = Date.now();
    while (Date.now() - start < timeoutMs) {
        const status = await API.getJobStatus(jobId);
        if (status.job_status === 'succeeded' || status.job_status === 'failed') {
            return status;
        }
        await sleep(800);
    }
    throw new Error(`waitForTerminal timeout for ${jobId}`);
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

async function testResultNotReady() {
    const created = await API.createScanJob(buildValidPayload());
    assert.equal(created.status, 'accepted');

    let gotExpected = false;
    try {
        await API.getJobResult(created.job_id);
    } catch (error) {
        assert.equal(error.code, 'JOB_NOT_SUCCEEDED');
        assert.equal(error.status, 409);
        gotExpected = true;
    }

    assert.equal(gotExpected, true, 'result-not-ready should return JOB_NOT_SUCCEEDED');
    return created.job_id;
}

async function testJobFailed() {
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

    const created = await API.createScanJob(badPayload);
    assert.equal(created.status, 'accepted');

    const terminal = await waitForTerminal(created.job_id, 90000);
    assert.equal(terminal.job_status, 'failed');
    return created.job_id;
}

async function testQueueFull() {
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
    const acceptedIds = [];
    for (let i = 0; i < 60; i += 1) {
        try {
            const created = await API.createScanJob(heavyPayload);
            acceptedIds.push(created.job_id);
        } catch (error) {
            if (error.code === 'QUEUE_FULL' && error.status === 429) {
                queueFullHit = true;
                break;
            }
            throw error;
        }
    }

    assert.equal(queueFullHit, true, 'queue-full case should be reproducible');
    return acceptedIds.length;
}

async function run() {
    const healthy = await API.checkHealth();
    assert.equal(healthy, true, 'server health should be OK before integration tests');

    const resultNotReadyJobId = await testResultNotReady();
    const failedJobId = await testJobFailed();
    const acceptedBeforeQueueFull = await testQueueFull();

    console.log('api.scan.integration.test.mjs: PASS');
    console.log(`result-not-ready job: ${resultNotReadyJobId}`);
    console.log(`failed job: ${failedJobId}`);
    console.log(`accepted jobs before queue-full: ${acceptedBeforeQueueFull}`);
}

run().catch((error) => {
    console.error('api.scan.integration.test.mjs: FAIL', error);
    process.exit(1);
});
