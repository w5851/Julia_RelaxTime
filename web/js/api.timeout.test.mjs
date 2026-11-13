import assert from 'node:assert/strict';
import { API, validate_scan_request } from './api.js';

const originalFetch = globalThis.fetch;
const originalSetTimeout = globalThis.setTimeout;
const originalClearTimeout = globalThis.clearTimeout;

try {
    // Simulate headers arriving immediately, followed by a stalled JSON body.
    globalThis.setTimeout = (fn) => originalSetTimeout(fn, 5);
    let requests = 0;
    globalThis.fetch = async (_url, { signal }) => {
        requests += 1;
        return {
            ok: true,
            json: () => new Promise((_resolve, reject) => {
                signal.addEventListener('abort', () => reject(new DOMException('aborted', 'AbortError')));
            }),
        };
    };
    await assert.rejects(API.getJobStatus('stalled-body'), { name: 'AbortError' });
    assert.equal(requests, 2);

    requests = 0;
    await assert.rejects(API.computeScattering({}), { name: 'AbortError' });
    assert.equal(requests, 1, 'scattering POST is not automatically retried');

    globalThis.fetch = async () => ({ ok: true, json: async () => ({ success: true, data: { value: 42 } }) });
    assert.deepEqual(await API.computeScattering({}), { value: 42 });

    globalThis.fetch = async () => ({ ok: true, json: async () => { throw new SyntaxError('invalid JSON'); } });
    await assert.rejects(API.getJobStatus('invalid-json'), { name: 'SyntaxError' });

    const pending = new Set();
    globalThis.setTimeout = (fn, ms) => {
        const id = originalSetTimeout(fn, ms);
        pending.add(id);
        return id;
    };
    globalThis.clearTimeout = (id) => { pending.delete(id); originalClearTimeout(id); };
    globalThis.fetch = async () => { throw new TypeError('offline'); };
    assert.equal(await API.checkHealth(), false);
    assert.equal(pending.size, 0);

    for (const invalid of [null, true, '', ' ', []]) {
        const checked = validate_scan_request({
            kind: 'tmu', params: { T_values: [150], mu_values: [invalid], xi: 0 },
        });
        assert.equal(checked.valid, false, `invalid numeric input: ${JSON.stringify(invalid)}`);
    }
} finally {
    globalThis.fetch = originalFetch;
    globalThis.setTimeout = originalSetTimeout;
    globalThis.clearTimeout = originalClearTimeout;
}
console.log('api.timeout.test.mjs: PASS');
